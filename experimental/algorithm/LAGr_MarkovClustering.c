#define LG_FREE_WORK                     \
    {                                    \
        GrB_free(&C);                    \
        GrB_free(&w);                    \
        GrB_free(&W);                    \
        GrB_free(&D);                    \
        GrB_free(&ones);                 \
        GrB_free(&C_temp_int64);          \
        GrB_free(&CC);                   \
        GrB_free(&C_scaled);             \
        GrB_free(&cs_row_vals);          \
        GrB_free(&vpc);                  \
        GrB_free(&col_seeds);            \
        GrB_free(&empty);                \
        LAGraph_Free((void *)&CI, NULL); \
        LAGraph_Free((void *)&CJ, NULL); \
        LAGraph_Free((void *)&CX, NULL); \
        GrB_free(&ROUND_HUN_FP64);       \
        GrB_free(&hash_unary_op);  \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(C_f); \
    }

#include "LG_internal.h"
#include <LAGraphX.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static void round_to_hundredths(void* z, const void* x) {
    double input = *((double*)x);
    double result = round(input * 100) / 100;
    *((double*)z) = result;
}

static void fp64_floor_to_int64(void* z, const void* x) {
    double input = *((double*)x);   // Convert input to double
    double floored = round(input);  // Apply the floor operation
    *((int64_t*)z) = (int64_t)floored; // Convert to uint64_t and store in output
}

int LAGr_MarkovClustering(
    // output:
    GrB_Matrix* C_f,
    // input
    int e,            // expansion coefficient
    int i,            // inflation coefficient
    double threshold, // threshold for pruning values
    int max_iter,     // maximum iterations
    LAGraph_Graph G,  // input graph
    char* msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    srand(time(NULL));

    // GxB_set(GxB_PRINT_1BASED, true);

    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // The newly computed cluster matrix at the end of each loop
    GrB_Matrix C_temp_int64 = NULL;
    GrB_Matrix C_scaled = NULL;
    GrB_Vector cs_row_vals = NULL; // Row reductions of C_scaled "c-scaled-values" 

    GrB_Vector w = NULL; // weight vector to normalize C matrix
    GrB_Matrix W = NULL;

    GrB_Matrix D = NULL;    // Diagonal workspace matrix
    GrB_Vector ones = NULL; // Vector of all 1's, used primarily to create identity

    GrB_Matrix CC = NULL;
    GrB_Vector* C_rows = NULL;

    GrB_Vector empty = NULL;

    GrB_UnaryOp ROUND_HUN_FP64 = NULL;
    GrB_UnaryOp hash_unary_op = NULL;

    GrB_Vector col_seeds = NULL;

    GrB_Vector vpc = NULL;

    GrB_Index* CI = NULL;
    GrB_Index* CJ = NULL;
    double* CX = NULL;

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG;

    GrB_Matrix A = G->A; // Adjacency matrix of G
    GrB_Index n;         // Dimension of A
    GRB_TRY(GrB_Matrix_nrows(&n, A));

    LG_ASSERT(C_f != NULL, GrB_NULL_POINTER);
    (*C_f) = NULL;
    LG_TRY(LAGraph_CheckGraph(G, msg));

    LG_ASSERT_MSG(G->out_degree != NULL, -106,
        "G->out_degree must be defined");

    //--------------------------------------------------------------------------
    // initializations
    //--------------------------------------------------------------------------


    GRB_TRY(GrB_Matrix_new(&C, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&w, GrB_FP64, n));
    GRB_TRY(GrB_Matrix_new(&W, GrB_FP64, n, n))
        GRB_TRY(GrB_Matrix_new(&D, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&vpc, GrB_INT16, n));
    GRB_TRY(GrB_Vector_new(&col_seeds, GrB_UINT64, n));
    GRB_TRY(GrB_Matrix_new(&C_scaled, GrB_INT64, n, n));
    GRB_TRY(GrB_Vector_new(&cs_row_vals, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&empty, GrB_FP64, n));

    GRB_TRY(GrB_Matrix_new(&C_temp_int64, GrB_INT64, n, n));

    GRB_TRY(GrB_UnaryOp_new(&ROUND_HUN_FP64, round_to_hundredths, GrB_FP64, GrB_FP64));
    GRB_TRY(GrB_UnaryOp_new(&hash_unary_op, fp64_floor_to_int64, GrB_INT64, GrB_FP64));


    // GxB_print(C_temp, GxB_COMPLETE);

    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&D, ones, 0));

    // Add self-edge to each vertex
    if (G->nself_edges != n)
    {
        GRB_TRY(GrB_assign(A, A, NULL, D, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;
        G->out_degree = NULL;
        G->in_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;
        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_InDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
    }

    GRB_TRY(GrB_Matrix_dup(&C_temp, A));
    GRB_TRY(GrB_Matrix_dup(&C, C_temp));

    double norm = 0;
    double last_norm = 0;
    GrB_Index streak = 0;
    GrB_Index iter = 0;
    while (true)
    {
        // Normalization step
        GRB_TRY(GrB_reduce(w, NULL, NULL, GrB_PLUS_MONOID_FP64, C_temp, GrB_DESC_RT0)); // Reduce across columns
        GRB_TRY(GrB_apply(w, NULL, NULL, GrB_MINV_FP64, w, GrB_DESC_R));
        GRB_TRY(GrB_Matrix_diag(&W, w, 0));
        // GxB_print(W, GxB_COMPLETE);
        GRB_TRY(GrB_mxm(C_temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C_temp, W, GrB_DESC_R));
        // printf("C_temp after normalization\n");
        // GxB_print(C_temp, GxB_COMPLETE);

        // Prune values less than some small threshold
        GRB_TRY(GrB_select(C_temp, NULL, NULL, GrB_VALUEGT_FP64, C_temp, threshold, NULL));

        bool res = NULL;
        LAGRAPH_TRY(LAGraph_Matrix_IsEqual(&res, C, C_temp, msg));
        if (res || iter > max_iter || streak >= 10)
        {
            printf("Terminated after %i iterations\n", iter);
            break;
        }

        GRB_TRY(GrB_Matrix_dup(&C, C_temp));

        // Expansion step
        for (int i = 0; i < e - 1; i++)
        {
            GRB_TRY(GrB_mxm(C_temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C_temp, C_temp, NULL));
        }

        // Inflation step
        GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP64(C_temp, NULL, NULL, GxB_POW_FP64, C_temp, (double)i, NULL));
        // printf("C_temp after inflation\n");
        // GxB_print(C_temp, GxB_COMPLETE);

        GRB_TRY(GrB_reduce(&norm, NULL, GrB_PLUS_MONOID_FP64, C_temp, NULL));

        if (norm == last_norm)
        {
            streak++;
        }
        else
        {
            streak = 0;
        }

        last_norm = norm;

        printf("End of iteration %i\nTotal = %f\nStreak = %i\n\n", iter, norm, streak);
        iter++;
    }

    // Post-processing searching for clusters

    for (int i = 0; i < n; i++)
    {
        uint64_t seed = (uint64_t)rand() % 1000 + 10;
        GRB_TRY(GrB_Vector_setElement(col_seeds, seed, i));
    }

    GRB_TRY(GrB_Vector_wait(col_seeds, GrB_MATERIALIZE));
    // GRB_TRY(GxB_Vector_sort(col_seeds, NULL, GrB_LT_INT64, col_seeds, NULL));
    // GRB_TRY(GrB_Vector_wait(col_seeds, GrB_MATERIALIZE));

    // GxB_print(col_seeds, GxB_SHORT);

    // GRB_TRY(GrB_apply(C_temp, NULL, NULL, ROUND_HUN_FP64, C_temp, NULL));
    GRB_TRY(GrB_Matrix_apply_BinaryOp1st_FP64(C_temp, NULL, NULL, GrB_TIMES_FP64, 100.0, C_temp, NULL));
    GRB_TRY(GrB_apply(C_temp_int64, NULL, NULL, hash_unary_op, C_temp, NULL));

    GxB_print(C_temp_int64, GxB_SHORT);


    GRB_TRY(GrB_Matrix_diag(&D, col_seeds, 0));
    GRB_TRY(GrB_mxm(C_scaled, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_INT64, C_temp_int64, D, NULL));

    // EXPERIMENTAL!
    // GxB_print(C_scaled, GxB_COMPLETE);
    // GRB_TRY(GrB_apply(C_scaled, NULL, NULL, hash_unary_op, C_scaled, NULL));
    // GxB_print(C_scaled, GxB_COMPLETE);

    GRB_TRY(GrB_reduce(cs_row_vals, NULL, NULL, GrB_PLUS_MONOID_INT64, C_scaled, NULL));

    GxB_print(cs_row_vals, GxB_COMPLETE);


    GrB_Index csrv_nvals;
    GRB_TRY(GrB_Vector_nvals(&csrv_nvals, cs_row_vals));

    GrB_Vector vals = NULL;
    GrB_Vector perm = NULL;
    GRB_TRY(GrB_Vector_new(&vals, GrB_INT64, n));
    GRB_TRY(GrB_Vector_new(&perm, GrB_INT64, n));
    GRB_TRY(GxB_Vector_sort(vals, perm, GrB_LT_INT64, cs_row_vals, NULL));

    GxB_print(C_temp, GxB_COMPLETE);
    GxB_print(C_scaled, GxB_SHORT);
    GxB_print(cs_row_vals, GxB_COMPLETE);
    GxB_print(vals, GxB_COMPLETE);
    GxB_print(perm, GxB_COMPLETE);

    int64_t* valsX;
    GrB_Index* permX;

    LAGRAPH_TRY(LAGraph_Malloc((void**)&permX, csrv_nvals, sizeof(GrB_Index), msg));
    LAGRAPH_TRY(LAGraph_Malloc((void**)&valsX, csrv_nvals, sizeof(int64_t), msg));

    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, valsX, &csrv_nvals, vals));
    GRB_TRY(GrB_Vector_extractTuples_INT64(NULL, permX, &csrv_nvals, perm));

    GrB_Index dupsIt = 0; // dups iterator
    GrB_Index* dupsIdx;
    LAGRAPH_TRY(LAGraph_Malloc((void**)&dupsIdx, csrv_nvals, sizeof(GrB_Index), msg));

    for (int i = 0; i < csrv_nvals; i++)
    {
        printf("  %i  ", valsX[i]);
    }
    printf("\n\n\n");

    double curVal = valsX[0];
    for (int i = 1; i < csrv_nvals; i++)
    {
        // printf("(valsX[%i]) %f =? %f (curVal)\n", i, valsX[i], curVal);
        if (valsX[i] == curVal)
        {
            // printf("ITERATION %i: Setting index %i to %i\n", i, dupsIt, permX[i]);
            dupsIdx[dupsIt] = permX[i];
            dupsIt++;
        }
        else
        {
            curVal = valsX[i];
        }
    }

    // for (int i = 0; i < csrv_nvals; i++)
    // {
    //     printf("  %i  ", indx[i]);
    // }
    // printf("\n\n\n");
    for (int i = 0; i < csrv_nvals; i++)
    {
        printf("  %i  ", dupsIdx[i]);
    }
    printf("\n\n\n");

    LAGRAPH_TRY(LAGraph_Realloc((void**)&dupsIdx, dupsIt, csrv_nvals, sizeof(GrB_Index), msg));

    int sum = 0;
    for (int i = 0; i < dupsIt; i++)
    {
        sum += dupsIdx[i];
        printf("  %i  ", dupsIdx[i]);
    }
    printf("\nsum = %i\n", sum);

    GxB_print(C_temp_int64, GxB_COMPLETE);


    // for (int i = 0; i < csrv_nvals; i++)
    // {
    //     if (keepX[i] == 0)
    //     {
    //         GRB_TRY(GrB_Row_assign(C_temp, NULL, NULL, empty, csrvI[i], GrB_ALL, n, GrB_DESC_R));
    //     }
    // }

    // CC matrix represents the actual clustering based on attractors. It can be
    // interpreted as CC[i][j] == 1 ==> vertex j is in cluster i
    GRB_TRY(GrB_Matrix_new(&CC, GrB_BOOL, n, n));


    GRB_TRY(GrB_assign(CC, C_temp, NULL, 1, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S));


    GRB_TRY(GrB_reduce(vpc, NULL, NULL, GrB_PLUS_MONOID_INT16, CC, NULL));

    // GxB_print(vpc, GxB_COMPLETE);
    // GxB_print(CC, GxB_COMPLETE);

    (*C_f) = C_temp; // Set output matrix

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}