#define LG_FREE_WORK                     \
    {                                    \
        GrB_free(&C);                    \
        GrB_free(&w);                    \
        GrB_free(&W);                    \
        GrB_free(&I);                    \
        GrB_free(&ones);                 \
        GrB_free(&CC);                   \
        GrB_free(&vpc);                   \
        LAGraph_Free((void *)&CI, NULL); \
        LAGraph_Free((void *)&CJ, NULL); \
        LAGraph_Free((void *)&CX, NULL); \
    }

#define LG_FREE_ALL    \
    {                  \
        LG_FREE_WORK;  \
        GrB_free(C_f); \
    }

#include "LG_internal.h"
#include <LAGraphX.h>

int LAGr_MarkovClustering(
    // output:
    GrB_Matrix* C_f,
    // input
    int e,            // expansion coefficient
    int i,            // inflation coefficient
    double threshold, // threshold for pruning values
    LAGraph_Graph G,  // input graph
    char* msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GxB_set(GxB_PRINT_1BASED, true);

    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // The newly computed cluster matrix at the end of each loop

    GrB_Vector w = NULL; // weight vector to normalize C matrix
    GrB_Matrix W = NULL;

    GrB_Matrix I = NULL;    // Identity workspace matrix
    GrB_Vector ones = NULL; // Vector of all 1's, used primarily to create identity

    GrB_Matrix CC = NULL;
    GrB_Vector* C_rows = NULL;

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
        GRB_TRY(GrB_Matrix_new(&I, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&ones, GrB_FP64, n));
    GRB_TRY(GrB_Vector_new(&vpc, GrB_INT16, n));

    // GxB_print(C_temp, GxB_COMPLETE);

    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&I, ones, 0));

    // Add self-edge to each vertex
    if (G->nself_edges != n)
    {
        GRB_TRY(GrB_assign(A, A, NULL, I, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;
        G->out_degree = NULL;
        G->in_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;
        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_InDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
    }

    // GxB_print(G->A, GxB_COMPLETE);
    // GxB_print(G->out_degree, GxB_COMPLETE);
    // GxB_print(G->in_degree, GxB_COMPLETE);

    GRB_TRY(GrB_Matrix_dup(&C_temp, A));
    GRB_TRY(GrB_Matrix_dup(&C, C_temp));

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
        if (res || iter > 1000)
        {
            printf("Terminated after %i iterations\n", iter);
            break;
        }

        GRB_TRY(GrB_Matrix_dup(&C, C_temp));

        // Expansion step
        for (int i = 0; i < e - 1; i++)
        {
            GRB_TRY(GrB_mxm(C_temp, NULL, NULL, GrB_PLUS_TIMES_SEMIRING_FP64, C_temp, C_temp, NULL));
            // printf("Expansion loop #%i\n", i);
        }

        // printf("C_temp after expansion\n");

        // GxB_print(C_temp, GxB_COMPLETE);

        // Inflation step
        GRB_TRY(GrB_Matrix_apply_BinaryOp2nd_FP64(C_temp, NULL, NULL, GxB_POW_FP64, C_temp, (double)i, NULL));
        // printf("C_temp after inflation\n");
        // GxB_print(C_temp, GxB_COMPLETE);

        // printf("End of iteration %i\n", iter);
        iter++;
    }

    // GxB_print(C_temp, GxB_COMPLETE);

    // Post-processing searching for clusters



    GrB_Index nvals_C;
    GRB_TRY(GrB_Matrix_nvals(&nvals_C, C_temp));

    LAGRAPH_TRY(LAGraph_Malloc((void**)&C_rows, n, sizeof(GrB_Vector), msg));

    bool to_remove[n];
    memset(to_remove, 0, sizeof(to_remove));

    for (int i = 0; i < n; i++)
    {
        GRB_TRY(GrB_Vector_new(&C_rows[i], GrB_FP64, n));

        GRB_TRY(GrB_Col_extract(C_rows[i], NULL, NULL, C_temp, GrB_ALL, n, i, GrB_DESC_T0));
        // GxB_print(C_rows[i], GxB_COMPLETE);
    }

    for (int i = 0; i < n; i++)
    {
        if (to_remove[i] == 1) continue;

        GrB_Vector cur = C_rows[i];
        for (int j = i + 1; j < n; j++)
        {
            if (to_remove[j] == 1) continue;

            bool res;
            LAGRAPH_TRY(LAGraph_Vector_IsEqual(&res, cur, C_rows[j], msg));
            if (res)
            {
                to_remove[j] = 1;
            }
        }
    }

    // Iterate over each row to clear the ones marked for removal
    for (GrB_Index i = 0; i < n; i++)
    {
        if (to_remove[i] == 1)
        {
            for (int j = 0; j < n; j++)
            {
                GRB_TRY(GrB_Matrix_removeElement(C_temp, i, j));
            }
        }
    }

    GxB_print(C_temp, GxB_SHORT);

    // for (int i = 0; i < n; i++)
    //     printf("  %i  ", to_remove[i]);

    // CC matrix represents the actual clustering based on attractors. It can be
    // interpreted as CC[i][j] == 1 ==> vertex j is in cluster i
    GRB_TRY(GrB_Matrix_new(&CC, GrB_BOOL, n, n));

    // LAGRAPH_TRY(LAGraph_Malloc((void**)&CI, nvals_C, sizeof(GrB_Index), msg));
    // LAGRAPH_TRY(LAGraph_Malloc((void**)&CJ, nvals_C, sizeof(GrB_Index), msg));
    // LAGRAPH_TRY(LAGraph_Malloc((void**)&CX, nvals_C, sizeof(double), msg));

    // GRB_TRY(GrB_Matrix_extractTuples_FP64(CI, CJ, CX, &nvals_C, C_temp));

    // for (int ii = 0; ii < nvals_C; ii++)
    // {
    //     GRB_TRY(GrB_Matrix_setElement_BOOL(CC, 1, CI[ii], CJ[ii]));
    // }

    GRB_TRY(GrB_assign(CC, C_temp, NULL, 1, GrB_ALL, n, GrB_ALL, n, GrB_DESC_S));

    // GRB_TRY(GrB_Matrix_wait(CC, GrB_MATERIALIZE));


    GxB_print(CC, GxB_COMPLETE);

    (*C_f) = C_temp; // Set output matrix

    LG_FREE_WORK;

    return (GrB_SUCCESS);
}