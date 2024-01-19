#define LG_FREE_WORK     \
    {                    \
        GrB_free(&C);    \
        GrB_free(&w);    \
        GrB_free(&I);    \
        GrB_free(&ones); \
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
    GrB_Matrix *C_f,
    // input
    int e,           // expansion coefficient
    int i,           // inflation coefficient
    LAGraph_Graph G, // input graph
    char *msg)
{
    char MATRIX_TYPE[LAGRAPH_MSG_LEN];

    GrB_Matrix C = NULL;      // Cluster workspace matrix
    GrB_Matrix C_temp = NULL; // The newly computed cluster matrix at the end of each loop

    GrB_Vector w = NULL; // weight vector to normalize C matrix

    GrB_Matrix I = NULL;    // Identity matrix
    GrB_Vector ones = NULL; // Vector of all 1's, used primarily to create identity

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

    printf("A\n");

    GRB_TRY(GrB_Matrix_new(&C, GrB_FP64, n, n));
    GRB_TRY(GrB_Matrix_new(&C_temp, GrB_FP64, n, n));
    GRB_TRY(GrB_Vector_new(&w, GrB_INT64, n));
    GRB_TRY(GrB_Matrix_new(&I, GrB_INT64, n, n));
    GRB_TRY(GrB_Matrix_new(&ones, GrB_INT64, n, n));

    GxB_print(C_temp, GxB_COMPLETE);

    printf("A\n");

    GRB_TRY(GrB_assign(ones, NULL, NULL, 1, GrB_ALL, n, NULL));
    GRB_TRY(GrB_Matrix_diag(&I, ones, 0));

    printf("A\n");

    // Add self-edge to each vertex
    if (G->nself_edges != n)
    {
        GRB_TRY(GrB_assign(A, A, NULL, I, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SC));
        G->A = A;
        G->out_degree = NULL;
        G->nself_edges = LAGRAPH_UNKNOWN;
        LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));
        LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));
    }

    printf("A\n");

    GxB_print(A, GxB_COMPLETE);

    // GRB_TRY(GrB_reduce(w, NULL, NULL, GrB_PLUS_MONOID_INT32, A, NULL));
    printf("APPLES %i\n", n); // Reduce along columns
    // GxB_print(w, GxB_COMPLETE);

    // GrB_Index iter = 0;
    // while (true)
    // {
    //     GrB
    // }

    printf("A\n");

    (*C_f) = C_temp; // Set output matrix

    GxB_print(*C_f, GxB_COMPLETE);


    LG_FREE_WORK;

    printf("B\n");

    return (GrB_SUCCESS);
}