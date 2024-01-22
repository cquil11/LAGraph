
#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_Xtest.h"

#define LG_FREE_ALL                     \
{                                       \
    LAGraph_Delete (&G, msg) ;          \
    GrB_free (&A) ;                     \
    GrB_free (&C_f) ;                     \
}

int main(int argc, char** argv)
{
    char msg[LAGRAPH_MSG_LEN];

    LAGraph_Graph G = NULL;
    GrB_Matrix A = NULL;
    GrB_Matrix C_f = NULL;         // Clustering result vector

    // start GraphBLAS and LAGraph
    bool burble = false;
    demo_init(burble);

    // FILE *f ;

    char* matrix_name = (argc > 1) ? argv[1] : "stdin";
    LAGRAPH_TRY(readproblem(&G, NULL,
        false, false, false, GrB_FP64, false, argc, argv));

    // compute G->out_degree
    LAGRAPH_TRY(LAGraph_Cached_OutDegree(G, msg));

    // compute G->in_degree, just to test it (not needed for any tests)
    LAGRAPH_TRY(LAGraph_Cached_InDegree(G, msg));

    // compute G->nself_edges
    LAGRAPH_TRY(LAGraph_Cached_NSelfEdges(G, msg));

    // GrB_Index n ;
    // GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    printf("Input Matrix:\n");
    // GxB_print (G->A, GxB_SHORT);

    // Run Markov Clustering algorithm
    GRB_TRY(LAGr_MarkovClustering(&C_f, 2, 2, 500, 0.0001, 1.0e-8, 250, G, msg));

    GxB_print(C_f, GxB_COMPLETE);

    // char *o_file = "pp_out.mtx";
    // f = fopen(o_file, "w");
    // if (f == NULL)
    // {
    //     printf("Error opening file %s\n", o_file);
    //     return -1;
    // }

    // // Write the C_f matrix to the file
    // LAGRAPH_TRY (LAGraph_MMWrite(C_f, f, NULL, msg));

    // // Close the file
    // fclose(f);
    // f = NULL;


    LG_FREE_ALL;
    LAGRAPH_TRY(LAGraph_Finalize(msg));
    return (GrB_SUCCESS);
}
