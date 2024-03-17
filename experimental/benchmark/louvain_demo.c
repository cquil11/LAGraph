#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_Xtest.h"

#define LG_FREE_ALL                     \
{                                       \
    LAGraph_Delete (&G, msg) ;          \
    GrB_free (&A) ;                     \
    GrB_free (&A_struct) ;              \
    GrB_free (&c_f) ;                   \
    GrB_free (&vpc) ;                   \
}

int main (int argc, char**argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix A = NULL ;
    GrB_Vector c_f, vpc = NULL;         // Clustering result vector
    GrB_Matrix A_struct = NULL;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    FILE *f1, *f2, *f3;

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, NULL,
        false, false, false, NULL, true, argc, argv)) ;

    // compute G->out_degree
    LAGRAPH_TRY (LAGraph_Cached_OutDegree (G, msg)) ;

    // compute G->in_degree, just to test it (not needed for any tests)
    LAGRAPH_TRY (LAGraph_Cached_InDegree (G, msg)) ;

    LAGRAPH_TRY(LAGraph_Cached_AT (G, msg));

    // compute G->nself_edges
    LAGRAPH_TRY (LAGraph_Cached_NSelfEdges (G, msg)) ;

    // GrB_Index n ;
    // GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;

    printf("Input Matrix:\n");
    GxB_print (G->A, GxB_SHORT);

    GRB_TRY(LAGraph_Matrix_Structure(&A_struct, G->A, msg));

    // Run Louvain Clustering algorithm
    printf("Running Louvain...\n");
    double tt = LAGraph_WallClockTime();
    int iters = 0;
    LAGRAPH_TRY(LAGraph_Louvain_LSMP(&c_f, G, 100, &iters, msg));
    tt = LAGraph_WallClockTime() - tt;
    printf("Louvain completed in %g sec\n", tt);

    GxB_print(c_f, GxB_SHORT);
    GxB_print(vpc, GxB_SHORT);

    GrB_Index nclusters, n;
    double avg_cluster_size = 0;
    GRB_TRY(GrB_Vector_nvals(&nclusters, vpc));
    GRB_TRY(GrB_Matrix_nrows(&n, G->A));
    avg_cluster_size = 1.0 * n / nclusters;
    printf("Average size of cluster: %f\n", avg_cluster_size);
    
    double cov, perf, mod;
    LAGr_PartitionQuality(&cov, &perf, c_f, G->A, msg);
    LAGr_Modularity(&mod, (double)1, c_f, G->A, msg);
    printf("\nCoverage: %f\nPerformance: %f\nModularity: %f\n", cov, perf, mod);


    char *o_file = "pp_out.mtx";
    f1 = fopen(o_file, "w");
    if (f1 == NULL)
    {
        printf("Error opening file %s\n", o_file);
        return -1;
    }
    LAGRAPH_TRY (LAGraph_MMWrite(c_f, f1, NULL, msg));

    char *oA_file = "pp_A_san.mtx";
    f2 = fopen(oA_file, "w");
    if (f2 == NULL)
    {
        printf("Error opening file %s\n", oA_file);
        return -1;
    }
    LAGRAPH_TRY (LAGraph_MMWrite(G->A, f2, NULL, msg));

    char *oAstruct_file = "pp_A_struct.mtx";
    f3 = fopen(oAstruct_file, "w");
    if (f3 == NULL)
    {
        printf("Error opening file %s\n", oAstruct_file);
        return -1;
    }
    LAGRAPH_TRY (LAGraph_MMWrite(G->A, f3, NULL, msg));

    fclose(f1);
    fclose(f2);
    fclose(f3);
    f1, f2, f3 = NULL;

    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}