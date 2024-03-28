//------------------------------------------------------------------------------
// LAGraph/experimental/benchmark/lcc_demo.c:
// benchmark for community detection using label propagation
//------------------------------------------------------------------------------

// LAGraph, (c) 2023 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// For additional details (including references to third party source code and
// other files) see the LICENSE file or contact permission@sei.cmu.edu. See
// Contributors.txt for a full list of contributors. Created, in part, with
// funding and support from the U.S. Government (see Acknowledgments.txt file).
// DM22-0790

// Contributed by Pascal Costanza, Intel, Belgium
// Based on tcc_demo by Tim Davis, Texas A&M

//------------------------------------------------------------------------------

// Usage:  peer_pressure_demo < matrixmarketfile.mtx
//         peer_pressure_demo matrixmarketfile.mtx
//         peer_pressure_demo matrixmarketfile.grb

#include "../../src/benchmark/LAGraph_demo.h"
#include "LAGraphX.h"
#include "LG_Xtest.h"

#define NTHREAD_LIST 1
#define THREAD_LIST 0

#define LG_FREE_ALL                 \
{                                   \
    LAGraph_Delete (&G, NULL) ;     \
    GrB_free (&c) ;                 \
    GrB_free (&C) ;                 \
}

int main (int argc, char **argv)
{

    //--------------------------------------------------------------------------
    // initialize LAGraph and GraphBLAS
    //--------------------------------------------------------------------------

    char msg [LAGRAPH_MSG_LEN] ;

    GrB_Vector c = NULL ;
    GrB_Matrix C = NULL ;
    LAGraph_Graph G = NULL ;

    // start GraphBLAS and LAGraph
    bool burble = false ;
    demo_init (burble) ;

    //--------------------------------------------------------------------------
    // read in the graph
    //--------------------------------------------------------------------------

    char *matrix_name = (argc > 1) ? argv [1] : "stdin" ;
    LAGRAPH_TRY (readproblem (&G, NULL,
                              false, true, true, NULL, false, argc, argv)) ;

    GrB_Index n, nvals ;
    GRB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GRB_TRY (GrB_Matrix_nvals (&nvals, G->A)) ;

    //--------------------------------------------------------------------------
    // initializations 
    //--------------------------------------------------------------------------

    GRB_TRY(GrB_Matrix_new(&C, GrB_BOOL, n, n));


    //--------------------------------------------------------------------------
    // run peer pressure clustering algorithm
    //--------------------------------------------------------------------------

    // compute check result
    double tt = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_PeerPressureClustering(&c, true, false, 0.0001, 50, G, msg)) ;
    tt = LAGraph_WallClockTime() - tt ;
    printf ("peer pressure run time %g sec\n", tt) ;

    // GxB_print(c, GxB_SHORT);

    double cov, perf, mod ;
    tt = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_PartitionQuality(&cov, &perf, c, G->A, msg)) ;
    tt = LAGraph_WallClockTime() - tt ;
    printf ("\npartition quality run time %g sec\n\tcoverage    = %f\n\tperformance = %f\n", tt, cov, perf) ;

    tt = LAGraph_WallClockTime ( ) ;
    LAGRAPH_TRY (LAGr_Modularity(&mod, (double)1, c, G->A, msg)) ;
    tt = LAGraph_WallClockTime() - tt ;
    printf ("modularity run time %g sec\n\tmodularity  = %f\n", tt, mod) ;

    //--------------------------------------------------------------------------
    // write cluster vector and adjacency matrix to files
    //--------------------------------------------------------------------------

    char *fn = "data/cluster_matrix.mtx";
    FILE *f = fopen(fn, "w");
    if (f == NULL) {
        fprintf(stderr, "Error opening file '%s' for writing\n", fn);
        return -1;
    }
    
    LAGRAPH_TRY(LAGraph_MMWrite((GrB_Matrix)c, f, NULL, msg));

    fclose(f);


    LG_FREE_ALL ;
    LAGRAPH_TRY (LAGraph_Finalize (msg)) ;
    return (GrB_SUCCESS) ;
}
