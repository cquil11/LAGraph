//------------------------------------------------------------------------------
// LG_ndiag: count the # of diagonal entries in a matrix
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

#define LAGraph_FREE_ALL    \
{                           \
    GrB_free (&M) ;         \
    GrB_free (&D) ;         \
    GrB_free (&d) ;         \
}

#include "LG_internal.h"

int LG_ndiag                // returns 0 if successful, < 0 if failure
(
    // output
    int64_t *ndiag,         // # of entries 
    // input
    GrB_Matrix A,           // matrix to count
    GrB_Type atype,         // type of A
    char *msg               // error message
)
{

    //--------------------------------------------------------------------------
    // extract the diagonal and count its entries
    //--------------------------------------------------------------------------

    GrB_Matrix D = NULL, M = NULL ;
    GrB_Vector d = NULL ;
    LG_CHECK (ndiag == NULL, -1, "ndiag is NULL") ;
    (*ndiag) = LAGRAPH_UNKNOWN ;

    GrB_Index nrows, ncols ;
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;
    GrB_Index n = LAGraph_MIN (nrows, ncols) ;

    #if defined ( GxB_SUITESPARSE_GRAPHBLAS )

        #if ( GxB_IMPLEMENTATION >= GxB_VERSION (5,0,2) )

            //------------------------------------------------------------------
            // SuiteSparse:GraphBLAS v5.0.2: use GxB_Vector_diag
            //------------------------------------------------------------------

            GrB_TRY (GrB_Vector_new (&d, atype, n)) ;
            GrB_TRY (GxB_Vector_diag (d, A, 0, NULL)) ;
            GrB_TRY (GrB_Vector_nvals (ndiag, d)) ;

        #else

            //------------------------------------------------------------------
            // SuiteSparse:GraphBLAS v4.x: use GxB_Select
            //------------------------------------------------------------------

            GrB_TRY (GrB_Matrix_new (&D, atype, nrows, ncols)) ;
            GrB_TRY (GxB_select (D, NULL, NULL, GxB_DIAG, A, NULL, NULL)) ;
            GrB_TRY (GrB_Matrix_nvals (ndiag, D)) ;

        #endif

    #else

        //----------------------------------------------------------------------
        // pure GrB version with no GxB extensions
        //----------------------------------------------------------------------

        GrB_TRY (GrB_Matrix_new (&M, GrB_BOOL, nrows, ncols)) ;
        GrB_TRY (GrB_Matrix_new (&D, atype, nrows, ncols)) ;
        for (int64_t i = 0 ; i < n ; i++)
        {
            // M (i,i) = true
            GrB_TRY (GrB_Matrix_setElement (M, (bool) true, i, i)) ;
        }

        // D<M,struct> = A
        GrB_TRY (GrB_assign (D, M, NULL, A, GrB_ALL, nrows, GrB_ALL, ncols,
            GrB_DESC_S)) ;
        GrB_TRY (GrB_Matrix_nvals (ndiag, D)) ;

    #endif

    LAGraph_FREE_ALL ;
    return (0) ;
}

