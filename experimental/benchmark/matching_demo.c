#include "../../src/benchmark/LAGraph_demo.h"
#include "LG_internal.h"
#include "LAGraphX.h"

int main (int argc, char** argv)
{
    char msg [LAGRAPH_MSG_LEN] ;

    LAGraph_Graph G = NULL ;
    GrB_Matrix E = NULL ;
    GrB_Vector matching = NULL ;

    bool burble = false;
    demo_init (burble) ;

    LAGRAPH_TRY (readproblem (&G, NULL,
        true, true, true, NULL, false, argc, argv)) ;
    
    GRB_TRY (LAGraph_A_to_E (&E, G, msg)) ;
    GrB_Index num_edges ;
    GRB_TRY (GrB_Matrix_nrows (&num_edges, E)) ;
    // GRB_TRY (GrB_Vector_new (&matching, GrB_BOOL, num_edges)) ;
    // LAGRAPH_TRY (LAGraph_Graph_Print (G, LAGraph_SHORT, stdout, msg)) ;
    printf("printing E now: \n");
    LAGRAPH_TRY (LAGraph_Matrix_Print (E, LAGraph_SHORT, stdout, msg)) ;
    printf("running max matching now...\n");
    // LAGRAPH_TRY (LAGraph_MaximalMatching (&matching, E, 0, 5, msg)) ;
    return (GrB_SUCCESS) ;
}