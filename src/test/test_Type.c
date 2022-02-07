//------------------------------------------------------------------------------
// LAGraph/src/test/test_Type.c:  test LAGraph_*Type* methods
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

//------------------------------------------------------------------------------

#include "LAGraph_test.h"

//------------------------------------------------------------------------------
// global variables
//------------------------------------------------------------------------------

GrB_Type type = NULL ;
char name [LAGRAPH_MAX_NAME_LEN] ;
char msg [LAGRAPH_MSG_LEN] ;
GrB_Scalar s = NULL ;

typedef int myint ;

//------------------------------------------------------------------------------
// test_TypeName :  test LAGraph_NameOfType
//------------------------------------------------------------------------------

void test_TypeName  (void)
{
    OK (LAGraph_Init (msg)) ;

    OK (LAGraph_NameOfType  (name, GrB_BOOL, msg)) ;
    OK (strcmp (name, "bool")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT8, msg)) ;
    OK (strcmp (name, "int8_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT16, msg)) ;
    OK (strcmp (name, "int16_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT32, msg)) ;
    OK (strcmp (name, "int32_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_INT64, msg)) ;
    OK (strcmp (name, "int64_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT8, msg)) ;
    OK (strcmp (name, "uint8_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT16, msg)) ;
    OK (strcmp (name, "uint16_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT32, msg)) ;
    OK (strcmp (name, "uint32_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_UINT64, msg)) ;
    OK (strcmp (name, "uint64_t")) ;

    OK (LAGraph_NameOfType  (name, GrB_FP32, msg)) ;
    OK (strcmp (name, "float")) ;

    OK (LAGraph_NameOfType  (name, GrB_FP64, msg)) ;
    OK (strcmp (name, "double")) ;

    char typename [LAGRAPH_MAX_NAME_LEN] ;
    OK (GrB_Scalar_new (&s, GrB_INT32)) ;
    OK (LAGraph_ScalarTypeName (name, s, msg)) ;
    OK (strcmp (name, "int32_t")) ;
    TEST_CHECK (LAGraph_ScalarTypeName (NULL, s, msg) == GrB_NULL_POINTER) ;
    TEST_CHECK (LAGraph_ScalarTypeName (name, NULL, msg) == GrB_NULL_POINTER) ;

    name [0] = '\0' ;
    OK (GrB_Type_new (&type, sizeof (myint))) ;
    int result = LAGraph_NameOfType (name, type, msg) ;
    #if LG_SUITESPARSE
    printf ("\nSuiteSparse knows the type name: [%s]\n", name) ;
    TEST_CHECK (result == GrB_SUCCESS) ;
    OK (strcmp (name, "myint")) ;
    #else
    TEST_CHECK (result == GrB_NOT_IMPLEMENTED) ;
    printf ("\nmsg: %s\n", msg) ;
    #endif

    TEST_CHECK (LAGraph_NameOfType (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("\nmsg: %s\n", msg) ;

    TEST_CHECK (LAGraph_NameOfType (name, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_NameOfType (NULL, GrB_BOOL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    GrB_free (&s) ;
    GrB_free (&type) ;
    OK (LAGraph_Finalize (msg)) ;
}

//------------------------------------------------------------------------------
// test_TypeSize :  test LAGraph_SizeOfType
//------------------------------------------------------------------------------

void test_TypeSize (void)
{
    OK (LAGraph_Init (msg)) ;
    size_t size ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_BOOL, msg)) ;
    TEST_CHECK (size == sizeof (bool)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_INT8, msg)) ;
    TEST_CHECK (size == sizeof (int8_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_INT16, msg)) ;
    TEST_CHECK (size == sizeof (int16_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_INT32, msg)) ;
    TEST_CHECK (size == sizeof (int32_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_INT64, msg)) ;
    TEST_CHECK (size == sizeof (int64_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_UINT8, msg)) ;
    TEST_CHECK (size == sizeof (uint8_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_UINT16, msg)) ;
    TEST_CHECK (size == sizeof (uint16_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_UINT32, msg)) ;
    TEST_CHECK (size == sizeof (uint32_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_UINT64, msg)) ;
    TEST_CHECK (size == sizeof (uint64_t)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_FP32, msg)) ;
    TEST_CHECK (size == sizeof (float)) ;

    size = 0 ;
    OK (LAGraph_SizeOfType  (&size, GrB_FP64, msg)) ;
    TEST_CHECK (size == sizeof (double)) ;

    size = 0 ;
    OK (GrB_Type_new (&type, sizeof (myint))) ;
    int result = LAGraph_SizeOfType (&size, type, msg) ;
    #if LG_SUITESPARSE
    printf ("\nSuiteSparse knows the type size: [%g]\n", (double) size) ;
    TEST_CHECK (result == GrB_SUCCESS) ;
    TEST_CHECK (size == sizeof (myint)) ;
    #else
    TEST_CHECK (result == GrB_NOT_IMPLEMENTED) ;
    printf ("\nmsg: %s\n", msg) ;
    #endif

    TEST_CHECK (LAGraph_SizeOfType (NULL, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("\nmsg: %s\n", msg) ;

    TEST_CHECK (LAGraph_SizeOfType (&size, NULL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    TEST_CHECK (LAGraph_SizeOfType (NULL, GrB_BOOL, msg) == GrB_NULL_POINTER) ;
    printf ("msg: %s\n", msg) ;

    GrB_free (&type) ;
    OK (LAGraph_Finalize (msg)) ;
}

//-----------------------------------------------------------------------------
// TEST_LIST: the list of tasks for this entire test
//-----------------------------------------------------------------------------

TEST_LIST =
{
    { "TypeName", test_TypeName  },
    { "TypeSize", test_TypeSize  },
    // no brutal test needed
    { NULL, NULL }
} ;
