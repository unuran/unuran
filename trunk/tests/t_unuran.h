/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests                                                                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <malloc.h>

#include <unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/
/* enable/disable tests                                                      */

/* methods                                                                   */
#define T_SROU

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

extern int test_ok;               /* all tests ok (boolean)                  */
extern int test_ok_local;         /* running test ok (boolean)               */
extern FILE *TESTLOG;             /* test log file                           */

/*---------------------------------------------------------------------------*/
/* True and false                                                            */

#ifndef TRUE
#define TRUE   (1)
#endif

#ifndef FALSE
#define FALSE  (0)
#endif

/*---------------------------------------------------------------------------*/
/* supportong routines                                                       */

/* compare error code */
#define compare_errno(errorcode) \
  do {do_compare_errno(__LINE__,(errorcode));} while (0)
inline void do_compare_errno( int line, unsigned errno );

/* check for expected NULL pointer */
#define check_expected_NULL(ptr) \
  do {unur_errno = 0; do_check_expected_NULL(__LINE__,(ptr)); } while(0)
inline void do_check_expected_NULL( int line, void *ptr );

/*---------------------------------------------------------------------------*/
