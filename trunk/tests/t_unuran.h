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

#include <prng.h>
#include <unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

#if UNUR_URNG_TYPE != UNUR_URNG_PRNG
#error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h
#endif

/*---------------------------------------------------------------------------*/

#define CHI_TEST_INTERVALS 100

/*---------------------------------------------------------------------------*/
/* enable/disable tests                                                      */

/* methods                                                                   */
#define T_SROU

/* distributions                                                             */
#define D_NORMAL
#define D_GAMMA

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

/* tests */
extern int test_ok;                 /* all tests ok (boolean)                */
extern int test_failed;             /* failed tests                          */
extern FILE *TESTLOG;               /* test log file                         */

/* tresholds */
#define PVAL_LIMIT 1e-3             /* treshold for p-value for stat. test   */

/* distributions */
struct list_distr {
  UNUR_DISTR *distr;                /* pointer to distribution object        */
  unsigned    type;                 /* type of distribution                  */
  double      c_max;                /* maximal value for c to be T_c concave */
};

#define T_TYPE_TDR 0x00000001

extern struct list_distr *list_of_distr;  /* pointer to list of distributions */
extern int n_distr;                 /* number of distributions               */

/*---------------------------------------------------------------------------*/
/* True and false                                                            */

#ifndef TRUE
#define TRUE   (1)
#endif

#ifndef FALSE
#define FALSE  (0)
#endif

#ifndef M_LN10
#define M_LN10  2.302585092994046
#endif

/*---------------------------------------------------------------------------*/
/* make list of distributions                                                */
void make_list_of_distributions( void );

/*---------------------------------------------------------------------------*/
/* supporting routines                                                       */

/* compare error code */
#define check_errorcode(errorcode) \
  do {do_check_errorcode(__LINE__,(errorcode));} while (0)
inline void do_check_errorcode( int line, unsigned errno );

/* check for expected NULL pointer */
#define check_expected_NULL(ptr) \
  do {unur_errno = 0; do_check_expected_NULL(__LINE__,(ptr)); } while(0)
inline void do_check_expected_NULL( int line, void *ptr );

/* check for "set failed" */
#define check_expected_setfailed(ok) \
  do {unur_errno = 0; do_check_expected_setfailed(__LINE__,(ok)); } while(0)
inline void do_check_expected_setfailed( int line, int ok );

/* compare two sequences */
#define compare_sequences(a,b,n) \
  do {do_compare_sequences(__LINE__,(a),(b),(n)); } while(0)
void do_compare_sequences( int line, double *a, double *b, int n );

/* check p-value of statistical test */
#define check_pval(gen,pval,trial) \
  do {do_check_pval(__LINE__,(gen),(pval),(trial)); } while(0)
void do_check_pval( int line, UNUR_GEN *gen, double pval, int trial );

/* print p-value of statistical test */
void print_pval( double pval, int trial );

/*---------------------------------------------------------------------------*/
