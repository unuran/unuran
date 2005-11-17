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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <unuran.h>
#include <unuran_tests.h>
#include "testdistributions/testdistributions.h"

#include <config.h>
#ifdef WITH_DMALLOC
#  include <dmalloc.h>
#endif

/*---------------------------------------------------------------------------*/

#define CHI_TEST_INTERVALS 100  /* number of intervals for chi^2 test        */
#define CHI_TEST_VERBOSITY 0    /* verbosity level for chi^2 test: 0 | 1 | 2 */

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

/* tresholds */
#define PVAL_LIMIT (1.e-3)                  /* treshold for p-value for stat. test */
#define DEFAULT_CHI2_FAILURES_TOLERATED (2) /* tolerated number of failed tests    */

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
/* compare two floating point numbers                                        */
/* (this is copied from ./src/utils/unur_fp_source.h)                        */

#define _unur_FP_equal(a,b) \
 ((a)==(b) || \
  fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b)) * 100. * FLT_EPSILON)


/*---------------------------------------------------------------------------*/
/* check for invalid NULL pointer, that should not happen in this program */
void abort_if_NULL( FILE *LOG, int line, const void *ptr );

/* compare error code */
int check_errorcode( FILE *LOG, int line, unsigned cherrno );

/* check for expected NULL pointer */
int check_expected_NULL( FILE *LOG, int line, const void *ptr );

/* check for "set failed" */
int check_expected_setfailed( FILE *LOG, int line, int rcode );

/* check for expected zero (int 0) */
int check_expected_zero( FILE *LOG, int line, int k );

/* check for INFINITY */
int check_expected_INFINITY( FILE *LOG, int line, double x );
int check_expected_negINFINITY( FILE *LOG, int line, double x );
int check_expected_INTMAX( FILE *LOG, int line, int k );

/* check for reinit */
int check_expected_reinit( FILE *LOG, int line, int rcode );

/* check for non existing reinit */
int check_expected_no_reinit( FILE *LOG, int line, int rcode );

/* compare double sequences generated by generator */
int compare_double_sequence_par_start( FILE *LOG, int line, UNUR_PAR *par, int sample_size );
int compare_double_sequence_urng_start( FILE *LOG, int line, int sample_size );
int compare_double_sequence_par( FILE *LOG, int line, UNUR_PAR *par, int sample_size );
int compare_double_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );
int compare_double_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );

/* compare int sequences generated by generator    */
int compare_int_sequence_par_start( FILE *LOG, int line, UNUR_PAR *par, int sample_size );
int compare_int_sequence_par( FILE *LOG, int line, UNUR_PAR *par, int sample_size );
int compare_int_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );
int compare_int_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );

/* compare sequences of double vectors generated by generator */
int compare_cvec_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );
int compare_cvec_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );

/* compare sequences of double matrices generated by generator */
int compare_matr_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );
int compare_matr_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size );

/* free memory used for comparing sequences */
void compare_free_memory( void );

/* print name of distribution */
void print_distr_name( FILE *LOG, const UNUR_DISTR *distr, const char *genid );

/* check p-value of statistical test and print result */
int print_pval( FILE *LOG, UNUR_GEN *gen, const UNUR_DISTR *distr, double pval, int trial, char todo );

/* run chi2 test */
int run_validate_chi2( FILE *LOG, int line, UNUR_GEN *gen, const UNUR_DISTR *distr, int testtype, char todo );

/* run verify hat test */
int run_validate_verifyhat( FILE *LOG, int line, UNUR_GEN *gen, const UNUR_DISTR *distr, char todo );

/* print result of verify hat test */
int print_verifyhat_result( FILE *LOG, UNUR_GEN *gen, const UNUR_DISTR *distr, int failed, char todo );

/* print result of timings */
void print_timing_results( FILE *LOG, int line, const UNUR_DISTR *distr,
			   double *timing_setup, double *timing_marginal, int n_results );

/*---------------------------------------------------------------------------*/
