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
#include <malloc.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <prng.h>
#include <unuran.h>
#include <unuran_tests.h>

#include <config.h>
#ifdef WITH_DMALLOC
#  include <dmalloc.h>
#endif

/*---------------------------------------------------------------------------*/

#if UNUR_URNG_TYPE != UNUR_URNG_PRNG
#error UNUR_URNG_TYPE must be set to UNUR_URNG_PRNG in unuran_config.h
#endif

/*---------------------------------------------------------------------------*/

#define CHI_TEST_INTERVALS 100  /* number of intervals for chi^2 test        */
#define CHI_TEST_VERBOSITY 0    /* verbosity level for chi^2 test: 0 | 1 | 2 */

/*---------------------------------------------------------------------------*/
/* enable/disable tests                                                      */
/* (comment out the methods or distributions that are not tested)            */

/* methods                                                                   */
#define T_DISTR
#define T_STRINGPARSER
#define T_CORDER
#define T_TIMING

#define T_AROU
#define T_CSTD
#define T_DARI
#define T_DAU
#define T_DGT
#define T_DSTD 
#define T_NINV
#define T_SROU
#define T_SSR
#define T_TABL
#define T_TDR
#define T_UNIF
#define T_UTDR

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

/* tresholds */
#define PVAL_LIMIT 1e-3             /* treshold for p-value for stat. test   */

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
/* check for invalid NULL pointer, that should not happen in this program */
void abort_if_NULL( FILE *LOG, int line, void *ptr );

/* compare error code */
int check_errorcode( FILE *LOG, int line, unsigned errno );

/* check for expected NULL pointer */
int check_expected_NULL( FILE *LOG, int line, void *ptr );

/* check for "set failed" */
int check_expected_setfailed( FILE *LOG, int line, int ok );

/* check for INFINITY */
int check_expected_INFINITY( FILE *LOG, int line, double x );

/* check for reinit */
int check_expected_reinit( FILE *LOG, int line, int ok );

/* check for non existing reinit */
int check_expected_no_reinit( FILE *LOG, int line, int ok );

/* compare double sequences generated by generator */
int compare_double_sequence_par_start( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size );
int compare_double_sequence_urng_start( FILE *LOG, int line, struct prng *urng, int sample_size );
int compare_double_sequence_par( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size );

int compare_double_sequence_gen_start( FILE *LOG, int line, struct prng *urng, UNUR_GEN *gen, int sample_size );
int compare_double_sequence_gen( FILE *LOG, int line, struct prng *urng, UNUR_GEN *gen, int sample_size );

/* compare int sequences generated by generator */
int compare_int_sequence_par_start( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size );
int compare_int_sequence_par( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size );
int compare_int_sequence_gen_start( FILE *LOG, int line, struct prng *urng, UNUR_GEN *gen, int sample_size );
int compare_int_sequence_gen( FILE *LOG, int line, struct prng *urng, UNUR_GEN *gen, int sample_size );

/* print name of distribution */
void print_distr_name( FILE *LOG, UNUR_DISTR *distr, const char *genid );

/* check p-value of statistical test and print result */
int print_pval( FILE *LOG, UNUR_GEN *gen, double pval, int trial, char todo );

/* run chi2 test */
int run_validate_chi2( FILE *LOG, int line, UNUR_GEN *gen, char todo );

/* run verify hat test */
int run_validate_verifyhat( FILE *LOG, int line, UNUR_GEN *gen, char todo );

/* print result of verify hat test */
int print_verifyhat_result( FILE *LOG, UNUR_GEN *gen, int failed, char todo );

/* print result of timings */
void print_timing_results( FILE *LOG, int line, UNUR_DISTR *distr, double *timing_result, int n_results );

/*---------------------------------------------------------------------------*/
