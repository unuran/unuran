/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran_tests.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for testing routines.      *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef UNURAN_TESTS_H_SEEN
#define UNURAN_TESTS_H_SEEN
/*---------------------------------------------------------------------------*/

/*
   =NODE  Testing Testing

   =UP TOP [75]

   =DESCRIPTION
     The following routines can be used to test the performance of the
     implemented generators and can be used to verify the implementions.
     They are declared in @file{unuran_tests.h} which has to be included.

   =END

*/

/*---------------------------------------------------------------------------*/
/* possible tests                                                            */
#define UNUR_TEST_ALL      (~0u)     /* make all possible tests */
#define UNUR_TEST_TIME     0x001u    /* estimate time */
#define UNUR_TEST_N_URNG   0x002u    /* count number of urng (needs compiler switch) */
#define UNUR_TEST_CHI2     0x004u    /* run chi^2 test for goodness of fit */
#define UNUR_TEST_SAMPLE   0x008u    /* print a sample file */

/*---------------------------------------------------------------------------*/
/* =ROUTINES */

/*---------------------------------------------------------------------------*/

void unur_run_tests( UNUR_PAR *parameters, unsigned tests);
/* 
   Run a battery of tests.
   The following tests are available (use @code{|} to combine these
   tests):
   @table @code
   @item UNUR_TEST_ALL
   run all possible tests.
   @item UNUR_TEST_TIME
   estimate generation times.
   @item UNUR_TEST_N_URNG
   count number of uniform random numbers
   @item UNUR_TEST_CHI2
   run chi^2 test for goodness of fit
   @item UNUR_TEST_SAMPLE
   print a small sample.
   @end table
   All these tests can be started individually (see below).
*/

/*---------------------------------------------------------------------------*/
/* particular tests                                                          */

void unur_test_printsample( UNUR_GEN *generator, int n_rows, int n_cols, FILE *out );
/* 
   Print a small sample with @var{n_rows} rows and @var{n_cols} columns.
   @var{out} is the output stream to which all results are written.
*/

UNUR_GEN *unur_test_timing( UNUR_PAR *parameters, int log_samplesize, 
			    double *time_setup, double *time_sample,
			    int verbosity, FILE *out );
/* 
   Timing. @var{parameters} is an parameter object for which setup
   time and marginal generation times have to be measured. The results
   are written into @var{time_setup} and @var{time_sample},
   respectively. @var{log_samplesize} is the common logarithm of the
   sample size that is used for timing. 

   If @var{verbosity} is TRUE then a small table is printed to
   the @code{stdout} with setup time, marginal generation time and
   average generation times for generating 10, 100, @dots{} random
   variates. All times are given in micro seconds and relative to 
   the generation times for the underlying uniform random number
   (using the UNIF interface) and an exponential distributed 
   random variate using the inversion method.

   The created generator object is returned.
   If a generator object could not be created successfully, then NULL
   is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.

   Notice: All timing results are subject to heavy changes. Reruning
   timings usually results in different results. Minor changes in 
   the source code can cause changes in such timings up to 25 percent.
*/

double unur_test_timing_uniform( const UNUR_PAR *parameters, int log_samplesize );
/* */

double unur_test_timing_exponential( const UNUR_PAR *parameters, int log_samplesize );
/* 
   Marginal generation times for the underlying uniform random number
   (using the UNIF interface) and an exponential distributed 
   random variate using the inversion method. These times are used in
   unur_test_timing() to compute the relative timings results.
*/

double unur_test_timing_total( const UNUR_PAR *parameters, int samplesize, double max_duration );
/* 
   Timing. @var{parameters} is an parameter object for which average
   times a sample of size @var{samplesize} (including setup) are
   estimated. Thus sampling is repeated and the median of these timings 
   is returned (in micro seconds). The number of iterations is computed
   automatically such that the total amount of time necessary for the
   test does not exceed @var{max_duration} (given in seconds).
   
   If an error occurs then @code{-1} is returned.

   Notice: All timing results are subject to heavy changes. Reruning
   timings usually results in different results. Minor changes in 
   the source code can cause changes in such timings up to 25 percent.
*/

int unur_test_count_urn( UNUR_GEN *generator, int samplesize, int verbosity, FILE *out );
/* 
   Count used uniform random numbers. It returns the total number of
   uniform random numbers required for a sample of non-uniform random
   variates of size @var{samplesize}. Counting uniform random numbers
   might not work for the chosen @code{UNUR_URNG_TYPE} in
   @file{unuran_config.h}. In this case @code{-1} is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

double unur_test_chi2( UNUR_GEN *generator, int intervals, int samplesize, int classmin,
		       int verbosity, FILE *out );
/* 
   Run a Chi^2 test with the @var{generator}. 
   The resulting p-value is returned.

   It works with discrete und continuous univariate distributions.
   For the latter the CDF of the distribution is required.

   @var{intervals} is the number of intervals that is used for
   continuous univariate distributions. @var{samplesize} is the size
   of the sample that is used for testing. If it is set to @code{0}
   then a sample of size @var{intervals}^2 is used (bounded to some
   upper bound).

   @var{classmin} is the minimum number of expected entries per
   class. If a class has to few entries then some classes are joined.

   @var{verbosity} controls the output of the routine. If it is set
   to @code{1} then the result is written to the output stream
   @var{out}. If it is set to @code{2} additionally the list of
   expected and observed data is printed. There is no output when it
   is set to @code{0}.
*/

int unur_test_moments( UNUR_GEN *generator, double *moments, int n_moments, int samplesize,
		       int verbosity, FILE *out );
/* 
   Computes the first @var{n_moments} central moments for a sample of
   size @var{samplesize}. The result is stored into the array
   @var{moments}.
   @var{n_moments} must be an integer between @code{1} and @code{4}.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

double unur_test_correlation( UNUR_GEN *generator1, UNUR_GEN *generator2,
			      int samplesize, int verbosity, FILE *out );
/* 
   Compute the correlation coefficient between streams from
   @var{generator1} and @var{generator2} for two samples of size
   @var{samplesize}.
   The resultung correlation is returned.

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

int unur_test_quartiles( UNUR_GEN *generator,
			 double *q0, double *q1, double *q2, double *q3, double *q4, 
			 int samplesize, int verbosity, FILE *out );
/* 
   Estimate quartiles of sample of size @var{samplesize}. 
   The resulting quantiles are stored in the variables @var{q}:
   @table @var
   @item q0
   minimum
   @item q1
   25%
   @item q2
   median (50%)
   @item q3
   75%
   @item q4
   maximum
   @end table

   If @var{verbosity} is TRUE the result is written to the output
   stream @var{out}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_TESTS_H_SEEN */
/*---------------------------------------------------------------------------*/

