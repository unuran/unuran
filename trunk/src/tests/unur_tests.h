/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_tests.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for testing routines.      *
 *                                                                           *
 *   USAGE:                                                                  *
 *         included in all test source files.                                *
 *         required for every application of UNURAN test routines.           *
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
#ifndef __UNUR_TESTS_H_SEEN
#define __UNUR_TESTS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>
#include <unur_methods.h>

/*---------------------------------------------------------------------------*/
/* run battery of tests                                                      */
void unur_run_tests( struct unur_par *par, unsigned long tests,
		     double (*cdf)(double x, double *fparam, int n_fparam));

/*---------------------------------------------------------------------------*/
/* particular tests                                                          */

/* print a sample                                                            */
void unur_test_printsample( struct unur_gen *gen, int n_rows, int n_cols );

/* timing                                                                    */
struct unur_gen *unur_test_timing( struct unur_par *par, int log_samplesize );

/* count used uniform random numbers                                         */
int unur_test_count_urn( struct unur_gen *gen, int samplesize );

/* Chi^2 tests                                                               */
double unur_test_chi2( struct unur_gen *gen, int intervals, int samplesize, int classmin, int output );

/* make scatterplot of generated numbers                                     */
int unur_make_scatterplot( struct unur_gen *gen );

/* possible tests                                                            */
#define UNUR_TEST_ALL      ~(0UL)     /* make all possible tests */
#define UNUR_TEST_TIME     0x0001UL   /* estimate time */
#define UNUR_TEST_N_URNG   0x0002UL   /* count number of urng (needs compiler switch) */
#define UNUR_TEST_CHI2     0x0004UL   /* run chi^2 test for goodness of fit */
#define UNUR_TEST_SAMPLE   0x0008UL   /* print a sample file */
#define UNUR_TEST_SCATTER  0x0010UL   /* make a scatter plot */

/*---------------------------------------------------------------------------*/
#endif  /* end of #ifdef __UNUR_TESTS_H_SEEN */
/*---------------------------------------------------------------------------*/

