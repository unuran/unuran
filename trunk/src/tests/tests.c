/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      testw.c                                                      *
 *                                                                           *
 *   run various statistical tests, timing, ...                              *
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

#include <unuran_tests.h>
#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

#define TEST_TIMING_LOG_SAMPLESIZE 5    /* common log of sample size for timing (>= 2!!) */

#define TEST_COUNTER_SAMPLESIZE 100000  /* sample size for counting number of URNs */

#define TEST_SAMPLE_ROWS    3     /* number of rows and columns for sample   */
#define TEST_SAMPLE_COLS    10


#define TEST_CHI2_VERBOSE   1     /* output switch for chi^2 test            */
#define TEST_CHI2_INTERVALS 100   /* number of intervals for chi^2 test      */

/*---------------------------------------------------------------------------*/
/* some tools                                                                */
static int _unur_print_method( struct unur_par *par );

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Run different tests                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

void 
unur_run_tests( struct unur_par *par, unsigned tests)
     /*----------------------------------------------------------------------*/
     /* test generator                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to paramters for building generator object       */
     /*   tests ... list of tests (stored as bit array)                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen = NULL;

  /* check arguments */
  _unur_check_NULL("Tests",par,/*void*/);

  /* print info about method */
  if (! _unur_print_method(par))
    return;  /* unknown method */

  /* init generator object */
  if (tests & UNUR_TEST_TIME)
    /* evaluate setup time and generation time */
    gen = unur_test_timing(par,TEST_TIMING_LOG_SAMPLESIZE);
  else
    gen = _unur_init(par);

  /* init successful ? */
  if (!gen) return;

  /* count number of uniform random numbers */
  if (tests & UNUR_TEST_N_URNG )
    unur_test_count_urn(gen,TEST_COUNTER_SAMPLESIZE);

  /* print a sample */
  if (tests & UNUR_TEST_SAMPLE )
    unur_test_printsample(gen,TEST_SAMPLE_ROWS,TEST_SAMPLE_COLS);

  /* run chi2-test*/
  if (tests & UNUR_TEST_CHI2)
    unur_test_chi2(gen,TEST_CHI2_INTERVALS,0,0,TEST_CHI2_VERBOSE);

  /* make scatterplot */
  if (tests & UNUR_TEST_SCATTER)
     unur_make_scatterplot(gen);

  /* free generator */
  _unur_free(gen);

} /* end of unur_run_tests() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Some tools                                                              *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

static int 
_unur_print_method( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* print name of method                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* first print type */
  switch (par->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    printf("\nTYPE:\t\tdiscrete univariate distribution\n");
    break;
  case UNUR_METH_CONT:
    printf("\nTYPE:\t\tcontinuous univariate distribution\n");
    break;
  case UNUR_METH_VEC:
    printf("\nTYPE:\t\tcontinuous multivariate distribution\n");
    break;
  default: /* unknown ! */
    _unur_warning("Tests",UNUR_ERR_GENERIC,"type of method unknown!");
    return 0;
  }

  /* print method description */
  switch (par->method) {

    /* discrete, univariate */
  case UNUR_METH_DAU:
    COOKIE_CHECK(par,CK_DAU_PAR,0);
    printf("METHOD:\t\talias and alias-urn method (DAU)\n");
    break;
  case UNUR_METH_DIS:
    COOKIE_CHECK(par,CK_DIS_PAR,0);
    printf("METHOD:\t\tindexed search (DIS)\n");
    break;
  case UNUR_METH_DSTD:
    COOKIE_CHECK(par,CK_DSTD_PAR,0);
    printf("METHOD:\t\tspecial (DSTD)\n");
    break;

    /* continuous, univariate */
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    printf("METHOD:\t\tautomatic ratio-of-uniforms method\n");
    break;
  case UNUR_METH_NINV:
    COOKIE_CHECK(par,CK_NINV_PAR,0);
    printf("METHOD:\t\tnumerical inversion of c.d.f.\n");
    break;
  case UNUR_METH_SROU:
    COOKIE_CHECK(par,CK_SROU_PAR,0);
    printf("METHOD:\t\tsimple universal ratio-of-uniforms method\n");
    break;
  case UNUR_METH_STDR:
    COOKIE_CHECK(par,CK_STDR_PAR,0);
    printf("METHOD:\t\tsimple transformed density rejection with universal bounds\n");
    break;
  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    printf("METHOD:\t\trejection from piecewise constant hat\n");
    break;
  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    printf("METHOD:\t\ttransformed density rejection\n");
    break;
  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    printf("METHOD:\t\ttransformed density rejection, 3-point method (UTDR)\n");
    break;
  case UNUR_METH_RECT:
    COOKIE_CHECK(par,CK_RECT_PAR,0);
    printf("METHOD:\t\tuniformly distributed in hypercube (RECT)\n");
    break;
  case UNUR_METH_CSTD:
    COOKIE_CHECK(par,CK_CSTD_PAR,0);
    printf("METHOD:\t\tspecial (CSTD)\n");
    break;
  default: /* unknown ! */
    _unur_error("Tests",UNUR_ERR_GENERIC,"method unknown!");
    return 0;
  }
  
  /* everything o.k. */
  return 1;

} /* end of _unur_print_method() */

/*---------------------------------------------------------------------------*/


