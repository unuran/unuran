/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      moments.c                                                    *
 *                                                                           *
 *   compute moments of samples                                              *
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

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/
static char test_name[] = "Moments";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_moments( struct unur_gen *gen, int n_moments, double *moments, int samplesize )
     /*----------------------------------------------------------------------*/
     /*  compute moments of samples.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   n_moments  ... number of moments to be calculated                  */
     /*   moments    ... array for storing moments                           */
     /*   samplesize ... maximal sample size                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  double x,xm;
  int i, mom;

  /* check parameter */
  _unur_check_NULL(test_name,gen,0);
  CHECK_NULL(moments,0);
  if (n_moments <= 0) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"number of moments < 1");
    return 0;
  }
  if (samplesize <= 0) 
    samplesize = 10;

  /* clear array of moments */
  moments[0] = 1.;  /* dummy field */
  for (mom = 1; mom <= n_moments; mom++ )
    moments[mom] = 0.;

  /* sample */

  /* which type of distribution */
  switch (gen->method & UNUR_MASK_TYPE) {

  case UNUR_METH_DISCR:
    for (i=0; i<samplesize; i++) {
      xm = x = (double) _unur_sample_discr(gen);
      for (mom = 1; mom <= n_moments; mom++ ) {
	moments[mom] += xm;
	xm *= x;
      }
    }
    break;

  case UNUR_METH_CONT:
    for (i=0; i<samplesize; i++) {
      xm = x = _unur_sample_cont(gen);
      for (mom = 1; mom <= n_moments; mom++ ) {
	moments[mom] += xm;
	xm *= x;
      }
    }
    break;

  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute moments for distribution");
    return 0;
  }

  /* compute moments */
  for (mom = 1; mom <= n_moments; mom++ )
    moments[mom] /= samplesize;

  /* now print results */
  printf("\nMOMENTS:\n");
  for (mom = 1; mom <= n_moments; mom++ )
    printf("\t[%d] =\t%g\n",mom,moments[mom]);
  printf("\n");

  return 1;

} /* end of unur_test_moments() */

/*---------------------------------------------------------------------------*/
