/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      quantiles.c                                                  *
 *                                                                           *
 *   compute estimate for quartiles and mean of samples                      *
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
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] R. Jain, I. Chlamtac (1985): The P^2 Algorithm for Dynamic          *
 *       Calculation of Quantiles and Histograms Without Storing             *
 *       Observations, Comm. ACM 28(19), pp. 1076-1085.                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/
static char test_name[] = "Quantiles";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_test_quartiles( UNUR_GEN *gen, double *q0 ,double *q1, double *q2, double *q3, double *q4, int samplesize )
     /*----------------------------------------------------------------------*/
     /*  compute estimate for quartiles and mean of samples.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   q0         ... minimum                                             */
     /*   q1         ... 1st quartile (25%)                                  */
     /*   q2         ... mean (2nd quartile, 50%)                            */
     /*   q3         ... 3rd quartile (75%)                                  */
     /*   q4         ... maximum                                             */
     /*   samplesize ... sample size                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  double x = 0.;
  int n;

  /* check parameter */
  _unur_check_NULL(test_name,gen,0);
  /* type of distribution */
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute quartiles for distribution");
    return 0;
  }

  /* sample size >= 10 */
  if (samplesize < 10) 
    samplesize = 10;

  /* sampling */
  /* estimate quartiles using algorithm [1]. */
  
  for (n=1; n<=samplesize; n++) {

    /* which type of distribution */
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = (double)(_unur_sample_discr(gen)); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(gen); break;
    }

    /* ....................... */

  }

  /* now print results */
  printf("\nQuartiles:\n");
  printf("\tmin = %g\n",*q0);
  printf("\t25%% =\t%g\n",*q1);
  printf("\t50%% =\t%g\n",*q2);
  printf("\t75%% =\t%g\n",*q3);
  printf("\tmax = %g\n",*q4);

  return 1;

} /* end of unur_test_quartiles() */

/*---------------------------------------------------------------------------*/





