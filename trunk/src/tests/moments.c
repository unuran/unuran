/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      moments.c                                                    *
 *                                                                           *
 *   compute central moments of samples                                      *
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
 *   [1] Spicer C.C. (1972): Algorithm AS 52: Calculation of Power Sums of   *
 *       Deviations about the mean, Applied Statistics 21(2), pp. 226-227.   *
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
unur_test_moments( UNUR_GEN *gen, double *moments, int n_moments, int samplesize )
     /*----------------------------------------------------------------------*/
     /*  compute central moments of samples.                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   moments    ... array for storing moments                           */
     /*   n_moments  ... number of moments to be calculated (at most 4)      */
     /*   samplesize ... sample size                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  double x = 0.;
  double an, an1, dx, dx2;
  int n, mom;

  /* check parameter */
  _unur_check_NULL(test_name,gen,0);
  /* type of distribution */
  if (! ( ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((gen->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"dont know how to compute moments for distribution");
    return 0;
  }
  /* array for storing moments */
  CHECK_NULL(moments,0);
  if (n_moments <= 0 || n_moments > 4) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"number of moments < 1 or > 4");
    return 0;
  }

  /* sample size >= 10 */
  if (samplesize < 10) 
    samplesize = 10;

  /* clear array of moments */
  moments[0] = 1.;  /* dummy field */
  for (mom = 1; mom <= n_moments; mom++ )
    moments[mom] = 0.;

  /* sampling */
  /* compute moments: we use a recurrence relation by Spicer [1]. */
  
  for (n=1; n<=samplesize; n++) {

    /* which type of distribution */
    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = (double)(_unur_sample_discr(gen)); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(gen); break;
    }

    an = (double)n;
    an1 = an-1.;
    dx = (x - moments[1]) / an;
    dx2 = dx * dx;
   
    switch (n_moments) {
    case 4:
      moments[4] -= dx * (4.*moments[3] - dx * (6.*moments[2] + an1*(1. + an1*an1*an1)*dx2));
    case 3:
      moments[3] -= dx * (3.*moments[2] - an*an1*(an-2.)*dx2);
    case 2:
      moments[2] += an * an1 * dx2;
    case 1:
      moments[1] += dx;
    }
  }

  /* compute moments */
  for (mom = 1; mom <= n_moments; mom++ )
    moments[mom] /= samplesize;

  /* now print results */
  printf("\nCentral MOMENTS:\n");
  for (mom = 1; mom <= n_moments; mom++ )
    printf("\t[%d] =\t%g\n",mom,moments[mom]);
  printf("\n");

  return 1;

} /* end of unur_test_moments() */

/*---------------------------------------------------------------------------*/

