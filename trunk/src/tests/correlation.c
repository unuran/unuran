/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      correlation.c                                                *
 *                                                                           *
 *   compute correlation coefficient of two samples                          *
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
static char test_name[] = "Correlation";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

double
unur_test_correlation( UNUR_GEN *genx, UNUR_GEN *geny, int samplesize )
     /*----------------------------------------------------------------------*/
     /*  compute correlation coefficient of two samples.                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genx       ... pointer to generator object                         */
     /*   geny       ... pointer to another generator object                 */
     /*   samplesize ... sample size                                         */
     /*                                                                      */
     /* return:                                                              */
     /*     correlation coefficient                                          */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  double x  =0., y =0.;  /* contains random numbers           */ 
  double mx =0., my=0.;  /* mean values (analog to s1 in [1]  */
  double dx =0., dy=0.;  /* analog to delta in [1]            */
  double sx =0., sy=0.;  /* analog to s2 in [1]               */
  double sxy=0.;         /* analog to s2 in [1]               */
  double factor;
  int n;

  /* check parameter */
  _unur_check_NULL(test_name,genx,0);
  _unur_check_NULL(test_name,geny,0);
  /* type of distribution */
  if (! ( ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((genx->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return 0;
  }
  if (! ( ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_DISCR) ||
	  ((geny->method & UNUR_MASK_TYPE) == UNUR_METH_CONT) )) {
    _unur_error(test_name,UNUR_ERR_GENERIC,
         "dont know how to compute correlation coefficient for distribution");
    return 0;
  }

  /* sample size >= 10 */
  if (samplesize < 10) 
    samplesize = 10;


  /* sampling */  
  for (n=1; n<=samplesize; n++) {

    /* which type of distribution */
    switch (genx->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      x = (double)(_unur_sample_discr(genx)); break;
    case UNUR_METH_CONT:
      x = _unur_sample_cont(genx); break;
    }
    switch (geny->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      y = (double)(_unur_sample_discr(geny)); break;
    case UNUR_METH_CONT:
      y = _unur_sample_cont(geny); break;
    }

    factor = (double) ( n*(n-1) );

    dx = (x - mx) / n;
    dy = (y - my) / n;
    mx += dx;
    my += dy;

    sx  += factor * dx*dx;
    sy  += factor * dy*dy;
    sxy += factor * dx*dy;
    //printf( "sx:%f sy:%f sxy:%f, c:%f\n",sx,sy,sxy, sxy/(sqrt(sx*sy)) );

  }

  /* now print results */
  printf("\nCorrelation coefficient:\n");
    printf( "%f", sxy/(sqrt(sx*sy)) );
  printf("\n");

  return ( sxy/(sqrt(sx*sy)) );

} /* end of unur_test_correlation() */

/*---------------------------------------------------------------------------*/










