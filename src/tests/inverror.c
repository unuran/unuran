/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file: inverror.c                                                        *
 *                                                                           *
 *   Estimate U-error for inversion methods                                  *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/x_gen_source.h>
#include <methods/hinv.h>

#include "unuran_tests.h"

/*---------------------------------------------------------------------------*/

static char test_name[] = "InvError";

/*---------------------------------------------------------------------------*/

double
unur_test_inverror( const UNUR_GEN *gen, 
		    double *max_error, double *MAE, double threshold,
		    int samplesize, int randomized, 
		    int verbosity, FILE *out )
     /*----------------------------------------------------------------------*/
     /* Estimate maximal u-error and mean absolute error (MAE) by means of   */
     /* (Quasi-) Monte-Carlo simulation.                                     */
     /* In addition a penalty value is computed and returned.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   max_error  ... pointer for storing maximal u-error                 */
     /*   MEA        ... pointer for storing MA u-error                      */
     /*   threshold  ... maximum allowed error                               */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*   randomized ... use pseudo-random (TRUE) or quasi-random (FALSE)    */
     /*   verbosity  ... verbosity level, 0 = no output, 1 = output          */
     /*   out        ... output stream                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   total number of used uniform random numbers                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define DISTR   gen->distr->data.cont

  double a = 1./samplesize;  /* multiplicator for QRNG */ 

  UNUR_FUNCT_CONT *cdf;      /* pointer to CDF */
  double CDFmin, CDFmax;     /* minimum and maximum of CDF in given domain */

  double (*quantile)(const UNUR_GEN *, double);  /* pointer to quantile function */

  double U, X;               /* uniform and non-uniform (Q)RN */
  double cdfX;               /* CDF at X */
  double uerror, umax, usum; /* last, maximum, sum of U-error(s) */
  double Xmax = INFINITY;    /* x-value with largest U-error */
  double penalty = 0;        /* score for error */

  int j;                     /* aux variable */

  /* check arguments */
  _unur_check_NULL(test_name,gen,-1.);
  if (verbosity) { _unur_check_NULL(test_name,out,-1.); }

  /* get pointer to function that approximates quantiles */
  switch (gen->method) {
  case UNUR_METH_HINV:
    quantile = unur_hinv_eval_approxinvcdf;
    break;

  default:
    _unur_error(test_name,UNUR_ERR_GENERIC,"inversion method required");
    return -1.;
  }

  /* CDF required */
  if (DISTR.cdf == NULL) {
    _unur_error(test_name,UNUR_ERR_GENERIC,"CDF required");
    return -2.;
  }
  cdf = DISTR.cdf;

  /* range of CDF */
  CDFmin = (DISTR.trunc[0] > -INFINITY) ? _unur_cont_CDF((DISTR.trunc[0]),(gen->distr)) : 0.;
  CDFmax = (DISTR.trunc[1] < INFINITY)  ? _unur_cont_CDF((DISTR.trunc[1]),(gen->distr)) : 1.;
  if (! _unur_FP_less(CDFmin,CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return -1;
  }

  /* initialize variables */
  umax = 0.;
  usum = 0.;

  /* get sample */
  for(j=0;j<samplesize;j++) {

    /* uniform random number */
    U = (randomized) ? _unur_call_urng(gen->urng) : (j+0.5)*a;

    /* compute inverse CDF */
    X = quantile(gen,U);

    /* compute (rescaled) CDF at X */
    cdfX = (_unur_cont_CDF(X,gen->distr) - CDFmin) / (CDFmax - CDFmin);


    /* compute U-error */
    uerror = fabs(U - cdfX);

    /* update error estimates */
    usum += uerror;
    if (uerror > umax) {
      umax = uerror;
      Xmax = X;
    }
    /* printf("j %d uerror %e maxerror %e average %e\n",j,uerror,max,average/(j+1)); */

    /* update penalty */
    if (_unur_FP_less(threshold,uerror))
      penalty += 1. + 10.*(uerror - threshold) / threshold;
  }

  /*
    printf("maximal error occured at x= %.16e percentage outside interval %e \n",
    errorat,(double)outside_interval/(double)samplesize);
  */

  /* save data */
  *max_error = umax;
  *MAE = usum/samplesize;

  /* return penalty */
  return penalty/samplesize;

#undef DISTR
} /* end of unur_test_estimate_inverror() */

/*---------------------------------------------------------------------------*/
