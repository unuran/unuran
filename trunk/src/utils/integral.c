/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: integral.c                                                        *
 *                                                                           *
 *   Routines for integral computations.                                     *
 *                                                                           *
 *****************************************************************************
 
 *****************************************************************************
 *                                                                           *
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
				                                                                                    
/*--------------------------------------------------------------------------*/

#include <unur_source.h>
#include "integral_source.h"
#include "../methods/x_gen.h"

/*---------------------------------------------------------------------------*/

#define UNUR_INTEGRAL_DEFAULT_N_SAMPLES (1e4)
#define UNUR_INTEGRAL_DEFAULT_N_REPETITIONS (100)

/*---------------------------------------------------------------------------*/

int
_unur_integral_evaluate_single (int dim, struct unur_funct_vgeneric *f, UNUR_GEN *gen, 
                                long n_samples, double *integral)
     /*----------------------------------------------------------------------*/
     /* Integrating f(x) in R^dim with x having a distribution given by the  */
     /* generator gen                                                        */
     /*                                                                      */
     /* input:								     */
     /*    dim  ... dimension of the arguments                               */
     /*    f    ... structure containing function to be integrated           */
     /*    gen  ... generator of the arguments                               */
     /*    n_samples ... number of samples to use in the integral calculation*/
     /*                                                                      */
     /* output:								     */
     /*   integral ... approximate value of the integral using n_samples     */
     /*                random points                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double *x;
  double sum=0;
  long i; 
  
  /* checking input */
  if (n_samples <= 0) {
    n_samples = UNUR_INTEGRAL_DEFAULT_N_SAMPLES; 
    _unur_warning(gen->genid , UNUR_ERR_PAR_SET, "n_samples<=0");
  }

  if (dim <= 0) {
    _unur_error(gen->genid, UNUR_ERR_PAR_SET, "dim<=0");
    return UNUR_ERR_PAR_SET;
  }
    
  CHECK_NULL( gen, UNUR_ERR_NULL );
  CHECK_NULL( f, UNUR_ERR_NULL );
  CHECK_NULL( f->f, UNUR_ERR_NULL );
    
  /* allocating memory for the arguments */
  x = _unur_xmalloc(dim*sizeof(double)); 
  
  /* calculating the integral-sum */
  for (i=0; i<n_samples; i++) {
    unur_sample_vec(gen, x);
    sum += f->f(x, f->params);
  }
  
  *integral = sum / n_samples;
  
  free(x);
  
  return UNUR_SUCCESS;  
} /* end of _unur_integral_evaluate() */

/*---------------------------------------------------------------------------*/

int
_unur_integral_evaluate (int dim, struct unur_funct_vgeneric *f, UNUR_GEN *gen, 
			 long n_samples, long n_repetitions, double *exact_value,
			 double *mean, double *variance, double *rmse)
     /*----------------------------------------------------------------------*/
     /* Integrating f(x) in R^dim with x having a distribution given by the  */
     /* generator gen                                                        */
     /*                                                                      */
     /* input:								     */
     /*    dim  ... dimension of the arguments                               */
     /*    f    ... structure containing function to be integrated           */
     /*    gen  ... generator of the arguments                               */
     /*    n_samples     ... number of samples to use in the calculation     */
     /*    n_repetitions ... number of repetitions of the calculation        */
     /*    exact_value   ... of integral if known, NULL otherwise            */
     /*                                                                      */
     /* output:								     */
     /*    mean     ... mean value of integral values (of all repetitions)   */
     /*    variance ... variance of integral values   (of all repetitions)   */
     /*    rmse     ... rmse if exact_value of integral known, 0 otherwise   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  long i;
  double *integral;
  double diff, sum_diff;
  
  /* checking input */
  if (n_samples <= 0) {
    n_samples = UNUR_INTEGRAL_DEFAULT_N_SAMPLES; 
    _unur_warning(gen->genid , UNUR_ERR_PAR_SET, "n_samples<=0");
  }

  if (n_repetitions <= 1) {
    n_repetitions = UNUR_INTEGRAL_DEFAULT_N_REPETITIONS; 
    _unur_warning(gen->genid , UNUR_ERR_PAR_SET, "n_repetitions<=1");
  }
  
  if (dim <= 0) {
    _unur_error(gen->genid, UNUR_ERR_PAR_SET, "dim<0");
    return UNUR_ERR_PAR_SET;
  }
    
  CHECK_NULL( gen, UNUR_ERR_NULL );
  CHECK_NULL( f, UNUR_ERR_NULL );
  CHECK_NULL( f->f, UNUR_ERR_NULL );

  /* allocating memory for the calculated integrals */
  integral = _unur_xmalloc(n_repetitions*sizeof(double)); 
  
  /* repeating the calculations and obtaining the mean value */
  *mean = 0;
  for (i=0; i<n_repetitions; i++) {
    _unur_integral_evaluate_single (dim, f, gen, n_samples, &integral[i]);
    *mean += integral[i];
  }
  *mean /= n_repetitions;
  
  /* calculating the variance using the corrected two-pass algorithm (num. rec.) */
  sum_diff = 0;
  for (i=0; i<n_repetitions; i++) {
    diff = integral[i] - *mean;
    sum_diff += diff;
  }
  
  *variance = 0;
  for (i=0; i<n_repetitions; i++) {
    diff = integral[i] - *mean;
    *variance += diff * diff;
  }
  *variance -= sum_diff*sum_diff / n_repetitions;
  *variance /= (n_repetitions -1);
  
  /* calculating the RMSE if exact_value for integram is given */
  /* othervise, 0 is returned ...  */
  *rmse = 0;  
  if (exact_value != NULL) {
    *rmse = 0;
    for (i=0; i<n_repetitions; i++) {
      diff = integral[i] - *exact_value;
      *rmse += diff * diff;
    }
    *rmse = sqrt(*rmse/n_repetitions);
  }
  
  free(integral);
  
  return UNUR_SUCCESS;  
} /* end of _unur_integral_evaluate() */
