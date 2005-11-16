/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multiexponential.c                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Multiexponential distribution                                     *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [5] S. Kotz, N. Balakrishnan, and N.L. Johnson                          *
 *       Continuous Multivariate Distributions,                              *
 *       Volume 1: Models and Applications                                   *
 *       John Wiley & Sons, Inc., New York, 2000                             *
 *                                                                           *
 *       John E. Freund                                                      *
 *       A Bivariate Extension of the Exponential Distribution               *
 *       Journal of the American Statistical Association,                    *
 *       Vol. 56, p971-977 (1961)                                            *
 *                                                                           * 
 *       Asit P. Basu, Kai Sun                                               *
 *       Multivariate Exponential Distributions with Constant Failure Rates  *
 *       Journal of Multivariate Analysis,                                   *
 *       Vol. 61, p159-169 (1997)                                            *
 *                                                                           *
 *                                                                           *
 *  pdf:       f(x) = Prod_{i=0}^{i=dim-1}                                   *
 *             exp(-(dim-i) (x_{i}-x_{i-1} - (theta_i-theta_{i-1}) ) / sigma_i) *
 *             with x_{-1}=0 and thet_{i-1}=0                                *                  
 *  domain:    [0, inf)^(dim)                                                *
 *  constant:  Prod_{i=0}^{i=dim-1} 1/sigma_i                                *
 *                                                                           *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma   ... scale parameter     (default : (1,1,...,1)-vector)    *
 *     1:  theta   ... location parameter  (default : (0,0,...,0)-vector)    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp( - Sum_{i_0}^{i=dim-1} (dim-i) (x_{i}-x{i-1}) )    *
 *  domain:    [0, inf)^(dim)                                                *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *****************************************************************************
 
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

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cvec.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "multiexponential";

/*---------------------------------------------------------------------------*/
/* parameters */
#define INDEX_SIGMA 0 
#define INDEX_THETA 1

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

static double _unur_pdf_multiexponential( const double *x, UNUR_DISTR *distr );
static double _unur_logpdf_multiexponential( const double *x, UNUR_DISTR *distr );
static int _unur_dlogpdf_multiexponential( double *result, const double *x, UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_multiexponential( const double *x, UNUR_DISTR *distr )
{ 
  double fx;
  fx=_unur_logpdf_multiexponential( x, distr );
  if (_unur_isfinite(fx))  
    return exp( _unur_logpdf_multiexponential( x, distr ) );
  
  return 0.;
} /* end of _unur_pdf_multiexponential() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_multiexponential( const double *x, UNUR_DISTR *distr )
{ 
  int i, dim;
  
  double dx, sum; /* arguments used in the evaluation of the density */
  double *sigma, *theta;
  
  dim = distr->dim ;
  sum=0;
  
  sigma = DISTR.param_vecs[INDEX_SIGMA];
  theta = DISTR.param_vecs[INDEX_THETA];
  
  /* TODO : reconsider logic concerning NULL pointers ? */
  if (sigma == NULL || theta == NULL) {
    /* standard form */
    for (i=0; i<dim; i++) { 
      dx = (i==0) ? ((x[i]<0)? INFINITY: x[i]) : ((x[i]<x[i-1])? INFINITY: x[i]-x[i-1]);  
      if (_unur_isfinite(dx)) 
        sum -= (dim-i) * dx;  
      else
        sum = -INFINITY;
    }
  }
  else {
    for (i=0; i<dim; i++) {
    
      dx = (i==0) ? ((x[i]-theta[i]<0)? INFINITY: x[i]-theta[i]) : ((x[i]-theta[i]<x[i-1]-theta[i-1])? INFINITY: x[i]-x[i-1]-theta[i]+theta[i-1]); 
      
      /* sigma[i] is expected to be > UNUR_EPSILON here */
      dx /= sigma[i];
      if (_unur_isfinite(dx)) 
        sum -= (dim-i) * dx;  
      else
        sum = -INFINITY;      
      
/*printf(">> dx=%f sum=%f\n", dx, sum); */   
    }
  }
          
  return ( sum + LOGNORMCONSTANT);

} /* end of _unur_logpdf_multiexponential() */

/*---------------------------------------------------------------------------*/

int
_unur_dlogpdf_multiexponential( double *result, const double *x, UNUR_DISTR *distr )
{
  int i, j, dim;
  double *sigma, *theta;
  double dx;
    
  dim = distr->dim;
  
  sigma = DISTR.param_vecs[INDEX_SIGMA];
  theta = DISTR.param_vecs[INDEX_THETA];
    
  /* calculation of the dlogpdf components */
  for (i=0; i<dim; i++) {
    result[i] = 0.;
    
    for (j=0; j<dim; j++) {
      dx = 0; 
      if (i==j)   dx=1;
      if (i==j-1) dx=-1; 
      /* sigma[j] is expected to be > UNUR_EPSILON here */
      if (sigma) dx /= sigma[j];
      result[i] -= (dim-j) * dx;
    }
    
    result[i] += LOGNORMCONSTANT;
  }
  
  return UNUR_SUCCESS; 

} /* end of _unur_dlogpdf_multiexponential() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_multiexponential( UNUR_DISTR *distr, const double *sigma, const double *theta )
{
  int i;
  double *default_sigma=NULL;
  double *default_theta=NULL;
  
  /* copy the sigma parameters into their parameter vectors */
  if(sigma==NULL) {
    /* initializing vector with default values */
    default_sigma = _unur_xmalloc( distr->dim * sizeof(double));
    for (i=0; i<distr->dim; i++) default_sigma[i]=1.;
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_SIGMA, default_sigma, distr->dim );
  }
  else {
    /* check parameter sigma */
    for (i=0; i<distr->dim; i++) {
      if ( sigma[i] <= UNUR_EPSILON ) {
        _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma is too low");
        return UNUR_ERR_DISTR_DOMAIN;
      }
    }
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_SIGMA, sigma, distr->dim );
  }
    
  /* copy the theta parameters into their parameter vectors */
  if(theta==NULL) {
    /* initializing vector with default values */
    default_theta = _unur_xmalloc(distr->dim * sizeof(double) );
    for (i=0; i<distr->dim; i++) default_theta[i]=0.;
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_THETA, default_theta, distr->dim );
  }
  else {
    unur_distr_cvec_set_pdfparams_vec( distr, INDEX_THETA, theta, distr->dim ); 
  }

  /* store number of parameters */
  DISTR.n_params = 0.; /* we have only vector parameter here ... */

  return UNUR_SUCCESS;
} /* end of _unur_set_params_multiexponential() */


/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multiexponential( int dim, const double *sigma, const double *theta )
/*
   @code{UNUR_DISTR *unur_distr_multiexponential(int dim, const double *sigma, const double *theta)}
   creates a distribution object for the linear multivariate exponential distribution 
   with @var{dim} components. 
   @var{sigma} (the scale parameter) is an array of size @var{dim} and each element must be strictly > 0.
   @var{theta} (the location parameter) is an array of size @var{dim}.
   A NULL pointer for @var{sigma} OR @var{theta} chooses the standard form, where
   the scale parameter is (1,@dots{},1) and the location parameter is (0,@dots,0).
*/   
{
  struct unur_distr *distr;
  struct unur_distr *stdmarginal;
  
  int i;
  double sumsigma; /* used in the calculation of LOGNORMCONSTANT */
    
  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 1 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_MEXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;

  /* functions */
  DISTR.pdf     = _unur_pdf_multiexponential;       /* pointer to PDF */
  DISTR.logpdf  = _unur_logpdf_multiexponential;    /* pointer to logPDF */
  DISTR.dpdf    = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  /* pointer to gradient of PDF */
  DISTR.dlogpdf = _unur_dlogpdf_multiexponential;    /* pointer to gradient of logPDF */
  DISTR.pdpdf    = _unur_distr_cvec_eval_pdpdf_from_pdlogpdf;  /* pointer to part. deriv. of PDF */
  /*DISTR.pdlogpdf = _unur_pdlogpdf_multiexponential;*/  /* pointer to partial derivative of logPDF */

  /* set standardized marginal distributions */  
  stdmarginal = unur_distr_exponential(NULL, 0); 
  unur_distr_cvec_set_stdmarginals(distr, stdmarginal);
  unur_distr_free(stdmarginal);

  /* set parameters for distribution */
  if (_unur_set_params_multiexponential(distr, sigma, theta)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  
  /* domain */

  /* log of normalization constant */
  /*  constant:  Prod_{i=0}^{i=dim-1} 1/sigma_i */
  
  sumsigma = 0.;
  for (i=0; i<distr->dim; i++) {
    sumsigma += DISTR.param_vecs[INDEX_SIGMA][i];
  }
  LOGNORMCONSTANT = - 1. / sumsigma;   
  
  /* mode */
  DISTR.mode = _unur_xmalloc(distr->dim * sizeof(double) );
  for (i=0; i<distr->dim; i++)  DISTR.mode[i]=0.;
  
  /* volume below p.d.f. */
  DISTR.volume = 1.; 

  /* indicate which parameters are set (additional to mean and covariance) */
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_PDFVOLUME |
		  UNUR_DISTR_SET_MODE );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_multiexponential() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
