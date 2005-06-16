/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multistudent.c                                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Multivariate Student distribution [5; ch.45, p.219]               *
 *                                                                           *
 *  pdf:       f(x) = 1 / ( 1 + (x-mu)^t . Sigma^-1 . (x-mu) / m)^(dim+m)/2 )* 
 *  domain:    Reals^(dim)                                                   *
 *  constant:  Gamma((dim+m)/2)                                              *
 *          / ( Gamma(m/2) (m*pi)^(dim/2) * sqrt(det(Sigma)) )               *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  d.f.    ... m       (default : 1)                                 *
 *     1:  mean    ... mu      (default : 0-vector)                          *
 *     2:  "covar" ... Sigma   (default : identity matrix)                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = 1 / ( 1 + x^t . x )^(dim+m)/2                          *
 *  domain:    Reals^(dim)                                                   *
 *                                                                           *
 *  parameters:                                                              *
 *                                                                           *
 *     d.f.  = m                                                             *
 *     mean  = (0,...,0)  ... 0-vector                                       *
 *     covar = identity matrix                                               *
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
#include <utils/matrix_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "multistudent";

/*---------------------------------------------------------------------------*/
/* parameters */
#define nu  params[0]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

    static double _unur_pdf_multistudent( const double *x, UNUR_DISTR *distr );
    static double _unur_logpdf_multistudent( const double *x, UNUR_DISTR *distr );
    static int _unur_dlogpdf_multistudent( double *result, const double *x, UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_multistudent( const double *x, UNUR_DISTR *distr )
{ 
  return exp( _unur_logpdf_multistudent( x, distr ) );
} /* end of _unur_pdf_multistudent() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_multistudent( const double *x, UNUR_DISTR *distr )
{ 
#define idx(a,b) ((a)*dim+(b))

  int i,j, dim;
  double *mean;
  const double *covar_inv; 
  
  double xx; /* argument used in the evaluation of exp(-xx/2) */
  double cx; /* element of multiplication of covariance matrix and x */
  
  dim = distr->dim;
  
  if (DISTR.mean == NULL && DISTR.covar == NULL) {
    /* standard form */
    xx=0.;
    for (i=0; i<dim; i++) { xx += x[i]*x[i]; }
    return ( - (dim+DISTR.nu)/2. * log(1+xx/DISTR.nu) + LOGNORMCONSTANT);  
  }

  mean = DISTR.mean;

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return INFINITY;

  xx=0.; /* resetting exponential function argument */
  for (i=0; i<dim; i++) {
    cx=0.; 
    /* multiplication of inverse covariance matrix and (x-mean) */
    for (j=0; j<dim; j++) {
      cx += covar_inv[idx(i,j)] * (x[j]-mean[j]);
    }
    xx += (x[i]-mean[i])*cx;
  }
  
  return (- (dim+DISTR.nu)/2. * log(1+xx/DISTR.nu) + LOGNORMCONSTANT);

#undef idx
} /* end of _unur_logpdf_multistudent() */

/*---------------------------------------------------------------------------*/

int
_unur_dlogpdf_multistudent( double *result, const double *x, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))

  int i, j, dim;
  double xx, cx;
  double *mean;
  const double *covar_inv;
    
  dim = distr->dim;
  mean = DISTR.mean;

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return UNUR_FAILURE;

  /* calculating (x-mean)' Sigma^-1 (x-mean) */  
  xx=0.; 
  for (i=0; i<dim; i++) {
    cx=0.; 
    /* multiplication of inverse covariance matrix and (x-mean) */
    for (j=0; j<dim; j++) {
      cx += covar_inv[idx(i,j)] * (x[j]-mean[j]);
    }
    xx += (x[i]-mean[i])*cx;
  }    
    
  /* calculation of the dlogpdf components */
  for (i=0; i<dim; i++) {
    result[i] = 0.;
    
    for (j=0; j<dim; j++) 
      result[i] -= (x[j]-mean[j]) * (covar_inv[idx(i,j)]+covar_inv[idx(j,i)]);
  
    result[i] *= .5*(dim+DISTR.nu)/(DISTR.nu+xx);
  }
  
  return UNUR_SUCCESS; 

#undef idx
} /* end of _unur_dlogpdf_multistudent() */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
_unur_set_params_multistudent( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter m (degrees of freedom) */
  if ( nu <= 0 ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.nu = nu; 

  /* store number of parameters */
  DISTR.n_params = n_params;

  return UNUR_SUCCESS;
} /* end of _unur_set_params_multistudent() */


/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multistudent( int dim, const double df, const double *mean, const double *covar )
/*
   @code{UNUR_DISTR *unur_distr_multistudent(int dim, const double nu, const double *mean, const double *covar)}
   creates a distribution object for the multivariate Student t-distribution with
   @var{dim} components and @var{nu} degrees of freedom. 
   @var{mean} is an array of size @var{dim}.
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
   @var{covar} is an array of size @var{dim}x@var{dim} and holds the
   covariance matrix, where the rows of the matrix are stored
   consecutively in this array. The NULL pointer can be used
   instead the identity matrix.
   If @var{covar} is not a valid covariance matrix (i.e., not positive
   definite) then no distribution object is created and NULL is returned.

   For standard form of the distribution use the null vector for @var{mean} and 
   the identity matrix for @var{covar}.
*/   
{
  struct unur_distr *distr;
  struct unur_distr *stdmarginal;
  double det_covar; /* determinant of covariance matrix */
  double par;
  
  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 1 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_MSTUDENT;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;

  /* copy (and check) parameters */
  if ((unur_distr_cvec_set_mean(distr,mean)!=UNUR_SUCCESS) ||
      (unur_distr_cvec_set_covar(distr,covar)!=UNUR_SUCCESS) ) {
    unur_distr_free( distr );
    return NULL;
  }

  /* functions */
  DISTR.pdf     = _unur_pdf_multistudent;       /* pointer to PDF */
  DISTR.logpdf  = _unur_logpdf_multistudent;    /* pointer to logPDF */
  DISTR.dpdf    = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  /* pointer to derivative of PDF */
  DISTR.dlogpdf = _unur_dlogpdf_multistudent;    /* pointer to derivative of logPDF */

  par = df; /* degrees of freedom */
  
  /* set standardized marginal distributions */
  stdmarginal = unur_distr_student(&par,1);
  unur_distr_cvec_set_stdmarginals(distr,stdmarginal);
  unur_distr_free(stdmarginal);

  /* set parameters for distribution */
  if (_unur_set_params_multistudent(distr, &par, 1)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }
  
  /* domain */

  /* log of normalization constant */
  /* constant:  Gamma((dim+1)/2) / ( pi^((dim+1)/2) * sqrt(det(Sigma)) )  */
  
 /*  constant:  Gamma((dim+df)/2)                                          */
 /*          / ( Gamma(df/2) (df*pi)^(dim/2) * sqrt(det(Sigma)) )          */
  det_covar = (DISTR.covar == NULL) ? 1. : _unur_matrix_determinant(dim, DISTR.covar);
  LOGNORMCONSTANT = _unur_sf_ln_gamma((distr->dim+df)/2.) - _unur_sf_ln_gamma(df/2.)
                  - ( distr->dim * log(df*M_PI) + log(det_covar) ) / 2.;

  /* mode */
  DISTR.mode = _unur_xmalloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );

  /* volume below p.d.f. */
  DISTR.volume = 1.; 

  /* indicate which parameters are set (additional to mean and covariance) */
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_PDFVOLUME |
		  UNUR_DISTR_SET_MODE );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_multistudent() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
