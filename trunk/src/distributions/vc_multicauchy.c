/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multicauchy.c                                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  distr: Multivariate Cauchy distribution                                  *
 *                                                                           *
 *  pdf:       f(x) = 1 / ( 1 + (x-mu)^t . Sigma^-1 . (x-mu) )^(dim+1)/2     * 
 *  domain:    Reals^(dim)                                                   *
 *  constant:  (2 pi)^(dim/2) * sqrt(det(Sigma))                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mean    ... mu      (default : 0-vector)                          *
 *     1:  "covar" ... Sigma   (default : identity matrix)                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = 1 / ( 1 + x^t . x )^(dim+1)/2                          *
 *  domain:    Reals^(dim)                                                   *
 *  constant:  (pi)^(dim/2)*c[p)                                             *
 *             with c(p)=c(p-1)/p for odd dim=2p+1, and c(0)=sqrt(pi)        *
 *             and  c(p)=c(p-1)*2/(2p-1) for even dim=2p, and c(1)=2         *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
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

static const char distr_name[] = "multicauchy";

/*---------------------------------------------------------------------------*/
/* parameters */

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

    static double _unur_pdf_multicauchy( const double *x, UNUR_DISTR *distr );
    static double _unur_logpdf_multicauchy( const double *x, UNUR_DISTR *distr );
    static int _unur_dlogpdf_multicauchy( double *result, const double *x, UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_multicauchy( const double *x, UNUR_DISTR *distr )
{ 
  return exp( _unur_logpdf_multicauchy( x, distr ) );
} /* end of _unur_pdf_multicauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_multicauchy( const double *x, UNUR_DISTR *distr )
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
    return ( - (dim+1)/2. * log(1+xx) - LOGNORMCONSTANT);  
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
  
  return (- (dim+1)/2. * log(1+xx) - LOGNORMCONSTANT);

#undef idx
} /* end of _unur_logpdf_multicauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_dlogpdf_multicauchy( double *result, const double *x, UNUR_DISTR *distr )
{
#define idx(a,b) ((a)*dim+(b))

  int i, dim;
  double *mean;
  const double *covar_inv;
    
  dim = distr->dim;
  mean = DISTR.mean;

  /* get inverse of covariance matrix */
  covar_inv = unur_distr_cvec_get_covar_inv(distr);
  if (covar_inv==NULL) 
    /* inverse of covariance matrix not available */
    return UNUR_FAILURE;

  for (i=0; i<dim; i++) {
    result[i] = 0.;
    
#if 0    
    for (j=0; j<dim; j++) 
      /* TODO : */
      result[i] += -0.5 * (x[j]-mean[j]) * (covar_inv[idx(i,j)]+covar_inv[idx(j,i)]);
#endif
  }
  
  return UNUR_SUCCESS; 

#undef idx
} /* end of _unur_dlogpdf_multicauchy() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multicauchy( int dim, const double *mean, const double *covar )
{
  struct unur_distr *distr;
  struct unur_distr *stdmarginal;
  double det_covar; /* determinant of covariance matrix */
  double logc; /* used for LOGNORMCONSTANT calculations */
  int p;
  
  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 1 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_MCAUCHY;

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
  DISTR.pdf     = _unur_pdf_multicauchy;       /* pointer to PDF */
  DISTR.logpdf  = _unur_logpdf_multicauchy;    /* pointer to logPDF */
  DISTR.dpdf    = _unur_distr_cvec_eval_dpdf_from_dlogpdf;  /* pointer to derivative of PDF */
  DISTR.dlogpdf = _unur_dlogpdf_multicauchy;    /* pointer to derivative of logPDF */

  /* set standardized marginal distributions */
  /* TODO : correct stdmarginals */
  stdmarginal = unur_distr_normal(NULL,0);
  unur_distr_cvec_set_stdmarginals(distr,stdmarginal);
  unur_distr_free(stdmarginal);

  /* copy other parameters of distribution */
  /* none */

  /* number of other parameters */
  /* DISTR.n_params = 0;  ... default */

  /* domain */

  /* log of normalization constant */
  if (_unur_FP_equal(distr->dim/2, distr->dim/2.)) {
  /* even dimension */
    logc=log(2); /* for dim=2p=2 i.e. p=1 */
    for (p=2; p<=distr->dim/2; p++) logc += log(2)-log(2*p-1);
  }
  else {
  /* odd dimension */
    logc=log(M_PI)/2; /* for dim=2p+1=1 i.e. p=0 */
    for (p=1; p<=(distr->dim-1)/2; p++) logc -= log(p);
  }
  det_covar = (DISTR.covar == NULL) ? 1. : _unur_matrix_determinant(dim, DISTR.covar);
  LOGNORMCONSTANT = logc + ( distr->dim * log(M_PI) + log(det_covar) ) / 2.;

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

} /* end of unur_distr_multicauchy() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
