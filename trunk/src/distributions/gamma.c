/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gamma.c                                                      *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [2] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 1, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1994                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *  
 *  Gamma distribution [2; ch.17, p.337]                                     *
 *                                                                           *
 *  pdf:       f(x) = (x-gamma)^(alpha-1) * exp( -(x-gamma)/beta )           *
 *  domain:    x > gamma                                                     *
 *  constant:  beta^alpha * Gamma(alpha)                                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0  ... shape                                              *
 *     1:  beta > 0   ... scale                                              *
 *     2:  gamma      ... location                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(alpha-1) * exp(-x)                                  *
 *  domain:    x > 0                                                         *
 *  constant:  Gamma(alpha)                                                  *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0  ... shape                                              *
 *                                                                           *
 *     1:  beta  = 1                                                         *
 *     2:  gamma = 0                                                         *
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
#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "gamma";

/* parameters */
#define alpha  params[0]
#define beta   params[1]
#define gamma  params[2]

/* function prototypes                                                       */
static double _unur_pdf_gamma(double x, double *params, int n_params);
static double _unur_dpdf_gamma(double x, double *params, int n_params);
static double _unur_cdf_gamma(double x, double *params, int n_params);
static double _unur_mode_gamma(double *params, int n_params);
static double _unur_lognormconstant_gamma(double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gamma( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 3:  /* non standard */
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (alpha == 1. && x >= 0.)
      return exp( -x );
    if (x <= 0.)
      return 0.;
    return exp( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);
    /*    return ( pow(x,alpha-1.) * exp(-x) ); */
  }
} /* end of _unur_pdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gamma( double x, double *params, int n_params )
{
  register double factor = 1.;

  switch (n_params) {
  case 3:  /* non standard */
    factor = 1./beta;
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (x <= 0.)
      return 0.;
    if (alpha == 1.)
      return( -exp(-x) * factor );
    return ( pow(x,alpha-2.) * exp(-x - LOGNORMCONSTANT) *  ((alpha-1.) -x) * factor ); 
  }
} /* end of _unur_dpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_gamma( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 3:  /* non standard */
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (x <= 0.)
      return 0.;

    return _unur_incgamma(x,alpha);
  }
} /* end of _unur_cdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_mode_gamma( double *params, int n_params )
{
  register double mode;

  mode = (alpha >= 1.) ? (alpha - 1.) : 0.;

  switch (n_params) {
  case 3:  /* non standard */
    return (mode * beta) + gamma;
  case 1: default: /* standard */
    return mode;
  }
} /* end of _unur_mode_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_gamma( double *params, int n_params )
{
  switch (n_params) {
  case 3:  /* non standard */
    return ( _unur_gammaln(alpha) + log(beta)*alpha );
  case 1: default: /* standard */
    return (_unur_gammaln(alpha));
  }
} /* end of _unur_lognormconstant_gamma() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gamma( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params < 1 || n_params > 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
    return NULL;
  }

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GAMMA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_gamma_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_gamma;    /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_gamma;   /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_gamma;    /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.beta  = 1.;
  DISTR.gamma = 0.;

  /* copy parameters */
  DISTR.alpha = alpha;
  switch (n_params) {
  case 3:
    DISTR.gamma = gamma;
  case 2:
    DISTR.beta = beta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameters alpha and beta */
  if (DISTR.alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= 0.");
    free( distr ); return NULL;
  }
  if (DISTR.beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  DISTR.mode = _unur_mode_gamma(DISTR.params,DISTR.n_params);
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.gamma;  /* left boundary  */
  DISTR.domain[1] = INFINITY;         /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_gamma() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
/*---------------------------------------------------------------------------*/
