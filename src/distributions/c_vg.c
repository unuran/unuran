/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_vg.c                                                       *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  distr: Variance gamma distribution                                       *
 *                                                                           *
 *  pdf:   f(x) = |x-mu|^(lambda-1/2) * exp(beta*(x-mu))                     *
 *                    * K_{lambda-1/2}(alpha*|x-mu|)}                        *
 *                                                                           *
 *  domain:   infinity < x < infinity                                        *
 *                                                                           *
 *  constant: (alpha^2 - beta^2)^lambda                                      *
 *                / (sqrt(pi) * (2*alpha)^(lambda-1/2) * Gamma(lambda))      *
 *                                                                           *
 *             [K_theta(.) ... modified Bessel function of second kind]      *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0 : lambda > 0    ... shape                                           *
 *     1 : alpha >|beta| ... shape                                           *
 *     2 : beta          ... shape (asymmetry)                               *
 *     3 : mu            ... location                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   We use the Rmath library for computing the Bessel function K_n.         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2011 Wolfgang Hoermann and Josef Leydold                  *
 *   Institute for Statistics and Mathematics, WU Wien, Austria              *
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
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "vg";

/* parameters */
#define lambda  params[0]    /* shape */
#define alpha   params[1]    /* shape */
#define beta    params[2]    /* shape (asymmetry) */
#define mu      params[3]    /* location */

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
#ifdef _unur_SF_bessel_k
static double _unur_pdf_vg( double x, const UNUR_DISTR *distr );
static double _unur_logpdf_vg( double x, const UNUR_DISTR *distr );
/* static double _unur_dpdf_vg( double x, const UNUR_DISTR *distr ); */
/* static double _unur_cdf_vg( double x, const UNUR_DISTR *distr ); */
#endif

static int _unur_upd_center_vg( UNUR_DISTR *distr );
static double _unur_lognormconstant_vg( const double *params, int n_params );
static int _unur_set_params_vg( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

#ifdef _unur_SF_bessel_k
double
_unur_pdf_vg(double x, const UNUR_DISTR *distr)
{
  /* Original implementation by Kemal Dingic */
  /* f(x) = |x-mu|^(lambda-1/2) * exp(beta*(x-mu)) * K_{lambda-1/2}(alpha*|x-mu|)} */

  return exp(_unur_logpdf_vg(x,distr));
} /* end of _unur_pdf_vg() */

/*---------------------------------------------------------------------------*/

double
_unur_logpdf_vg(double x, const UNUR_DISTR *distr)
{
  /* Original implementation by Kemal Dingic */
  /* f(x) = |x-mu|^(lambda-1/2) * exp(beta*(x-mu)) * K_{lambda-1/2}(alpha*|x-mu|)} */

  const double *params = DISTR.params;
  double nu = lambda - 0.5;   /* order of modified Bessel function K()       */
  double res;                 /* result of computation                       */
  double y, absy;             /* auxiliary variables                         */

  res = LOGNORMCONSTANT;

  y = x - mu;
  if (_unur_iszero(y)) {
    res += -M_LN2 + _unur_SF_ln_gamma(nu) + nu*log(2./alpha);
  }
  else {
    absy = fabs(y);
    res += log(absy)*nu + beta*y;
    if (lambda < 50) 
      /* threshold value 50 is selected by experiments */
      res += _unur_SF_ln_bessel_k(alpha*absy, nu);
    else
      res += _unur_SF_bessel_k_nuasympt(alpha*absy, nu, TRUE, FALSE);
  }

  return res;
} /* end of _unur_logpdf_vg() */
#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_center_vg( UNUR_DISTR *distr )
{
  const double *params = DISTR.params;

  /* we simply use parameter 'mu' */
  DISTR.center = mu;

  /* an alternative approach would be the mean of the distribution:          */
  /* double gamma = sqrt(alpha*alpha-beta*beta);                             */
  /* DISTR.center = mu + 2*beta*lambda / (gamma*gamma);                      */

  /* center must be in domain */
  if (DISTR.center < DISTR.domain[0])
    DISTR.center = DISTR.domain[0];
  else if (DISTR.center > DISTR.domain[1])
    DISTR.center = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_center_vg() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_vg(const double *params, int n_params ATTRIBUTE__UNUSED)
{
  /*
    (alpha^2 - beta^2)^lambda 
    / (sqrt(pi) * (2*alpha)^(lambda-1/2) * Gamma(lambda))
  */

  return (lambda*log(alpha*alpha - beta*beta) - 0.5*M_LNPI 
	  - (lambda-0.5)*log(2*alpha) - _unur_SF_ln_gamma(lambda));
} /* end of _unur_normconstant_vg() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_vg( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 4) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter omega */
  if (lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  if (alpha <= fabs(beta)) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"alpha <= |beta|");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.lambda = lambda;
  DISTR.alpha = alpha;
  DISTR.beta = beta;
  DISTR.mu = mu;

  /* default parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;   /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_vg() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_vg( const double *params, int n_params)
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_VG;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  /* DISTR.init = _unur_stdgen_vg_init; */

  /* functions */
#ifdef _unur_SF_bessel_k
  DISTR.pdf     = _unur_pdf_vg;     /* pointer to PDF                  */
  DISTR.logpdf  = _unur_logpdf_vg;  /* pointer to log-PDF              */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_CENTER |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* set parameters for distribution */
  if (_unur_set_params_vg(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_vg(params,n_params);

  /* we need the center of the distribution */
  if (_unur_upd_center_vg(distr)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* mode and area below p.d.f. */
  /* DISTR.mode = ? */
  DISTR.area = 1;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_vg;

  /* function for updating derived parameters */
  /* DISTR.upd_mode  = _unur_upd_mode_vg; /\* funct for computing mode *\/ */
  /* DISTR.upd_area  = _unur_upd_area_vg; /\* funct for computing area *\/ */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_vg() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef alpha
#undef beta
#undef lambda
#undef DISTR
/*---------------------------------------------------------------------------*/
