/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_powerexponential.c                                         *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  distr: power-exponential distribution                                    *
 *    (also called Subbotin distribution)                                    *
 *                                                                           *
 *  (see also [3; ch.24, p.195] for a different definition;                  *
 *  also called Subbotin distribution)                                       *
 *                                                                           *
 *  pdf:       exp(-abs(x)^tau)                                              *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  1 / (2 * Gamma(1+1/tau))                                      *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  tau > 0 ... shape                                                 *
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

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "powerexponential";

/*---------------------------------------------------------------------------*/
/* parameters */
#define tau  params[0]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#ifdef HAVE_UNUR_SF_INCOMPLETE_GAMMA
#  define HAVE_CDF
#else
#  undef  HAVE_CDF
#endif

/* can we compute the area below the pdf ? */
#ifdef HAVE_UNUR_SF_LN_GAMMA
#  define HAVE_AREA
#else
#  undef  HAVE_AREA
#endif

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_powerexponential( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_powerexponential( double x, const UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_powerexponential( double x, const UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_powerexponential( UNUR_DISTR *distr );
#ifdef HAVE_AREA
static int _unur_upd_area_powerexponential( UNUR_DISTR *distr );
#endif
static int _unur_set_params_powerexponential( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_powerexponential( double x, const UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  return exp( - pow( fabs(x), tau ) - LOGNORMCONSTANT);
} /* end of _unur_pdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_powerexponential( double x, const UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  register double tmp;

  if (x == 0.)    /* derivative is not defined, but ...        */
    return 0.;    /* a tangent parallel to x-axis is possible. */

  tmp = exp( -pow(fabs(x),tau) - LOGNORMCONSTANT + (tau-1.)*log(fabs(x)) ) * tau;

  /* sign ! */
  return ( (x<0.) ? tmp : -tmp );
} /* end of _unur_dpdf_powerexponential() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_powerexponential( double x, const UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double cdf;

  /* compute cdf(abs(x)) - cdf(0) */
  cdf = _unur_sf_incomplete_gamma(pow(fabs(x),tau),1./tau) / 2.;
  return ((x<0.) ? 0.5 - cdf : 0.5 + cdf);

} /* end of _unur_cdf_powerexponential() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_powerexponential( UNUR_DISTR *distr )
{
  DISTR.mode = 0;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_powerexponential() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_AREA

int
_unur_upd_area_powerexponential( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_sf_ln_gamma(1. + 1./DISTR.tau) + M_LN2;
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return UNUR_SUCCESS;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_powerexponential( DISTR.domain[1],distr) 
		 - _unur_cdf_powerexponential( DISTR.domain[0],distr) );
  return UNUR_SUCCESS;
#else
  return UNUR_ERR_DISTR_REQUIRED;
#endif

} /* end of _unur_upd_area_powerexponential() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_set_params_powerexponential( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return UNUR_ERR_DISTR_NPARAMS; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,UNUR_ERR_NULL);

  /* check parameter tau */
  if (tau <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"tau <= 0");
    return UNUR_ERR_DISTR_DOMAIN;
  }

  /* copy parameters for standard form */
  DISTR.tau = tau;

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return UNUR_SUCCESS;
} /* end of _unur_set_params_powerexponential() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_powerexponential( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_POWEREXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_powerexponential_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_powerexponential;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_powerexponential; /* pointer to derivative of PDF */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_powerexponential;  /* pointer to CDF               */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );

  /* set parameters for distribution */
  if (_unur_set_params_powerexponential(distr,params,n_params)!=UNUR_SUCCESS) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
#ifdef HAVE_AREA
  LOGNORMCONSTANT = _unur_sf_ln_gamma(1. + 1./DISTR.tau) + M_LN2;
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and area below p.d.f. */
  DISTR.mode = 0;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_powerexponential;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_powerexponential; /* funct for computing mode */
#ifdef HAVE_AREA
  DISTR.upd_area  = _unur_upd_area_powerexponential; /* funct for computing area */
#endif

  /* return pointer to object */
  return distr;

} /* end of unur_distr_powerexponential() */

/*---------------------------------------------------------------------------*/
#undef tau
#undef DISTR
/*---------------------------------------------------------------------------*/
