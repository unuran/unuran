/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gamma.c                                                    *
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
 *  distr: Gamma distribution [2; ch.17, p.337]                              *
 *                                                                           *
 *  pdf:       f(x) = ((x-gamma)/beta)^(alpha-1) * exp( -(x-gamma)/beta )    *
 *  domain:    x > gamma                                                     *
 *  constant:  beta * Gamma(alpha)                                           *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0       ... shape                                         *
 *     1:  beta > 0   (1)  ... scale                                         *
 *     2:  gamma      (0)  ... location                                      *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "gamma";

/*---------------------------------------------------------------------------*/
/* parameters */
#define alpha  params[0]   /* shape */
#define beta   params[1]   /* scale */
#define gamma  params[2]   /* location */

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
static double _unur_pdf_gamma(double x, UNUR_DISTR *distr);
static double _unur_dpdf_gamma(double x, UNUR_DISTR *distr);
#ifdef HAVE_CDF
static double _unur_cdf_gamma(double x, UNUR_DISTR *distr);
#endif

static int _unur_upd_mode_gamma( UNUR_DISTR *distr );
#ifdef HAVE_AREA
static int _unur_upd_area_gamma( UNUR_DISTR *distr );
static double _unur_lognormconstant_gamma(double *params, int n_params);
#endif

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gamma( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
  case 2:
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (alpha == 1. && x >= 0.)
      return exp( -x - LOGNORMCONSTANT);
    if (x <= 0.)
      return 0.;
    return exp( (alpha-1.)*log(x) - x - LOGNORMCONSTANT);
    /*    return ( pow(x,alpha-1.) * exp(-x) ); */
  }
} /* end of _unur_pdf_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gamma( double x, UNUR_DISTR *distr )
{
  register double factor = 1.;
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
  case 2:
    factor = 1./beta;
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (alpha == 1. && x >= 0.)
      return( -exp(-x - LOGNORMCONSTANT) * factor );
    if (x <= 0.)
      return 0.;

    return ( exp( log(x) * (alpha-2.) - x - LOGNORMCONSTANT) *  ((alpha-1.) -x) * factor ); 
  }
} /* end of _unur_dpdf_gamma() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_gamma( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
  case 2:
    x = (x-gamma) / beta;     /* standardize */
  case 1: default: /* standard */
    if (x <= 0.)
      return 0.;

    return _unur_sf_incomplete_gamma(x,alpha);
  }
} /* end of _unur_cdf_gamma() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_gamma( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  DISTR.mode = (alpha >= 1.) ? (alpha - 1.) : 0.;

  if (DISTR.n_params > 1)
    DISTR.mode = DISTR.mode * beta + gamma;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_gamma() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_AREA

int
_unur_upd_area_gamma( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_gamma( DISTR.domain[1],distr) 
		 - _unur_cdf_gamma( DISTR.domain[0],distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_area_gamma() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_gamma( double *params, int n_params )
{
  switch (n_params) {
  case 3:  /* non standard */
  case 2:
    return ( _unur_sf_ln_gamma(alpha) + log(beta) );
  case 1: default: /* standard */
    return (_unur_sf_ln_gamma(alpha));
  }
} /* end of _unur_lognormconstant_gamma() */

#endif

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
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GAMMA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_gamma_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_gamma;    /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_gamma;   /* pointer to derivative of PDF */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_gamma;    /* pointer to CDF               */
#endif

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

  /* domain */
  DISTR.domain[0] = DISTR.gamma;  /* left boundary  */
  DISTR.domain[1] = INFINITY;     /* right boundary */

  /* log of normalization constant */
#ifdef HAVE_AREA
  LOGNORMCONSTANT = _unur_lognormconstant_gamma(DISTR.params,DISTR.n_params);
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and area below p.d.f. */
  _unur_upd_mode_gamma( distr );
  DISTR.area = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_gamma; /* funct for computing mode */
#ifdef HAVE_AREA
  DISTR.upd_area  = _unur_upd_area_gamma; /* funct for computing area */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_gamma() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
#undef DISTR
/*---------------------------------------------------------------------------*/
