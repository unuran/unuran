/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_beta.c                                                     *
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
 *  distr: Beta distribution [3; ch.25, p.210]                                *
 *                                                                           *
 *  pdf:       f(x) = (x-a)^(p-1) * (b-x)^(q-1)                              *
 *  domain:    a < x < b                                                     *
 *  constant:  Beta(p,q) * (b-a)^(p+q-1)                                     *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0:  p > 0        ... shape                                            *
 *     1:  q > 0        ... shape                                            *
 *     2:  a      (0)   ... location                                         *
 *     3:  b >a   (1)   ... location                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(p-1) * (1-x)^(q-1)                                  *
 *  domain:    0 < x < 1                                                     *
 *  constant:  Beta(p,q)                                                     *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  p > 0    ... shape                                                *
 *     1:  q > 0    ... shape                                                *
 *                                                                           *
 *     2:  a = 0                                                             *
 *     3:  b = 1                                                             *
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

static const char distr_name[] = "beta";

/*---------------------------------------------------------------------------*/
/* parameters */
#define p  params[0]
#define q  params[1]
#define a  params[2]
#define b  params[3]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#ifdef HAVE_UNUR_SF_INCOMPLETE_BETA
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
static double _unur_pdf_beta( double x, UNUR_DISTR *distr );
static double _unur_dpdf_beta( double x, UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_beta( double x, UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_beta( UNUR_DISTR *distr );
#ifdef HAVE_AREA
static int _unur_upd_area_beta( UNUR_DISTR *distr );
inline static double _unur_lognormconstant_beta( double *params, int n_params );
#endif
static int _unur_set_params_beta( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_beta(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (DISTR.n_params > 2)
    /* standardize */
    x = (x-a) / (b-a);

  /* standard form */

  if (x > 0. && x < 1.)
    return exp((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);

  if (x < 0. || x > 1.)
    /* out of support */
    return 0.;

  if ( (x==0. && p==1.) || (x==1. && q==1.) )
    return exp(-LOGNORMCONSTANT);

  /* else */
  return 0.;

} /* end of _unur_pdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_beta(double x, UNUR_DISTR *distr)
{ 
  register double factor = 1.;
  register double *params = DISTR.params;

  if (DISTR.n_params > 2) {
    /* standardize */
    factor = 1./(b-a);
    x = (x-a) / (b-a);
  }

  /* standard form */

  if (x > 0. && x < 1.)
    return (exp((p-2.)*log(x) + (q-2.)*log(1.-x) - LOGNORMCONSTANT) * ( (p-1.)*(1.-x) - (q-1.)*x ) * factor );

  if (x < 0. || x > 1.)
    return 0.;      /* out of support */

  if (x==0. && p==1.)
    return (-(q-1.) * exp(-LOGNORMCONSTANT));

  if (x==1. && q==1.)
    return ((p-1.) * exp(-LOGNORMCONSTANT));

  /* else */
  return 0.;

} /* end of _unur_dpdf_beta() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_beta(double x, UNUR_DISTR *distr)
{
  register double *params = DISTR.params;

  if (DISTR.n_params > 2)
    /* standardize */
    x = (x-a) / (b-a);

  /* standard form */

  /* out of support of p.d.f.? */
  if (x <= 0.) return 0.;
  if (x >= 1.) return 1.;

  return _unur_sf_incomplete_beta(x,p,q);

} /* end of _unur_cdf_beta() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_beta( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  if (p <= 1. && q > 1.)
    DISTR.mode = 0.;              /* left limit of domain */

  else if (p > 1. && q <= 1.)
    DISTR.mode = 1.;              /* right limit of domain */

  else if (p > 1. && q > 1.)
    DISTR.mode = (p - 1.) / (p + q - 2.);
  
  else {
    /* p.d.f. is not unimodal */
    DISTR.mode = INFINITY;
    return 0;
  }

  if (DISTR.n_params > 2)
    DISTR.mode = DISTR.mode * (b - a) + a;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_beta() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_AREA

int
_unur_upd_area_beta( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);
  
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_beta( DISTR.domain[1],distr) 
		 - _unur_cdf_beta( DISTR.domain[0],distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_area_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_beta(double *params, int n_params)
{ 
  if (n_params > 2)
    /* non-standard form */
    /* log( Beta(p,q) * (b-a) ) */
    return (_unur_sf_ln_gamma(p) + _unur_sf_ln_gamma(q) - _unur_sf_ln_gamma(p+q) + log(b-a) );

  else
    /* standard form */
    /* log( Beta(p,q) ) */
    return (_unur_sf_ln_gamma(p) + _unur_sf_ln_gamma(q) - _unur_sf_ln_gamma(p+q));

} /* end of _unur_lognormconstant_beta() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_set_params_beta( UNUR_DISTR *distr, double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params == 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
    n_params = 2; }
  if (n_params > 4) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 4; }
  CHECK_NULL(params,0);


  /* check parameters p and q */
  if (p <= 0. || q <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 or q <= 0");
    return 0;
  }

  /* check parameters a and b */
  if (n_params > 2 && a >= b) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    return 0;
  }

  /* copy parameters for standard form */
  DISTR.p = p;
  DISTR.q = q;

  /* copy optional parameters */
  if (n_params > 2) {
    DISTR.a = a;
    DISTR.b = b;
  }
  else { /* or use defaults */
    DISTR.a = 0.;      /* default for a */
    DISTR.b = 1.;      /* default for b */
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.a; /* left boundary  */
    DISTR.domain[1] = DISTR.b; /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_beta() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_beta( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_BETA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_beta_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_beta;    /* pointer to PDF                  */
  DISTR.dpdf = _unur_dpdf_beta;   /* pointer to derivative of PDF    */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_beta;    /* pointer to CDF                  */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );

  /* set parameters for distribution */
  if (!_unur_set_params_beta(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
#ifdef HAVE_AREA
  LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and area below p.d.f. */
  _unur_upd_mode_beta( distr );
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_beta;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_beta; /* funct for computing mode */
#ifdef HAVE_AREA
  DISTR.upd_area  = _unur_upd_area_beta; /* funct for computing area */
#endif

  /* return pointer to object */
  return distr;

} /* end of unur_distr_beta() */

/*---------------------------------------------------------------------------*/
#undef p
#undef q
#undef a
#undef b
#undef DISTR
/*---------------------------------------------------------------------------*/
