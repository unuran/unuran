/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      beta.c                                                       *
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
 *  Beta distribution [3; ch.25, p.210]                                      *
 *                                                                           *
 *  pdf:       f(x) = (x-a)^(p-1) * (b-x)^(q-1)                              *
 *  domain:    a < x < b                                                     *
 *  constant:  Beta(p,q) * (b-a)^(p+q-1)                                     *
 *                                                                           *
 *  parameters: 4                                                            *
 *     0:  p > 0    ... shape                                                *
 *     1:  q > 0    ... shape                                                *
 *     2:  a        ... location                                             *
 *     3:  b (>a)   ... location                                             *
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

#include <source_unuran.h>
#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "beta";

/* parameters */
#define p  params[0]
#define q  params[1]
#define a  params[2]
#define b  params[3]

/* function prototypes                                                       */
static double _unur_pdf_beta(double x, double *params, int n_params);
static double _unur_dpdf_beta(double x, double *params, int n_params);
static double _unur_cdf_beta(double x, double *params, int n_params);
inline static double _unur_mode_beta(double *params, int n_params);
inline static double _unur_lognormconstant_beta(double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_beta(double x, double *params, int n_params)
{ 
  switch (n_params) {
  case 4:                /* non standard */
    x = (x-a) / (b-a);   /* -> standardize */
  case 2: default:       /* standard */
    if (x <= 0. || x >= 1.)
      return 0.;         /* out of support */
    return exp((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);
  }
} /* end of _unur_pdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_beta(double x, double *params, int n_params)
{ 
  register double factor = 1.;

  switch (n_params) {
  case 4:                /* non standard */
    x = (x-a) / (b-a);   /* -> standardize */
    factor = 1./(b-a);
  case 2: default:       /* standard */
    if (x <= 0. || x >= 1.)
      return 0.;         /* out of support */

    return (exp((p-2.)*log(x) + (q-2.)*log(1.-x) - LOGNORMCONSTANT) * ( (p-1.)*(1.-x) - (q-1.)*x ) * factor );
  }
} /* end of _unur_dpdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_beta(double x, double *params, int n_params)
{
  switch (n_params) {
  case 4:                /* non standard */
    x = (x-a) / (b-a);   /* -> standardize */
  case 2: default:       /* standard */
    /* out of support of p.d.f.? */
    if (x <= 0.) return 0.;
    if (x >= 1.) return 1.;

    return _unur_incbeta(x,p,q);
  }
} /* end of _unur_cdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_mode_beta(double *params, int n_params)
{ 
  double mode = INFINITY;

  if (p <= 1. && q > 1.)
    mode = 0.;              /* left limit of domain */

  else if (p > 1. && q <= 1.)
    mode = 1.;              /* right limit of domain */

  else if (p > 1. && q > 1.)
    mode = (p - 1.) / (p + q - 2.);
  
  /* else p.d.f. is not unimodal */

  return( (n_params==2) ? mode : mode * (b-a) + a );
  /** TODO: possible overflow if mode == INFINITY **/
    
} /* end of _unur_mode_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_beta(double *params, int n_params)
{ 
  switch (n_params) {
  case 4:                /* non standard */
    /* log( Beta(p,q) * (b-a)^(p+q-1) ) */
    return (_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q) + (b-a)*(p+q-1.) );
  case 2: default:       /* standard */
    /* log( Beta(p,q) ) */
    return (_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q));
  }
} /* end of _unur_lognormconstant_beta() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_beta( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params == 3)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
  if (n_params > 4)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_BETA;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_beta_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_beta;    /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_beta;   /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_beta;    /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.a = 0.;           /* default for a */
  DISTR.b = 1.;           /* default for b */

  /* copy parameters */
  DISTR.p = p;
  DISTR.q = q;
  if (n_params == 4) {
    DISTR.a = a;
    DISTR.b = b;
  }

  /* check parameters p and q */
  if (DISTR.p <= 0. || DISTR.q <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 or q <= 0");
    free( distr ); return NULL;
  }
  /* check parameters a and b */
  if (DISTR.a >= DISTR.b) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT = _unur_lognormconstant_beta(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  DISTR.mode = _unur_mode_beta(DISTR.params,DISTR.n_params);
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.a; /* left boundary  */
  DISTR.domain[1] = DISTR.b; /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_beta() */

/*---------------------------------------------------------------------------*/
#undef p
#undef q
#undef a
#undef b
/*---------------------------------------------------------------------------*/
