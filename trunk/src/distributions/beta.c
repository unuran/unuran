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

#include <float.h>
#include <stdlib.h>

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "beta";

/* parameters */
#define p (params[0])
#define q (params[1])
#define a (params[2])
#define b (params[3])

/*---------------------------------------------------------------------------*/

double
_unur_pdf_beta(double x, double *params, int n_params)
{ 
  switch (n_params) {
  case 4:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 2:  /* standard */
    if (x <= 0. || x >= 1.)
      return 0.;
    else
      return exp((p-1.)*log(x) + (q-1.)*log(1.-x) - LOGNORMCONSTANT);

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return INFINITY;
  }

} /* end of _unur_pdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_beta(double x, double *params, int n_params)
{ 
  register double factor = 1.;

  switch (n_params) {
  case 4:  /* non standard */
    /* standardize */
    factor = 1./(b-a);
    x = (x-a) / (b-a);

  case 2:  /* standard */
    if (x <= 0. || x >= 1.)
      return 0.;
    else
      return (exp((p-2.)*log(x) + (q-2.)*log(1.-x) - LOGNORMCONSTANT) * ( (p-1.)*(1.-x) - (q-1.)*x ) * factor );

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return INFINITY;
  }

} /* end of _unur_dpdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_beta(double x, double *params, int n_params)
{
  switch (n_params) {
  case 4:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 2:  /* standard */
    if (x <= 0.) return 0.;
    if (x >= 1.) return 1.;

    return _unur_cdf_beta_ext(x,p,q);

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return INFINITY;
  }

} /* end of _unur_cdf_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_mode_beta(double *params, int n_params)
{ 
  double mode = INFINITY;

  if (p <= 1. && q > 1.)
    mode = 0.;              /* left limit of domain */

  if (p > 1. && q <= 1.)
    mode = 1.;              /* right limit of domain */

  if (p > 1. && q > 1.)
    mode = (p - 1.) / (p + q - 2.);

  switch (n_params) {

  case 2:  /* standard */
    return mode;

  case 4:  /* non standard */
    if (mode <= 1.) 
      mode = mode * (b-a) + a;
    /* else p.d.f. is not unimodal */

    return mode;

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return INFINITY;
  }
  
} /* end of _unur_mode_beta() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_beta(double *params, int n_params)
{ 
  switch (n_params) {
  case 2:  /* standard */
    /* log( Beta(p,q) ) */
    return (_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q));

  case 4:  /* non standard */
    /* log( Beta(p,q) * (b-a)^(p+q-1) ) */
    return (_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q) + (b-a)*(p+q-1.) );

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return INFINITY;
  }

} /* end of _unur_lognormconstant_beta() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_beta( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params != 2 && n_params != 4) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }

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
  DISTR.params[2] = 0.;           /* default for a */
  DISTR.params[3] = 1.;           /* default for b */

  /* copy parameters */
  DISTR.params[0] = p;
  DISTR.params[1] = q;
  switch (n_params) {
  case 4:
    DISTR.params[3] = b;
  case 3:
    DISTR.params[2] = a;
    n_params = 4;              /* number of parameters for non-standard form */
  default:
  }

  /* check parameters p and q */
  if (DISTR.params[0] <= 0. || DISTR.params[1] <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
    free( distr ); return NULL;
  }
  /* check parameters a and b */
  if (DISTR.params[2] >= DISTR.params[3]) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"invalid domain: a >= b!");
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
  DISTR.domain[0] = DISTR.params[2]; /* left boundary  */
  DISTR.domain[1] = DISTR.params[3]; /* right boundary */

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
