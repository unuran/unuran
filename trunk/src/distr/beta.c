/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      beta.c                                                       *
 *                                                                           *
 *   Normalization constants for pdf OMITTED!                                *
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
 *  pdf:     f(x) = (x-a)^(p-1) * (b-x)^(q-1)                                *
 *  domain:  a < x < b                                                       *
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
 *  pdf:     f(x) = x^(p-1) * (1-x)^(q-1)                                    *
 *  domain:  0 < x < 1                                                       *
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

#include <unur_distr.h>

#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static char distr_name[] = "Beta distribution";

#define p (param[0])
#define q (param[1])
#define a (param[2])
#define b (param[3])
/*---------------------------------------------------------------------------*/

double
unur_pdf_beta(double x, double *param, int n_param)
{ 
  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 4:  /* non standard */
#if CHECKARGS
    if (a >= b) {
      _unur_error(distr_name,UNUR_ERR_DISTR,"invalid domain: a >= b!");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-a) / (b-a);

  case 2:  /* standard */
#if CHECKARGS
    if (p <= 0. || q <= 0.) {
      _unur_error(distr_name,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
      return 0.;
    }
#endif
    if (x <= 0. || x >= 1.)
      return 0.;
    
    if (p == 1.)
      return pow(1.-x,q-1.);
    
    if (q == 1.)
      return pow(x,p-1.);
    
    return (pow(x,p-1.) * pow(1.-x,q-1.));
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_beta() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_beta(double x, double *param, int n_param)
{ 
  register double factor = 1.;

  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 4:  /* non standard */
#if CHECKARGS
    if (a >= b) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"invalid domain: a >= b!");
      return 0.;
    }
#endif
    /* standardize */
    factor = 1./(b-a);
    x = (x-a) / (b-a);

  case 2:  /* standard */
#if CHECKARGS
    if (p <= 0. || q <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
      return 0.;
    }
#endif
    if (x <= 0. || x >= 1.)
      return 0.;
    
    if (p == 1.)
      return ( -(q-1.) * pow(1.-x,q-2.) * factor );

    if (q == 1.)
      return ((p-1.) * pow(x,p-2.) * factor);

    return (pow(x,p-2.) * pow(1.-x,q-2.) * ( (p-1.)*(1.-x) - (q-1.)*x ) * factor);
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_beta() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_beta(double x, double *param, int n_param)
{
  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 4:  /* non standard */
#if CHECKARGS
    if (a >= b) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"invalid domain: a >= b!");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-a) / (b-a);

  case 2:  /* standard */
#if CHECKARGS
    if (p <= 0. || q <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
      return 0.;
    }
#endif
    if (x <= 0.) return 0.;
    if (x >= 1.) return 1.;

    return _unur_cdf_beta_ext(x,p,q);

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_beta() */

/*---------------------------------------------------------------------------*/

double
unur_mode_beta(double *param, int n_param)
{ 
  double mode = NOT_UNIMODAL;

  CHECK_NULL(param,RETURN_NULL);

#if CHECKARGS
  if (p <= 0. || q <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
    return 0.;
  }
#endif

  if (p <= 1. && q > 1.)
    mode = 0.;              /* left limit of domain */

  if (p > 1. && q <= 1.)
    mode = 1.;              /* right limit of domain */

  if (p > 1. && q > 1.)
    mode = (p - 1.) / (p + q - 2.);

  switch (n_param) {

  case 2:  /* standard */
    return mode;

  case 4:  /* non standard */
#if CHECKARGS
    if (a >= b) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"invalid domain: a >= b!");
      return 0.;
    }
#endif
    if (mode <= 1.) 
      mode = mode * (b-a) + a;
    /* else p.d.f. is not unimodal */

    return mode;

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }
  
} /* end of unur_mode_beta() */

/*---------------------------------------------------------------------------*/

double
unur_area_beta(double *param, int n_param)
{ 

  CHECK_NULL(param,RETURN_NULL);

#if CHECKARGS
    if (p <= 0. || q <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"p <= 0 or q <= 0.");
      return 0.;
    }
#endif

  switch (n_param) {

  case 2:  /* standard */
    /* Beta(p,q) */
    return exp(_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q));

  case 4:  /* non standard */
#if CHECKARGS
    if (a >= b) {
      _unur_error("Beta distribution",UNUR_ERR_DISTR,"invalid domain: a >= b!");
      return 0.;
    }
#endif
    /* Beta(p,q) * (b-a)^(p+q-1) */
    return exp(_unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q) + (b-a)*(p+q-1.) );

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_beta() */

/*---------------------------------------------------------------------------*/
#undef p
#undef q
#undef a
#undef b
/*---------------------------------------------------------------------------*/
