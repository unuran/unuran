/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_burr.c                                                     *
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
 *  Burr family of distributions [2; ch.12, p.54]                            *
 *                                                                           *
 *  first parameter must be type of distribution (1 ... 12)!                 *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type I                                          *
 *  first parameter is 1                                                     *
 *                                                                           *
 *  cdf:       F(x) = x                                                      *
 *  pdf:       f(x) = 1                                                      *
 *  domain:    0 <= x <= 1                                                   *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type II                                         *
 *  first parameter is 2                                                     *
 *                                                                           *
 *  cdf:       F(x) = (exp(-x) + 1)^(-k)                                     *
 *  pdf:       f(x) = k * exp(-x) * (exp(-x)+1)^(-k-1)                       *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  k > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type III                                        *
 *  first parameter is 3                                                     *
 *                                                                           *
 *  cdf:       F(x) = (x^(-c)+1)^(-k)                                        *
 *  pdf:       f(x) =  c * k * x^(-c-1) * (x^(-c) + 1)^(-k-1)                *
 *  domain:    0 < x < infinity                                              *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type IV                                         *
 *  first parameter is 4                                                     *
 *                                                                           *
 *  cdf:       F(x) = (((c-x)/x)^(1/c) + 1)^(-k)                             *
 *  pdf:       f(x) =                                                        *
 *  domain:    0 < x < c                                                     *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type V                                          *
 *  first parameter is 5                                                     *
 *                                                                           *
 *  cdf:       F(x) = (c * exp(-tan(x)) + 1)^(-k)                            *
 *  pdf:       f(x) =                                                        *
 *  domain:    -pi/2 < x < pi/2                                              *
 *  constant:                                                                *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type VI                                         *
 *  first parameter is 6                                                     *
 *                                                                           *
 *  cdf:       F(x) = (c * exp(-k*sinh(x)) + 1)^(-k)                         *
 *  pdf:       f(x) =                                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:                                                                *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type VII                                        *
 *  first parameter is 7                                                     *
 *                                                                           *
 *  cdf:       F(x) = 2^(-k) * (1 + tanh(x))^k                               *
 *  pdf:       f(x) =                                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:                                                                *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type VIII                                       *
 *  first parameter is 8                                                     *
 *                                                                           *
 *  cdf:       F(x) = (2/pi * arctan(exp(x)))^k                              *
 *  pdf:       f(x) =                                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:                                                                *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  k > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type IX                                         *
 *  first parameter is 9                                                     *
 *                                                                           *
 *  cdf:       F(x) = 1 - 2 / (2 + c * ((1+exp(x))^k - 1))                   *
 *  pdf:       f(x) =                                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:                                                                *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type X                                          *
 *  first parameter is 10                                                    *
 *                                                                           *
 *  cdf:       F(x) = (1 - exp(-x^2))^k                                      *
 *  pdf:       f(x) =                                                        *
 *  domain:    0 < x < infinity                                              *
 *  constant:                                                                *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  k > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type XI                                         *
 *  first parameter is 11                                                    *
 *                                                                           *
 *  cdf:       F(x) = (x - (pi/2) * sin( 2*pi*x))^k                          *
 *  pdf:       f(x) =                                                        *
 *  domain:    0 < x < 1                                                     *
 *  constant:                                                                *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0                                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  distr: Burr distribution Type XII                                        *
 *  first parameter is 12                                                    *
 *                                                                           *
 *  cdf:       F(x) = 1 - (1 + x^c)^(-k)                                     *
 *  pdf:       f(x) =                                                        *
 *  domain:    0 < x < infinity                                              *
 *  constant:                                                                *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0                                                             *
 *     1:  c > 0                                                             *
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
static const char distr_name[] = "burr";

/* parameters */
#define burr_type  params[0]
#define k          params[1]
#define c          params[2]

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
/*  static double _unur_pdf_burr(double x, UNUR_DISTR *distr);    */
/*  static double _unur_dpdf_burr(double x, UNUR_DISTR *distr);   */
static double _unur_cdf_burr(double x, UNUR_DISTR *distr);

/*---------------------------------------------------------------------------*/

double
_unur_cdf_burr( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  switch ((int) (burr_type + 0.5)) {

  case  1: /* Type I:   F(x) = x                                             */
    if (x<=0.) return 0.;
    if (x>=1.) return 1.;
    return x;

  case  2: /* Type II:  F(x) = (exp(-x) + 1)^(-k)                            */
    return pow(exp(-x) + 1., -k);

  case  3: /* Type III: F(x) = (x^(-c)+1)^(-k)                               */
    if (x<=0.) return 0.;
    return pow( pow(x, -c) + 1., -k);

  case  4: /* Type IV:  F(x) = (((c-x)/x)^(1/c) + 1)^(-k)                    */
    if (x<=0.) return 0.;
    if (x>=c)  return 1.;
    return pow( pow( (c-x)/x, 1/c ) + 1., -k);

  case  5: /* Type V:   F(x) = (c * exp(-tan(x)) + 1)^(-k)                   */
    if (x<=-M_PI/2.) return 0.;
    if (x>= M_PI/2.) return 1.;
    return pow( c * exp(-tan(x)) + 1., -k );

  case  6: /* Type VI:  F(x) = (c * exp(-k*sinh(x)) + 1)^(-k)                */
    return pow( c * exp(-k*sinh(x)) + 1., -k );

  case  7: /* Type VII: F(x) = 2^(-k) * (1 + tanh(x))^k                      */
    return pow( (1. + tanh(x))/2, k );

  case  8: /* Type VIII:F(x) = (2/pi * arctan(exp(x)))^k                     */
    return pow( 2./M_PI * atan(exp(x)), k );

  case  9: /* Type IX:  F(x) = 1 - 2 / (2 + c * ((1+exp(x))^k - 1))          */
    return (1. - 2. / (2. + c * (pow(1.+exp(x), k) - 1.)));

  case 10: /* Type X:   F(x) = (1 - exp(-x^2))^k                             */
    if (x<=0.) return 0.;
    return pow( 1. - exp(-x*x), k );

  case 11: /* Type XI:  F(x) = (x - (pi/2) * sin( 2*pi*x))^k                 */
    if (x<=0.) return 0.;
    if (x>=1.) return 1.;
    return pow( x - M_PI/2. * sin( 2. * M_PI * x), k );

  case 12: /* Type XII: F(x) = 1 - (1 + x^c)^(-k)                            */
    if (x<=0.) return 0.;
    return (1. - pow( 1 + pow(x, c), -k ) );

  default:
    _unur_error(distr_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0.;
  }

} /* end of _unur_cdf_burr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_burr( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* get type of distribution and 
     set distribution id           */
  switch ((int) (burr_type + 0.5)) {
  case  1:  distr->id = UNUR_DISTR_BURR_I;    break;
  case  2:  distr->id = UNUR_DISTR_BURR_II;   break;
  case  3:  distr->id = UNUR_DISTR_BURR_III;  break;
  case  4:  distr->id = UNUR_DISTR_BURR_IV;   break;
  case  5:  distr->id = UNUR_DISTR_BURR_V;    break;
  case  6:  distr->id = UNUR_DISTR_BURR_VI;   break;
  case  7:  distr->id = UNUR_DISTR_BURR_VII;  break;
  case  8:  distr->id = UNUR_DISTR_BURR_VIII; break;
  case  9:  distr->id = UNUR_DISTR_BURR_IX;   break;
  case 10:  distr->id = UNUR_DISTR_BURR_X;    break;
  case 11:  distr->id = UNUR_DISTR_BURR_XI;   break;
  case 12:  distr->id = UNUR_DISTR_BURR_XII;  break;
  default:
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"type < 1 || type > 12");
    free( distr ); return NULL;
  }

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_burr_init;

  /* functions */
  /* DISTR.pdf  = _unur_pdf_burr;  pointer to PDF                  */
  /* DISTR.dpdf = _unur_dpdf_burr; pointer to derivative of PDF    */
  DISTR.cdf  = _unur_cdf_burr;  /* pointer to CDF                  */

  /* check number of parameters for == 3 */
  switch (distr->id) {
  case UNUR_DISTR_BURR_III:
  case UNUR_DISTR_BURR_IV:
  case UNUR_DISTR_BURR_V:
  case UNUR_DISTR_BURR_VI:
  case UNUR_DISTR_BURR_IX:
  case UNUR_DISTR_BURR_XII:
    if (n_params < 3) {
      _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few");
      free( distr ); return NULL;
    }
  default: /* all other cases */
    if (n_params == 3) {
      _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
      n_params = 2;
    }
  }

  /* copy parameters */
  DISTR.burr_type = burr_type;
  switch (n_params) {
  case 3:
    DISTR.c = c;
  case 2:
    DISTR.k = k;
  default:
  }

  /* check parameters */
  if (DISTR.k <= 0. || (DISTR.c <= 0. && n_params == 3) ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0 || c <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
/*    DISTR.mode = 0.; */
/*    DISTR.area = 1.; */

  /* domain */
  DISTR.domain[0] = -INFINITY;  /* left boundary  */
  DISTR.domain[1] = INFINITY;   /* right boundary */

  switch (distr->id) {
  case UNUR_DISTR_BURR_I:
    DISTR.domain[0] = 0.;       /* left boundary  */
    DISTR.domain[1] = 1.;       /* right boundary */
    break;
  case UNUR_DISTR_BURR_III:
    DISTR.domain[0] = 0.;       /* left boundary  */
    break;
  case UNUR_DISTR_BURR_IV:
    DISTR.domain[0] = 0.;       /* left boundary  */
    DISTR.domain[1] = DISTR.c;  /* right boundary */
    break;
  case UNUR_DISTR_BURR_V:
    DISTR.domain[0] = -M_PI/2.; /* left boundary  */
    DISTR.domain[1] = M_PI/2.;  /* right boundary */
    break;
  case UNUR_DISTR_BURR_X:
    DISTR.domain[0] = 0.;       /* left boundary  */
    break;
  case UNUR_DISTR_BURR_XI:
    DISTR.domain[0] = 0.;       /* left boundary  */
    DISTR.domain[1] = 1.;       /* right boundary */
    break;
  case UNUR_DISTR_BURR_XII:
    DISTR.domain[0] = 0.;       /* left boundary  */
    break;
  }

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN );
		 /* UNUR_DISTR_SET_MODE   | */
		 /* UNUR_DISTR_SET_PDFAREA ); */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_burr() */

/*---------------------------------------------------------------------------*/
#undef burr_type
#undef k
#undef c
#undef DISTR
/*---------------------------------------------------------------------------*/
