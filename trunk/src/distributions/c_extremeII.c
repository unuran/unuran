/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_extremeII.c                                                *
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
 *  Extreme value type II distribution  [3; ch.22, p.2]                      *
 *  (also Frechet-type distribution)                                         *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *  Type II (also Frechet-type distribution)                                 *
 *                                                                           *
 *  cdf:       F(x) = exp( -((x-zeta)/theta)^(-k) )                          *
 *  pdf:       f(x) = exp( -((x-zeta)/theta)^(-k)) * ((x-zeta)/theta)^(-k-1) *
 *  domain:    zeta < x <infinity                                            *
 *  constant:  k/theta                                                       *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  k     > 0   ... shape                                             *
 *     1:  zeta        ... location                                          *
 *     2:  theta > 0   ... scale                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Standard Form:                                                            *
 *                                                                           *
 *  cdf:       F(x) = exp( -x^(-k) )                                         *
 *  pdf:       f(x) = exp( -x^(-k) ) * x^(-k-1)                              *
 *  domain:    0 < x <infinity                                               *
 *  constant:  k                                                             *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  k           ... shape                                             *
 *                                                                           *
 *     1:  zeta  = 0                                                         *
 *     2:  theta = 1                                                         *
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

static const char distr_name[] = "extremeII";

/* parameters */
#define k      params[0]    /* shape */
#define zeta   params[1]    /* location */
#define theta  params[2]    /* scale */

/* function prototypes                                                       */
static double _unur_pdf_extremeII(double x, double *params, int n_params);
static double _unur_dpdf_extremeII(double x, double *params, int n_params);
static double _unur_cdf_extremeII(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_extremeII( double x, double *params, int n_params )
{ 
  register double xk;

  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x - zeta) / theta;
  case 1: default: /* standard */
    if (x<=0.) return 0.;

    xk = pow( x, -k - 1.);
    return ( exp( -xk * x) * xk * k / theta );
  }
} /* end of _unur_pdf_extremeII() */

/*---------------------------------------------------------------------------*/
double
_unur_dpdf_extremeII( double x, double *params, int n_params )
{ 
  register double factor = 1.;
  register double xk;

  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    factor = 1. / theta;
    x = (x - zeta) / theta;
  case 1: default: /* standard */
    if (x<=0.) return 0.;

    xk = pow( x, -k);
    return ( -exp(xk) * k * pow( x, -2.*(k+1.) ) * (xk + k*(xk-1.)) * factor );
  }
} /* end of unur_dpdf_extremeII() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_extremeII( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x - zeta) / theta;
  case 1: default: /* standard */
    if (x<=0.) return 0.;

    return ( exp( -pow( x, -k ) ) );
  }
} /* end of _unur_cdf_extremeII() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_extremeII( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_EXTREME_II;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_extremeII_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_extremeII;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_extremeII; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_extremeII;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;
  
  /* copy parameters */
  DISTR.k = k;
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
  case 2:
    DISTR.zeta = zeta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameters */
  if (DISTR.k <= 0 || DISTR.theta <= 0. ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0 || theta <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant: not required */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.zeta + pow( k/(k+1.), 1/k ) * DISTR.theta;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.zeta;      /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_extremeII() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
/*---------------------------------------------------------------------------*/
