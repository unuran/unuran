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
 *  distr: Extreme value type II distribution  [3; ch.22, p.2]               *
 *  (also Frechet-type distribution)                                         *
 *                                                                           *
 *                                                                           *
 *  Type II (also Frechet-type distribution)                                 *
 *                                                                           *
 *  cdf:       F(x) = exp( -((x-zeta)/theta)^(-k) )                          *
 *  pdf:       f(x) = exp( -((x-zeta)/theta)^(-k)) * ((x-zeta)/theta)^(-k-1) *
 *  domain:    zeta < x <infinity                                            *
 *  constant:  k / theta                                                     *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  k     > 0        ... shape                                        *
 *     1:  zeta        (0)  ... location                                     *
 *     2:  theta > 0   (1)  ... scale                                        *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "extremeII";

/* parameters */
#define k      params[0]    /* shape */
#define zeta   params[1]    /* location */
#define theta  params[2]    /* scale */

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_extremeII( double x, UNUR_DISTR *distr );
static double _unur_dpdf_extremeII( double x, UNUR_DISTR *distr );
static double _unur_cdf_extremeII( double x, UNUR_DISTR *distr );

static int _unur_upd_mode_extremeII( UNUR_DISTR *distr );
static int _unur_upd_area_extremeII( UNUR_DISTR *distr );
static int _unur_set_params_extremeII( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_extremeII( double x, UNUR_DISTR *distr )
{ 
  register double xk;
  register double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x - zeta) / theta;

  /* standard form */

  if (x<=0.) return 0.;

  xk = pow( x, -k - 1.);
  return ( exp( -xk * x - LOGNORMCONSTANT) * xk * k );

} /* end of _unur_pdf_extremeII() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_extremeII( double x, UNUR_DISTR *distr )
{ 
  register double factor = 1.;
  register double xk;
  register double *params = DISTR.params;

  if (DISTR.n_params > 1) {
    /* standardize */
    factor = 1. / (theta * theta);
    x = (x - zeta) / theta;
  }

  /* standard form */

  if (x<=0.) return 0.;

  xk = pow(x, k);
  return (- factor * exp(-1./xk) * k * (xk + k*(xk - 1.)) / pow(x,2. + 2.*k)) ;

} /* end of unur_dpdf_extremeII() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_extremeII( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  if (DISTR.n_params > 1)
    /* standardize */
    x = (x - zeta) / theta;

  /* standard form */

  if (x<=0.) return 0.;

  return ( exp( -pow( x, -k ) ) );

} /* end of _unur_cdf_extremeII() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_extremeII( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.zeta + pow( DISTR.k/(DISTR.k+1.), 1/DISTR.k ) * DISTR.theta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_extremeII() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_extremeII( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.theta);


  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }

  /* else */
  DISTR.area = ( _unur_cdf_extremeII( DISTR.domain[1],distr) 
		 - _unur_cdf_extremeII( DISTR.domain[0],distr) );
  return 1;
  
} /* end of _unur_upd_area_extremeII() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_extremeII( UNUR_DISTR *distr, double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,0);

  /* check parameter k */
  if (k <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"k <= 0");
    return 0;
  }

  /* check parameter theta */
  if (n_params > 2 && theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    return 0;
  }

  /* copy parameters for standard form */
  DISTR.k = k;

  /* default parameters */
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;

  /* copy optional parameters */
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
  case 2:
    DISTR.zeta = zeta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.zeta;      /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_extremeII() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_extremeII( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_EXTREME_II;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_extremeII_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_extremeII;  /* pointer to PDF                  */
  DISTR.dpdf = _unur_dpdf_extremeII; /* pointer to derivative of PDF    */
  DISTR.cdf  = _unur_cdf_extremeII;  /* pointer to CDF                  */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (!_unur_set_params_extremeII(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = log(DISTR.theta);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.zeta + pow( DISTR.k/(DISTR.k+1.), 1/DISTR.k ) * DISTR.theta;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_extremeII;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_extremeII; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_extremeII; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_extremeII() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
/*---------------------------------------------------------------------------*/
