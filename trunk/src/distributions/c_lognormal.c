/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_lognormal.c                                                *
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
 *  distr: Lognormal distribution [2; ch.14, p.207]                          *
 *                                                                           *
 *  pdf:       f(x) = 1/(x-theta) * exp( -(log(x-theta)-zeta)^2/(2 sigma^2) )*
 *  domain:    x > theta                                                     *
 *  constant:  1 / (sigma * sqrt(2*pi))                                      *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  zeta                                                              *
 *     1:  sigma > 0                                                         *
 *     2:  theta         ... location                                        *
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
static const char distr_name[] = "lognormal";

/* parameters */
#define zeta   params[0]
#define sigma  params[1]
#define theta  params[2]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_lognormal( double x, const UNUR_DISTR *distr );
static double _unur_dpdf_lognormal( double x, const UNUR_DISTR *distr );

static int _unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double z;

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;
  return ( 1./(x-theta) * exp( -z*z/(2.*sigma*sigma) ) / NORMCONSTANT );

} /* end of _unur_pdf_lognormal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_lognormal( double x, const UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double z, sigmasqu;

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;
  sigmasqu = sigma * sigma;

  return ( 1/((x-theta)*(x-theta)) * exp( -z*z/(2*sigmasqu) ) * (1.+z/sigmasqu) / NORMCONSTANT );
} /* end of _unur_dpdf_lognormal() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_logistic( UNUR_DISTR *distr, const double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,0);

  /* check parameter sigma */
  if (sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return 0;
  }

  /* copy parameters for standard form */
  DISTR.zeta = zeta;
  DISTR.sigma = sigma;

  /* default parameters */
  DISTR.theta = 0.;        /* default for theta */

  /* copy optional parameters */
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
  default:
    n_params = 3;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = DISTR.theta;     /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_logistic() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_lognormal( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOGNORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;         /* _unur_stdgen_lognormal_init; */

  /* functions */
  DISTR.pdf  = _unur_pdf_lognormal;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_lognormal; /* pointer to derivative of PDF */
  /* DISTR.cdf = _unur_cdf_lognormal; pointer to CDF                 */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
/*  		 UNUR_DISTR_SET_MODE   | */
  		 UNUR_DISTR_SET_PDFAREA );


  /* set parameters for distribution */
  if (!_unur_set_params_logistic(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* normalization constant */
  NORMCONSTANT = DISTR.sigma * sqrt(2.*M_PI);

  /* mode and area below p.d.f. */
  /* DISTR.mode = unur_mode_lognormal(DISTR.params,DISTR.n_params); */
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_logistic;

  /* function for updating derived parameters */
  /* DISTR.upd_mode  = _unur_upd_mode_logistic;    funct for computing mode */
  /* DISTR.upd_area  = _unur_upd_area_logistic;    funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_lognormal() */

/*---------------------------------------------------------------------------*/
#undef zeta 
#undef sigma
#undef theta
#undef DISTR
/*---------------------------------------------------------------------------*/
