/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lognormal.c                                                  *
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
 *  Lognormal distribution [2; ch.14, p.207]                                 *
 *                                                                           *
 *  pdf:       f(x) = 1/(x-theta) * exp( -(log(x-theta)-zeta)^2/(2 sigma^2) )*
 *  domain:    x > theta                                                     *
 *  constant:  sigma * sqrt(2*pi)                                            *
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
static const char distr_name[] = "lognormal";

/* parameters */
#define zeta   params[0]
#define sigma  params[1]
#define theta  params[2]

/* function prototypes                                                       */
static double _unur_pdf_lognormal(double x, double *params, int n_params);
static double _unur_dpdf_lognormal(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_lognormal( double x, double *params, int n_params )
{ 
  register double z;

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;
  return ( 1./(x-theta) * exp( -z*z/(2.*sigma*sigma) ) / NORMCONSTANT );

} /* end of _unur_pdf_lognormal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_lognormal( double x, double *params, int n_params )
{ 
  register double z, sigmasqu;

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;
  sigmasqu = sigma * sigma;

  return ( 1/((x-theta)*(x-theta)) * exp( -z*z/(2*sigmasqu) ) * (1.+z/sigmasqu) / NORMCONSTANT );
} /* end of _unur_dpdf_lognormal() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_lognormal( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 2 || n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOGNORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;         /* _unur_stdgen_lognormal_init; */

  /* functions */
  DISTR.pdf  = _unur_pdf_lognormal;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_lognormal; /* pointer to derivative of p.d.f. */
  /* DISTR.cdf = _unur_cdf_lognormal; pointer to c.d.f.               */

  /* default parameters */
  DISTR.theta = 0.;        /* default for theta */
  
  /* copy parameters */
  DISTR.zeta = zeta;
  DISTR.sigma = sigma;
  switch (n_params) {
  case 3:
    DISTR.theta = theta;
  default:
    n_params = 3;
  }

  /* check parameter sigma */
  if (DISTR.sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant */
  DISTR.NORMCONSTANT = DISTR.sigma * sqrt(2.*M_PI);

  /* mode and area below p.d.f. */
  /* DISTR.mode = unur_mode_lognormal(DISTR.params,DISTR.n_params); */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.theta;     /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
/*  		 UNUR_DISTR_SET_MODE   | */
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_lognormal() */

/*---------------------------------------------------------------------------*/
#undef zeta 
#undef sigma
#undef theta
/*---------------------------------------------------------------------------*/
