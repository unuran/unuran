/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lognormal.c                                                  *
 *                                                                           *
 *   Normalization constants for pdf OMITTED!                                *
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
 *  Lognormal distribution [2; ch.14, p.208]                                 *
 *                                                                           *
 *  pdf:     f(x) = 1/(x-theta) * exp( -(log(x-theta)-zeta)^2/(2 sigma^2) )  *
 *  domain:  x > theta                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  zeta                                                              *
 *     1:  sigma > 0                                                         *
 *     2:  theta                                                             *
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

static char distr_name[] = "Lognormal distribution";

#define zeta  (param[0])
#define sigma (param[1])
#define theta (param[2])
/*---------------------------------------------------------------------------*/

double
unur_pdf_lognormal( double x, double *param, int n_param )
{ 
  double z;

  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,3,0.);

#if CHECKARGS
  if (sigma <= 0. ) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"sigma <= 0.");
    return 0.;
  }
#endif

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;

  return ( 1./(x-theta) * exp( -z*z/(2.*sigma*sigma) ) );

} /* end of unur_pdf_lognormal() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_lognormal( double x, double *param, int n_param )
{ 
  double z, sigmasqu;

  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,3,0.);

#if CHECKARGS
  if (sigma <= 0. ) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"sigma <= 0.");
    return 0.;
  }
#endif

  if (x <= theta)
    return 0.;

  z = log(x-theta)-zeta;
  sigmasqu = sigma * sigma;

  return ( 1/((x-theta)*(x-theta)) * exp( -z*z/(2*sigmasqu) ) * (1.+z/sigmasqu) );
} /* end of unur_dpdf_lognormal() */

/*---------------------------------------------------------------------------*/
#undef zeta 
#undef sigma
#undef theta
/*---------------------------------------------------------------------------*/



