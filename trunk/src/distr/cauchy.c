/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cauchy.c                                                     *
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
 *  Cauchy distribution [2; ch.16, p.299]                                    *
 *                                                                           *
 *  pdf:     f(x) = 1./( 1 + ((x-theta)/lambda)^2 )                          *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta       ... location                                          *
 *     1:  lambda > 0  ... scale                                             *
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

static char distr_name[] = "Cauchy distribution";

#define theta  (param[0])
#define lambda (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_cauchy(double x, double *param, int n_param)
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (lambda <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"lambda <= 0.");
    return 0.;
  }
#endif

  /* standardize */
  x = (x - theta) / lambda;

  return (1./(1+x*x));

} /* end of unur_pdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_cauchy(double x, double *param, int n_param)
{
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (lambda <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"lambda <= 0.");
    return 0.;
  }
#endif

  /* standardize */
  x = (x - theta) / lambda;

  return ( -2.*x/(lambda*(1+x*x)*(1+x*x)) );

} /* end of unur_dpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_cauchy(double x, double *param, int n_param)
{
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (lambda <= 0.) {
    _unur_error(distr_name , UNUR_ERR_DISTR, "lambda <= 0.");
    return 0.;
  }
#endif

  return ( 0.5 + atan( (x-theta)/lambda )/M_PI );

} /* end of unur_cdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_mode_cauchy(double *param, int n_param)
{
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (lambda <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"lambda <= 0.");
    return 0.;
  }
#endif

  return theta;

} /* end of unur_mode_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_area_cauchy(double *param, int n_param)
{
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (lambda <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"lambda <= 0.");
    return 0.;
  }
#endif

  return (M_PI*lambda);

} /* end of unur_area_cauchy() */

/*---------------------------------------------------------------------------*/
#undef theta 
#undef lambda
/*---------------------------------------------------------------------------*/


