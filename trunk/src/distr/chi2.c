/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      chi2.c                                                       *
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
 *  Chisquare distribution [2; ch.18, p.416]                                 *
 *                                                                           *
 *  pdf:     f(x) = x^((nu/2)-1) * exp( -x/2 )                               *
 *  domain:  0 <= x < infinity                                               *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  nu > 0   ... shape (degrees of freedom)                           *
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

static char distr_name[] = "Chi^2 distribution";

#define nu (param[0])

/*---------------------------------------------------------------------------*/

double
unur_pdf_chisquare(double x, double *param, int n_param)
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);

#if CHECKARGS
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"shape parameter nu <= 0.");
    return 0.;
  }
#endif

  if (x <= 0.)
    return 0.;

  if (nu == 2.)
    return exp(-x/2.);

  return (pow(x,nu/2. - 1.) * exp(-x/2.));

} /* end of unur_pdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_chisquare(double x, double *param, int n_param)
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);

#if CHECKARGS
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"shape parameter nu <= 0.");
    return 0.;
  }
#endif

  if (x <= 0.)
    return 0.;

  if (nu == 2.)
    return ( -exp(-x/2.) / 2. );

  return ( pow(x,nu/2. - 2.) * exp(-x/2.) * (nu - 2. - x)/2. );

} /* end of unur_dpdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_chisquare(double x, double *param, int n_param)
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);

#if CHECKARGS
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"shape parameter nu <= 0.");
    return 0.;
  }
#endif

  if (x <= 0.)
    return 0.;

  return _unur_cdf_chisquare_ext(x,nu);

} /* end of unur_cdf_chisquare() */

/*---------------------------------------------------------------------------*/
#undef nu
/*---------------------------------------------------------------------------*/

