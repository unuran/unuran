/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      normal.c                                                     *
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
 *  Normal (Gaussian) distribution [2; ch.13, p.80]                          *
 *                                                                           *
 *  pdf:     f(x) = exp( -1/2 * ((x-mu)/sigma)^2 )                           *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mu          ... location                                          *
 *     1:  sigma > 0   ... scale                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp( - x^2 / 2)                                          *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  mu    = 0.                                                        *
 *     1:  sigma = 1.                                                        *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Fri Dec 17 13:41:15 CET 1999                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
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

static char distr_name[] = "Normal distribution";

#define mu    (param[0])
#define sigma (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_normal( double x, double *param, int n_param )
{ 
  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (sigma <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    x = (x - mu) / sigma;

  case 0:  /* standard */
    return exp(-x*x/2.); 

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_normal( double x, double *param, int n_param )
{
  register double factor = 1.;

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (sigma <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    factor = 1./sigma;
    x = (x - mu) / sigma;

  case 0:  /* standard */
    return ( -x * exp(-x*x/2.) * factor );

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_normal( double x, double *param, int n_param ) 
{
  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (sigma <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    x = (x - mu) / sigma;

  case 0:  /* standard */
    return _unur_cdf_normal_ext(x);

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_mode_normal( double *param, int n_param )
{
  switch (n_param) {
  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
    return mu;
  case 0:  /* standard */
    return 0.;
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }
} /* end of unur_mode_normal() */

/*---------------------------------------------------------------------------*/

double
unur_area_normal( double *param, int n_param )
{
  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (sigma <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    return 2.506628274631*sigma;

  case 0:  /* standard */
    return 2.506628274631;

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_normal() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef sigma
/*---------------------------------------------------------------------------*/








