/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      exponential.c                                                *
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
 *  Exponential distribution [2; ch.19, p.494]                               *
 *                                                                           *
 *  pdf:     f(x) = exp( - (x-theta)/sigma )                                 *
 *  domain:  x >= theta                                                      *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma > 0  ... scale                                              *
 *     1:  theta      ... location                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp(-x)                                                  *
 *  domain:  x >= 0                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  sigma = 1                                                         *
 *     1:  theta = 0                                                         *
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

static char distr_name[] = "Exponential distribution";

#define sigma (param[0])
#define theta (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_exponential( double x, double *param, int n_param )
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
    x = (x-theta) / sigma;

  case 0:  /* standard */
    return ( (x<0.) ? 0. : exp(-x) );

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_exponential() */

/*---------------------------------------------------------------------------*/
  
double
unur_dpdf_exponential( double x, double *param, int n_param )
{
  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (sigma <= 0.) {
      _unur_error("Exponential distribution",UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
      return 0.;
    }
#endif
    return ( (x<theta) ? 0. : -exp( -(x-theta)/sigma ) / sigma);

  case 0:  /* standard */
    return ( (x<0.) ? 0. : -exp(-x) );

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_exponential( double x, double *param, int n_param )
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
    x = (x-theta) / sigma;

  case 0:  /* standard */
    return ( (x<0.) ? 0. : 1.-exp(-x) );
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_area_exponential( double *param, int n_param )
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
    return sigma;

  case 0:  /* standard */
    return 1.;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_mode_exponential( double *param, int n_param )
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
    return theta;

  case 0:  /* standard */
    return 0.;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_exponential() */

/*---------------------------------------------------------------------------*/
#undef sigma 
#undef theta 
/*---------------------------------------------------------------------------*/






