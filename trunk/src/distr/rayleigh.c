/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      rayleigh.c                                                   *
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
 *  Rayleigh distribution [2; ch.18, p.456]                                  *
 *                                                                           *
 *  pdf:     f(x) = x * exp( -1/2 * (x/sigma)^2 )                            *
 *  domain:  0 <= x < infinity                                               *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma > 0   ... scale                                             *
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

static char distr_name[] =  "Rayleigh distribution";

#define sigma (param[0])
/*---------------------------------------------------------------------------*/

double
unur_pdf_rayleigh( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);

#if CHECKARGS
  if (sigma <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    return 0.;
  }
#endif

  return ( (x<=0.) ? 0. : x * exp(-x*x/(2.*sigma*sigma) ) ); 

} /* end of unur_pdf_rayleigh() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_rayleigh(double x, double *param, int n_param)
{ 
  double z;

  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);

#if CHECKARGS
  if (sigma <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    return 0.;
  }
#endif

  z = x*x/(sigma*sigma);

  return ( (x<=0.) ? 0. : exp(-z/2) * (1-z) ); 

} /* end of unur_dpdf_rayleigh() */

/*---------------------------------------------------------------------------*/
#undef sigma
/*---------------------------------------------------------------------------*/
