/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      powerexpon.c                                                 *
 *                                                                           *
 *   Normalization constants for pdf OMITTED!                                *
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
 *  Power-exponential (Subbotin) distribution [3; ch.24, p.195]              *
 *                                                                           *
 *  pdf:     exp(-1/2 * abs((x-theta)/phi) ^ (2/delta) )                     *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta     ... location                                            *
 *     1:  delta > 0 ... shape                                               *
 *     2:  phi > 0   ... scale                                               *
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

static char distr_name[] = "Power-exponential distribution";

#define theta (param[0])
#define delta (param[1])
#define phi   (param[2])
/*---------------------------------------------------------------------------*/

double
unur_pdf_powerexponential( double x, double *param, int n_param )
{ 
  register double z;

  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,3,0.);

#if CHECKARGS
  if (delta <= 0. || phi <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"delta <= 0 or phi <= 0.");
    return 0.;
  }
#endif

  z = (x - theta) / phi;

  return exp( - pow( abs(z), 2./delta ) * 0.5 );

} /* end of unur_pdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_powerexponential( double x, double *param, int n_param )
{
  register double z, tmp;

  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,3,0.);

#if CHECKARGS
  if (delta <= 0. || phi <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"delta <= 0 or phi <= 0.");
    return 0.;
  }
#endif

  z = (x - theta) / phi;

  if (z == 0.)    /* derivative is not defined, but ...        */
    return 0.;    /* a tangent parallel to x-axis is possible. */

  tmp = exp( - pow( abs(z), 2./delta ) * 0.5 ) * pow(abs(z),2./delta-1.) / (delta*phi);

  /* sign ! */
  return ( (z<0.) ? tmp : -tmp );

} /* end of unur_dpdf_powerexponential() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef delta
#undef phi  
/*---------------------------------------------------------------------------*/

