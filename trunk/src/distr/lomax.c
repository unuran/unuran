/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lomax.c                                                      *
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
 *  Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]      *
 *                                                                           *
 *  pdf:     f(x) = (x+C)^(-(a+1))                                           *
 *  domain:  x >= 0                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  C > 0   ... location                                              *
 *     1:  a > 0   ... shape                                                 *
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

static char distr_name[] = "Lomax distribution";

#define C (param[0])
#define a (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_lomax( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (C <= 0. || a <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"C <= 0 or a <= 0.");
    return 0.;
  }
#endif

  return ( (x<0.) ? 0. : pow(x+C,-(a+1.)) );

} /* end of unur_pdf_lomax() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_lomax( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (C <= 0. || a <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"C <= 0 or a <= 0.");
    return 0.;
  }
#endif

  return ( (x<0.) ? 0. : -(a+1.) * pow(x+C,-(a+2.)) );

} /* end of unur_dpdf_lomax() */

/*---------------------------------------------------------------------------*/
#undef C
#undef a
/*---------------------------------------------------------------------------*/






