/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pareto.c                                                     *
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
 *  Pareto distribution (of first kind) [2; ch.20, p.574]                    *
 *                                                                           *
 *  pdf:     f(x) = x^(-(a+1))                                               *
 *  domain:  x >= k                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0   ... location, shape                                       *
 *     1:  a > 0   ... shape                                                 *
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

static char distr_name[] = "Pareto distribution";

#define k (param[0])
#define a (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_pareto( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (k <= 0. || a <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"k <= 0 or a <= 0.");
    return 0.;
  }
#endif

  return ( (x<k) ? 0. : pow(x,-(a+1.)) );

} /* end of unur_pdf_pareto() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_pareto( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,2,0.);

#if CHECKARGS
  if (k <= 0. || a <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"k <= 0 or a <= 0.");
    return 0.;
  }
#endif

  return ( (x<k) ? 0. : (1.-a) * pow(x,-(a+2.)) );
} /* end of unur_dpdf_pareto() */

/*---------------------------------------------------------------------------*/
#undef k
#undef a
/*---------------------------------------------------------------------------*/










