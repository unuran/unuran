/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      student.c                                                    *
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
 *  Student distribution or t-distribution [3; ch. 28; p. 362]               *
 *                                                                           *
 *  pdf:     f(x) = ( 1 + (x^2)/nu )^(-(nu+1)/2)                             *
 *  domain:  -infinity < x < infintiy                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0: a >= 1  ... shape                                                  *
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

static char distr_name[] = "Student distribution";

#define nu (param [0])
/*---------------------------------------------------------------------------*/

double
unur_pdf_student( double x, double *param, int n_param )
{
  CHECK_NULL(param, RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0.);
  
#if CHECKARGS
  if (nu <= 0.) {
    _unur_error(distr_name , UNUR_ERR_DISTR,"shape parameter nu <= 0.");
    return 0.;
  }
#endif

  return pow( (1. + x*x/nu), (-nu-1.)*0.5 );

}  /* end of unur_pdf_student() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_student( double x, double *param, int n_param )
{
  CHECK_NULL(param,RETURN_NULL);
  CHECK_N_PARAMS(n_param,1,0);

#if CHECKARGS
  if (nu <= 0.) {
    _unur_error(distr_name , UNUR_ERR_DISTR,"scale parameter nu <= 0.");
    return 0.;
  }
#endif

  return ( (-nu-1.)*x/nu * pow( (1. + x*x/nu), (-nu-3.)*0.5 ) );

} /* end of unur_dpdf_student() */

/*---------------------------------------------------------------------------*/
#undef nu
/*---------------------------------------------------------------------------*/




