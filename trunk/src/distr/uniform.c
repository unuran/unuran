/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      uniform.c                                                    *
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
 *  Uniform distribution [3; ch.26, p.276]                                   *
 *                                                                           *
 *  pdf:     f(x) = 1 / (b-a)                                                *
 *  domain:  a <= x <= b                                                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  a      ... location                                               *
 *     1:  b (>a) ... location                                               *
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

static char distr_name[] = "Uniform distribution";

#define a (param[0])
#define b (param[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_uniform( double x, double *param, int n_param )
{ 

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (b <= a) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"b <= a");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */

    return ((x < 0. || x > 1.) ? 0. : 1.);
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_uniform( double x, double *param, int n_param )
{ 

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (b <= a) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"b <= a");
      return 0.;
    }
#endif

  case 0:  /* standard */

    return 0.;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_uniform( double x, double *param, int n_param )
{ 

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (b <= a) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"b <= a");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */

    if (x<=0.) return 0.;
    if (x>=1.) return 1.;

    return x;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_mode_uniform( double *param, int n_param )
{ 

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (b <= a) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"b <= a");
      return 0.;
    }
#endif
    return (a+b)/2.;

  case 0:  /* standard */

    return 0.5;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_area_uniform(double *param, int n_param)
{ 

  switch (n_param) {

  case 2:  /* non standard */
    CHECK_NULL(param,RETURN_NULL);
#if CHECKARGS
    if (b <= a) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"b <= a");
      return 0.;
    }
#endif
    return b-a;

  case 0:  /* standard */

    return 1.;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_uniform() */

/*---------------------------------------------------------------------------*/
#undef a
#undef b
/*---------------------------------------------------------------------------*/




