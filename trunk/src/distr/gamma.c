/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gamma.c                                                      *
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
 *  Gamma distribution [2; ch.17, p.337]                                     *
 *                                                                           *
 *  pdf:     f(x) = (x-gamma)^(alpha-1) * exp( -(x-gamma)/beta )             *
 *  domain:  x > gamma                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0  ... shape                                              *
 *     1:  beta > 0   ... scale                                              *
 *     2:  gamma      ... location                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = x^(alpha-1) * exp(-x)                                    *
 *  domain:  x > 0                                                           *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  alpha > 0  ... shape                                              *
 *                                                                           *
 *     1:  beta  = 1                                                         *
 *     2:  gamma = 0                                                         *
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

static char distr_name[] = "Gamma distribution";

#define alpha (param[0])
#define beta  (param[1])
#define gamma (param[2])
/*---------------------------------------------------------------------------*/

double
unur_pdf_gamma( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 3:  /* non standard */
#if CHECKARGS
    if (beta <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale beta <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-gamma) / beta;

  case 1:  /* standard */
#if CHECKARGS
    if (alpha <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
      return 0.;
    }
#endif

    if (x <= 0.)
      return 0.;
    
    if (alpha == 1.)
      return exp( -x );
    
    return exp( (alpha-1.)*log(x) - x );
    /*    return ( pow(x,alpha-1.) * exp(-x) ); */

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_gamma( double x, double *param, int n_param )
{
  register double factor = 1.;

  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 3:  /* non standard */
#if CHECKARGS
    if (beta <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale beta <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    factor = 1./beta;
    x = (x-gamma) / beta;

  case 1:  /* standard */
#if CHECKARGS
    if (alpha <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
      return 0.;
    }
#endif

    if (x <= 0.)
      return 0.;

    if (alpha == 1.)
      return( -exp(-x) * factor );

    return ( pow(x,alpha-2.) * exp(-x) *  ((alpha-1.) -x) * factor ); 

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_gamma( double x, double *param, int n_param )
{ 
  CHECK_NULL(param,RETURN_NULL);

  switch (n_param) {

  case 3:  /* non standard */
#if CHECKARGS
    if (beta <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale beta <= 0.");
      return 0.;
    }
#endif
    /* standardize */
    x = (x-gamma) / beta;

  case 1:  /* standard */
#if CHECKARGS
    if (alpha <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
      return 0.;
    }
#endif

    if (x <= 0.)
      return 0.;

    return _unur_cdf_gamma_ext(x,alpha,1.);

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_mode_gamma( double *param, int n_param )
{
  register double mode;

  CHECK_NULL(param,RETURN_NULL);

#if CHECKARGS
  if (alpha <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
    return 0.;
  }
#endif

  mode = (alpha >= 1.) ? (alpha - 1.) : 0.;

  switch (n_param) {
  case 3:  /* non standard */
#if CHECKARGS
    if (beta <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale beta <= 0.");
      return 0.;
    }
#endif
    return (mode * beta) + gamma;
  case 1:  /* standard */
    return mode;

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_area_gamma( double *param, int n_param )
{
  CHECK_NULL(param,RETURN_NULL);

#if CHECKARGS
  if (alpha <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
    return 0.;
  }
#endif

  switch (n_param) {

  case 3:  /* non standard */
#if CHECKARGS
    if (beta <= 0.) {
      _unur_error(distr_name ,UNUR_ERR_DISTR,"scale beta <= 0.");
      return 0.;
    }
#endif
    return exp( _unur_gammaln(alpha) + beta*alpha );

  case 1:  /* standard */
    return exp(_unur_gammaln(alpha));

  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_gamma() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
/*---------------------------------------------------------------------------*/









