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

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static char distr_name[] = "gamma";

#define alpha (params[0])
#define beta  (params[1])
#define gamma (params[2])
/*---------------------------------------------------------------------------*/

double
unur_pdf_gamma( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x-gamma) / beta;
  case 1:  /* standard */
    if (x <= 0.)
      return 0.;
    if (alpha == 1.)
      return exp( -x );
    return exp( (alpha-1.)*log(x) - x );
    /*    return ( pow(x,alpha-1.) * exp(-x) ); */

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_gamma( double x, double *params, int n_params )
{
  register double factor = 1.;

  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    factor = 1./beta;
    x = (x-gamma) / beta;
  case 1:  /* standard */
    if (x <= 0.)
      return 0.;
    if (alpha == 1.)
      return( -exp(-x) * factor );
    return ( pow(x,alpha-2.) * exp(-x) *  ((alpha-1.) -x) * factor ); 

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_gamma( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x-gamma) / beta;
  case 1:  /* standard */
    if (x <= 0.)
      return 0.;
    return _unur_cdf_gamma_ext(x,alpha,1.);

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_mode_gamma( double *params, int n_params )
{
  register double mode;

  mode = (alpha >= 1.) ? (alpha - 1.) : 0.;

  switch (n_params) {
  case 3:  /* non standard */
    return (mode * beta) + gamma;
  case 1:  /* standard */
    return mode;

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_gamma() */

/*---------------------------------------------------------------------------*/

double
unur_area_gamma( double *params, int n_params )
{
  switch (n_params) {
  case 3:  /* non standard */
    return exp( _unur_gammaln(alpha) + beta*alpha );
  case 1:  /* standard */
    return exp(_unur_gammaln(alpha));

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_gamma() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gamma( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params < 1 || n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set distribution id */
  distr->id = UNUR_DISTR_GAMMA;

  /* name of distribution */
  distr->name = distr_name;

  /* functions */
  DISTR.pdf  = unur_pdf_gamma;    /* pointer to p.d.f.            */
  DISTR.dpdf = unur_dpdf_gamma;   /* pointer to derivative of p.d.f. */
  DISTR.cdf  = unur_cdf_gamma;    /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.params[1] = 1.;         /* default for beta  */
  DISTR.params[2] = 0.;         /* default for gamma */

  /* copy parameters */
  DISTR.params[0] = alpha;
  switch (n_params) {
  case 3:
    DISTR.params[2] = gamma;
  case 2:
    DISTR.params[1] = beta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameters alpha and beta */
  if (alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"shape parameter alpha <= 0.");
    free( distr ); return NULL;
  }
  if (beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"scale parameter beta <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = unur_mode_gamma(DISTR.params,DISTR.n_params);
  DISTR.area = unur_area_gamma(DISTR.params,DISTR.n_params);

  /* domain */
  DISTR.domain[0] = DISTR.params[2];  /* left boundary  */
  DISTR.domain[1] = INFINITY;         /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_gamma() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
#undef gamma
/*---------------------------------------------------------------------------*/

