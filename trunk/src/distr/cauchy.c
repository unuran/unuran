/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cauchy.c                                                     *
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
 *  Cauchy distribution [2; ch.16, p.299]                                    *
 *                                                                           *
 *  pdf:     f(x) = 1./( 1 + ((x-theta)/lambda)^2 )                          *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta       ... location                                          *
 *     1:  lambda > 0  ... scale                                             *
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

static char distr_name[] = "Cauchy distribution";

#define theta  (params[0])
#define lambda (params[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_cauchy(double x, double *params, int n_params)
{ 
  /* standardize */
  x = (x - theta) / lambda;

  return (1./(1+x*x));

} /* end of unur_pdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_cauchy(double x, double *params, int n_params)
{
  /* standardize */
  x = (x - theta) / lambda;

  return ( -2.*x/(lambda*(1+x*x)*(1+x*x)) );

} /* end of unur_dpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_cauchy(double x, double *params, int n_params)
{
  return ( 0.5 + atan( (x-theta)/lambda )/M_PI );
} /* end of unur_cdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_mode_cauchy(double *params, int n_params)
{
  return theta;
} /* end of unur_mode_cauchy() */

/*---------------------------------------------------------------------------*/

double
unur_area_cauchy(double *params, int n_params)
{
  return (M_PI*lambda);
} /* end of unur_area_cauchy() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cauchy( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params != 2) {
    _unur_warning(NULL,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );

  /* set magiv cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* functions */
  distr->data.cont.pdf  = unur_pdf_cauchy;   /* pointer to p.d.f.            */
  distr->data.cont.dpdf = unur_dpdf_cauchy;  /* pointer to derivative of p.d.f. */
  distr->data.cont.cdf  = unur_cdf_cauchy;   /* pointer to c.d.f.            */

  /* copy parameters */
  distr->data.cont.params[0] = params[0];    /* theta */
  distr->data.cont.params[1] = params[1];    /* lambda */

  /* check parameter lambda */
  if (distr->data.cont.params[1] <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"lambda <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  distr->data.cont.n_params = n_params;

  /* mode and area below p.d.f. */
  distr->data.cont.mode = unur_mode_cauchy(distr->data.cont.params,distr->data.cont.n_params);
  distr->data.cont.area = unur_area_cauchy(distr->data.cont.params,distr->data.cont.n_params);

  /* domain */
  distr->data.cont.domain[0] = -INFINITY;   /* left boundary  */
  distr->data.cont.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cauchy() */

/*---------------------------------------------------------------------------*/
#undef theta 
#undef lambda
/*---------------------------------------------------------------------------*/


