/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cauchy.c                                                     *
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
 *  pdf:       f(x) = 1./( 1 + ((x-theta)/lambda)^2 )                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  pi * lambda                                                   *
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

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "cauchy";

/* parameters */
#define theta  params[0]
#define lambda params[1]

/* function prototypes                                                       */
static double _unur_pdf_cauchy(double x, double *params, int n_params);
static double _unur_dpdf_cauchy(double x, double *params, int n_params);
static double _unur_cdf_cauchy(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_cauchy(double x, double *params, int n_params)
{ 
  /* standardize */
  x = (x - theta) / lambda;

  return (1./((1+x*x)*NORMCONSTANT));

} /* end of _unur_pdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_cauchy(double x, double *params, int n_params)
{
  /* standardize */
  x = (x - theta) / lambda;

  return ( -2.*x/(lambda*(1+x*x)*(1+x*x)*NORMCONSTANT) );

} /* end of _unur_dpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_cauchy(double x, double *params, int n_params)
{
  return ( 0.5 + atan( (x-theta)/lambda )/M_PI );
} /* end of _unur_cdf_cauchy() */

/*---------------------------------------------------------------------------*/

double 
_unur_normconstant_cauchy(double *params, int n_params)
{
  return (M_PI*lambda);
} /* end of _unur_normconstant_cauchy() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cauchy( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params != 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CAUCHY;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = NULL;           /* _unur_stdgen_cauchy_init */;

  /* functions */
  DISTR.pdf  = _unur_pdf_cauchy;   /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_cauchy;  /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_cauchy;   /* pointer to c.d.f.            */

  /* copy parameters */
  DISTR.theta  = theta;
  DISTR.lambda = lambda;

  /* check parameter lambda */
  if (DISTR.lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"lambda <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant */
  DISTR.NORMCONSTANT = M_PI * DISTR.lambda;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta; 
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;   /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_cauchy() */

/*---------------------------------------------------------------------------*/
#undef theta 
#undef lambda
/*---------------------------------------------------------------------------*/

