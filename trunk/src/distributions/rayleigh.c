/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      rayleigh.c                                                   *
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
 *  pdf:       f(x) = x * exp( -1/2 * (x/sigma)^2 )                          *
 *  domain:    0 <= x < infinity                                             *
 *  constant:  sigma^2                                                       *
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

#include <unur_defs.h>

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
static const char distr_name[] =  "rayleigh";

/* parameters */
#define sigma  params[0]

/* function prototypes                                                       */
static double _unur_pdf_rayleigh(double x, double *params, int n_params);
static double _unur_dpdf_rayleigh(double x, double *params, int n_params);
static double _unur_pdf_rayleigh(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_rayleigh( double x, double *params, int n_params )
{ 
  return ( (x<=0.) ? 0. : x * exp(-x*x/(2.*sigma*sigma) - LOGNORMCONSTANT ) ); 
} /* end of _unur_pdf_rayleigh() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_rayleigh( double x, double *params, int n_params )
{ 
  register double z;
  z = x*x/(sigma*sigma);
  return ( (x<=0.) ? 0. : exp(-z/2 - LOGNORMCONSTANT) * (1-z) ); 
} /* end of _unur_dpdf_rayleigh() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_rayleigh( double x, double *params, int n_params )
{ 
  return ( (x<=0.) ? 0. : 1. - exp(-x*x/(2.*sigma*sigma)) );
} /* end of _unur_cdf_rayleigh() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_rayleigh( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params != 1) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params > 0)
    CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_RAYLEIGH;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;    /* _unur_stdgen_rayleigh_init; */
                
  /* functions */
  DISTR.pdf  = _unur_pdf_rayleigh;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_rayleigh; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_rayleigh;  /* pointer to c.d.f.               */

  /* copy parameters */
  DISTR.sigma = sigma;

  /* check parameter sigma */
  if (DISTR.sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT =   2. * log(DISTR.sigma);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.sigma;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0.;              /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr; 

#undef DISTR
} /* end of unur_distr_rayleigh() */

/*---------------------------------------------------------------------------*/
#undef sigma
/*---------------------------------------------------------------------------*/

