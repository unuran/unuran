/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      laplace.c                                                    *
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
 *  Laplace distribution [3; ch.24, p.164]                                   *
 *                                                                           *
 *  pdf:       f(x) = exp(- abs(x-theta) / phi )                             *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  2 * phi                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta     ... location                                            *
 *     1:  phi > 0   ... scale                                               *
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
static const char distr_name[] = "laplace";

/* parameters */
#define theta  params[0]
#define phi    params[1]

/* function prototypes                                                       */
static double _unur_pdf_laplace(double x, double *params, int n_params);
static double _unur_dpdf_laplace(double x, double *params, int n_params);
static double _unur_cdf_laplace(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_laplace( double x, double *params, int n_params )
{ 
  register double z;
  z = (x>theta) ? (x-theta)/phi : (theta-x)/phi;
  return exp(-z) / (2.*phi); 
} /* end of _unur_pdf_laplace() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_laplace( double x, double *params, int n_params )
{ 
  register double z;
  z = (x>theta) ? (x-theta)/phi : (theta-x)/phi;

  if (z == 0.)   /* derivative is not defined, but ...                      */
    return 0.;   /* a tangent parallel to x-axis is possible.               */

  return ( ((x>theta) ? -exp(-z)/phi : exp(-z)/phi) / (2.*phi) );
} /* end of unur_cpdf_laplace() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_laplace( double x, double *params, int n_params )
{ 
  register double z;
  z = (x-theta)/phi;
  return ( (x>theta) ? 1.-0.5 * exp(-z) : 0.5*exp(z) );
} /* end of _unur_cdf_laplace() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_laplace( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0 || n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params > 0)
    CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LAPLACE;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = NULL;     /* _unur_stdgen_beta_init; */

  /* functions */
  DISTR.pdf  = _unur_pdf_laplace;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_laplace; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_laplace;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.theta = 0.;
  DISTR.phi   = 1.;
  
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.phi = phi;
  case 1:
    DISTR.theta = theta;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.phi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"scale parameter phi <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
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
} /* end of unur_distr_laplace() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef phi  
/*---------------------------------------------------------------------------*/
