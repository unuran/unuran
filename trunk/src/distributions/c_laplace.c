/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_laplace.c                                                  *
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
 *  constant:  1 / (2 * phi)                                                 *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta     ... location                                            *
 *     1:  phi > 0   ... scale                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp(- abs(x))                                          *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  2                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *     0:  theta = 0  ... location                                           *
 *     1:  phi   = 1  ... scale                                              *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
static const char distr_name[] = "laplace";

/* parameters */
#define theta  params[0]
#define phi    params[1]

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_laplace( double x, UNUR_DISTR *distr );
static double _unur_dpdf_laplace( double x, UNUR_DISTR *distr );
static double _unur_cdf_laplace( double x, UNUR_DISTR *distr );

static int _unur_upd_mode_laplace( UNUR_DISTR *distr );
static int _unur_upd_area_laplace( UNUR_DISTR *distr );
static int _unur_set_params_laplace( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_laplace( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double z;

  z = (x>theta) ? (x-theta)/phi : (theta-x)/phi;
  return exp(-z) / (2.*phi); 
} /* end of _unur_pdf_laplace() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_laplace( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double z;

  z = (x>theta) ? (x-theta)/phi : (theta-x)/phi;

  if (z == 0.)   /* derivative is not defined, but ...                      */
    return 0.;   /* a tangent parallel to x-axis is possible.               */

  return ( ((x>theta) ? -exp(-z)/phi : exp(-z)/phi) / (2.*phi) );
} /* end of unur_dpdf_laplace() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_laplace( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double z;

  z = (x-theta)/phi;
  return ( (x>theta) ? 1.-0.5 * exp(-z) : 0.5*exp(z) );
} /* end of _unur_cdf_laplace() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_laplace( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_laplace() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_laplace( UNUR_DISTR *distr )
{
  /* normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }

  /* else */
  DISTR.area = ( _unur_cdf_laplace( DISTR.domain[1],distr) 
		 - _unur_cdf_laplace( DISTR.domain[0],distr) );
  return 1;
  
} /* end of _unur_upd_area_laplace() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_laplace( UNUR_DISTR *distr, double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,0);

  /* check parameter phi */
  if (n_params == 2 && phi <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"phi <= 0");
    return 0;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.theta = 0.;
  DISTR.phi   = 1.;

  /* copy optional parameters */
  switch (n_params) {
  case 2:
    DISTR.phi = phi;
  case 1:
    DISTR.theta = theta;
  default:
    n_params = 2;           /* number of parameters for non-standard form */
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_laplace() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_laplace( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LAPLACE;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_laplace_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_laplace;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_laplace; /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_laplace;  /* pointer to CDF               */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA |
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (!_unur_set_params_laplace(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* normalization constant: none */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_laplace;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_laplace; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_laplace; /* funct for computing area */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_laplace() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef phi  
#undef DISTR
/*---------------------------------------------------------------------------*/
