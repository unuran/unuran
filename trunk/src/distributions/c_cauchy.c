/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_cauchy.c                                                   *
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
 *  distr: Cauchy distribution [2; ch.16, p.299]                             *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "cauchy";

/* parameters */
#define theta  params[0]
#define lambda params[1]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_cauchy(double x, UNUR_DISTR *distr);
static double _unur_dpdf_cauchy(double x, UNUR_DISTR *distr);
static double _unur_cdf_cauchy(double x, UNUR_DISTR *distr);

static int _unur_upd_mode_cauchy( UNUR_DISTR *distr );
static int _unur_upd_area_cauchy( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_cauchy(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:                      /* non standard */
    x = (x - theta) / lambda;  /* -> standardize */
  case 0: default:             /* standard */
    return (1./((1+x*x)*NORMCONSTANT));
  }
} /* end of _unur_pdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_cauchy(double x, UNUR_DISTR *distr)
{
  register double *params = DISTR.params;

  /* standardize */
  x = (x - theta) / lambda;

  return ( -2.*x/(lambda*(1+x*x)*(1+x*x)*NORMCONSTANT) );

} /* end of _unur_dpdf_cauchy() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_cauchy(double x, UNUR_DISTR *distr)
{
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:                      /* non standard */
    x = (x - theta) / lambda;  /* -> standardize */
  case 0: default:             /* standard */
    return ( 0.5 + atan(x)/M_PI );
  }
} /* end of _unur_cdf_cauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_cauchy( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.theta; 
  return 1;
} /* end of _unur_upd_mode_cauchy() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_cauchy( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = M_PI * DISTR.lambda;

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }
  
  else {
    DISTR.area = ( _unur_cdf_cauchy( DISTR.domain[1],distr) 
		   - _unur_cdf_cauchy( DISTR.domain[0],distr) );
    if (DISTR.area <= 0.) {
      /* this must not happen */
      _unur_warning(distr_name,UNUR_ERR_DISTR_SET,"upd area <= 0");
      DISTR.area = 1.;   /* 0 might cause a FPE */
      return 0.;
    }
    else
      return 1;
  }
} /* end of _unur_upd_area_cauchy() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cauchy( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CAUCHY;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_cauchy_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_cauchy;   /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_cauchy;  /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_cauchy;   /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.theta  = 0.;
  DISTR.lambda = 1.;

  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.lambda = lambda;
  case 1:
    DISTR.theta  = theta;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter lambda */
  if (DISTR.lambda <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"lambda <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* domain */
  DISTR.domain[0] = -INFINITY;   /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* normalization constant */
  NORMCONSTANT = M_PI * DISTR.lambda;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta; 
  DISTR.area = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_cauchy; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_cauchy; /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cauchy() */

/*---------------------------------------------------------------------------*/
#undef theta 
#undef lambda
#undef DISTR
/*---------------------------------------------------------------------------*/
