/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_normal.c                                                   *
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
 *  Normal (Gaussian) distribution [2; ch.13, p.80]                          *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * ((x-mu)/sigma)^2 )                         *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  sigma * sqrt(2 pi)                                            *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mu          ... location                                          *
 *     1:  sigma > 0   ... scale                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp( - x^2 / 2)                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  sqrt(2 pi)                                                    *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  mu    = 0.                                                        *
 *     1:  sigma = 1.                                                        *
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

static const char distr_name[] = "normal";

/* parameters */
#define mu    params[0]
#define sigma params[1]

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pdf_normal(double x, double *params, int n_params);
static double _unur_dpdf_normal(double x, double *params, int n_params);
static double _unur_cdf_normal(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_normal( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - mu) / sigma;
  case 0: default: /* standard */
    return exp(-x*x/2. - LOGNORMCONSTANT); 
  }
} /* end of _unur_pdf_normal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_normal( double x, double *params, int n_params )
{
  register double factor = 1.;

  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    factor = 1./sigma;
    x = (x - mu) / sigma;
  case 0: default: /* standard */
    return ( -x * exp(-x*x/2. - LOGNORMCONSTANT) * factor );
  }
} /* end of _unur_dpdf_normal() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_normal( double x, double *params, int n_params ) 
{
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - mu) / sigma;
  case 0: default: /* standard */
    return _unur_cdf_normal_ext(x);
  }
} /* end of _unur_cdf_normal() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_normal( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0) n_params = 0;
  if (n_params > 2)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  if (n_params > 0)
    CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_NORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_normal_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_normal;   /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_normal;  /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_normal;   /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.mu    = 0.;
  DISTR.sigma = 1.;

  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.sigma = sigma;
  case 1:
    DISTR.mu = mu;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.sigma <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT = log(M_SQRTPI * M_SQRT2 * DISTR.sigma);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.mu;
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
} /* end of unur_distr_normal() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef sigma
/*---------------------------------------------------------------------------*/
