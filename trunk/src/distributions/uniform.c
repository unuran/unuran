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
static const char distr_name[] = "uniform";

/* parameters */
#define a  params[0]
#define b  params[1]

/* function prototypes                                                       */
static double _unur_pdf_uniform(double x, double *params, int n_params);
static double _unur_dpdf_uniform(double x, double *params, int n_params);
static double _unur_cdf_uniform(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_uniform( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */
    return ((x < 0. || x > 1.) ? 0. : 1.);
    
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of _unur_pdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_uniform( double x, double *params, int n_params )
{ 
  return 0.;
} /* end of _unur_dpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_uniform( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */
    if (x<=0.) return 0.;
    if (x>=1.) return 1.;
    return x;
    
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of _unur_cdf_uniform() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_uniform( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params != 0 && n_params != 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params>0)
    CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_UNIFORM;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;    /* _unur_stdgen_uniform_init; ???? */

  /* functions */
  DISTR.pdf  = _unur_pdf_uniform;  /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_uniform; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_uniform;  /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.a = 0.;
  DISTR.b = 1.;

  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.a = a;
    DISTR.b = b;
  default:
  }

  /* check parameters a and b */
  if (DISTR.a >= DISTR.b) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"invalid domain: a >= b!");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.a;      /* left boundary  */
  DISTR.domain[1] = DISTR.b;      /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_uniform() */

/*---------------------------------------------------------------------------*/
#undef a
#undef b
/*---------------------------------------------------------------------------*/

