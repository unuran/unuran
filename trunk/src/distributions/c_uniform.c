/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_uniform.c                                                  *
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
 *  distr: Uniform distribution [3; ch.26, p.276]                            *
 *                                                                           *
 *  pdf:       f(x) = 1 / (b-a)                                              *
 *  domain:    a <= x <= b                                                   *
 *  constant:  1                                                             *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
static const char distr_name[] = "uniform";

/* parameters */
#define a  params[0]
#define b  params[1]

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_uniform(double x, UNUR_DISTR *distr);
static double _unur_dpdf_uniform(double x, UNUR_DISTR *distr);
static double _unur_cdf_uniform(double x, UNUR_DISTR *distr);

static int _unur_upd_mode_uniform( UNUR_DISTR *distr );
static int _unur_upd_area_uniform( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_uniform( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:  /* non standard */
    return ((x < a || x > b) ? 0. : 1./(b-a));
  case 0: default: /* standard */
    return ((x < 0. || x > 1.) ? 0. : 1.);
  }
} /* end of _unur_pdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_uniform( double x, UNUR_DISTR *distr )
{ 
  return 0.;
} /* end of _unur_dpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_uniform( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);
  case 0: default: /* standard */
    if (x<=0.) return 0.;
    if (x>=1.) return 1.;
    return x;
  }
} /* end of _unur_cdf_uniform() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_uniform( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_uniform() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_uniform( UNUR_DISTR *distr )
{
  /* nothing to do */
  return 1;
} /* end of _unur_upd_area_uniform() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_uniform( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0) n_params = 0;
  if (n_params == 1)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"");
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_UNIFORM;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_uniform_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_uniform;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_uniform; /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_uniform;  /* pointer to CDF               */

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
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a >= b");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* domain */
  DISTR.domain[0] = DISTR.a;      /* left boundary  */
  DISTR.domain[1] = DISTR.b;      /* right boundary */

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.a + DISTR.b) / 2.;
  DISTR.area = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_uniform; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_uniform; /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_uniform() */

/*---------------------------------------------------------------------------*/
#undef a
#undef b
#undef DISTR
/*---------------------------------------------------------------------------*/
