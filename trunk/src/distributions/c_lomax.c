/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_lomax.c                                                    *
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
 * distr: Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]*
 *                                                                           *
 *  pdf:       f(x) = (x+C)^(-(a+1))                                         *
 *  domain:    x >= 0                                                        *
 *  constant:  a * C^a                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  a > 0   ... shape                                                 *
 *     1:  C > 0   ... location                                              *
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
static const char distr_name[] = "lomax";

/* parameters */
#define a params[0]
#define C params[1]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_lomax(double x, UNUR_DISTR *distr);
static double _unur_dpdf_lomax(double x, UNUR_DISTR *distr);
static double _unur_cdf_lomax(double x, UNUR_DISTR *distr);

static int _unur_upd_mode_lomax( UNUR_DISTR *distr );
static int _unur_upd_area_lomax( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_lomax( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  return ( (x<0.) ? 0. : pow(x+C,-(a+1.)) / NORMCONSTANT );
} /* end of _unur_pdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_lomax( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  return ( (x<0.) ? 0. : -(a+1.) * pow(x+C,-(a+2.)) / NORMCONSTANT );
} /* end of _unur_dpdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_lomax( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  return ( (x<0.) ? 0. : 1. - pow((C/(x+C)),a) );
} /* end of _unur_cdf_lomax() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_lomax( UNUR_DISTR *distr )
{
  DISTR.mode = 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_lomax() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_lomax( UNUR_DISTR *distr )
{
  /* normalization constant */
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }

  /* else */
  DISTR.area = ( _unur_cdf_lomax( DISTR.domain[1],distr) 
		 - _unur_cdf_lomax( DISTR.domain[0],distr) );
  return 1;
  
} /* end of _unur_upd_area_lomax() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_lomax( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOMAX;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = NULL;            /* _unur_stdgen_lomax_init; */

  /* functions */
  DISTR.pdf  = _unur_pdf_lomax;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_lomax; /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_lomax;  /* pointer to CDF               */

  /* default parameters */
  DISTR.C = 1.; 
  
  /* copy parameters */
  DISTR.a = a;
  switch (n_params) {
  case 2:
    DISTR.C = C;
  default:
    n_params = 2;
  }

  /* check parameters */
  if (DISTR.a <= 0. || DISTR.C <= 0. ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"a <= 0 or C <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant */
  NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);

  /* domain */
  DISTR.domain[0] = 0;               /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_lomax; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_lomax; /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_lomax() */

/*---------------------------------------------------------------------------*/
#undef a
#undef C
#undef DISTR
/*---------------------------------------------------------------------------*/
