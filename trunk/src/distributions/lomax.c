/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      lomax.c                                                      *
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
 *  Lomax distribution (Pareto distr. of second kind) [2; ch.20, p.575]      *
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

#include <source_unuran.h>
#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
static const char distr_name[] = "lomax";

/* parameters */
#define a params[0]
#define C params[1]

/* function prototypes                                                       */
static double _unur_pdf_lomax(double x, double *params, int n_params);
static double _unur_dpdf_lomax(double x, double *params, int n_params);
static double _unur_cdf_lomax(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_lomax( double x, double *params, int n_params )
{ 
  return ( (x<0.) ? 0. : pow(x+C,-(a+1.)) / NORMCONSTANT );
} /* end of _unur_pdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_lomax( double x, double *params, int n_params )
{ 
  return ( (x<0.) ? 0. : -(a+1.) * pow(x+C,-(a+2.)) / NORMCONSTANT );
} /* end of _unur_dpdf_lomax() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_lomax( double x, double *params, int n_params )
{ 
  return ( (x<0.) ? 0. : 1. - pow((C/(x+C)),a) );
} /* end of _unur_cdf_lomax() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_lomax( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 2)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
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
  DISTR.pdf  = _unur_pdf_lomax;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_lomax; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_lomax;  /* pointer to c.d.f.               */

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
  DISTR.NORMCONSTANT = DISTR.a * pow(DISTR.C,DISTR.a);

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0;               /* left boundary  */
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
} /* end of unur_distr_lomax() */

/*---------------------------------------------------------------------------*/
#undef a
#undef C
/*---------------------------------------------------------------------------*/
