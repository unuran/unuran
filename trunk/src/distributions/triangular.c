/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      triangular.c                                                 *
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
 *  Triangular distribution [3; ch.26, p.297]                                *
 *                                                                           *
 *  pdf:              / 2*x / H           for 0 <= x <= H                    *
 *             f(x) = |                                                      *
 *                    \ 2*(1-x) / (1-H)   for H <= x <= 1                    *
 *  domain:    0 <= x <= 1                                                   *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  H    ... shape                                                    *
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
static const char distr_name[] = "triangular";

/* parameters */
#define H   params[0]    /* shape */

/* function prototypes                                                       */
static double _unur_pdf_triangular(double x, double *params, int n_params);
static double _unur_dpdf_triangular(double x, double *params, int n_params);
static double _unur_cdf_triangular(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_triangular( double x, double *params, int n_params )
{ 
  if (x <= 0.)    return 0.;
  if (x <= H)     return (2.*x/H);
  if (x <= 1.)    return (2.*(1.-x)/(1.-H));
  /* otherwise */ return 0.;
} /* end of _unur_pdf_triangular() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_triangular( double x, double *params, int n_params )
{ 
  if (x <= 0.)    return 0.;
  if (x <= H)     return (2./H);
  if (x <= 1.)    return (-2./(1.-H));
  /* otherwise */ return 0.;
} /* end of unur_dpdf_triangular() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_triangular( double x, double *params, int n_params )
{ 
  if (x <= 0.)    return 0.;
  if (x <= H)     return (x*x/H);
  if (x <= 1.)    return ((H + x * (x-2.))/(H-1.));
  /* otherwise */ return 1.;
} /* end of _unur_cdf_triangular() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_triangular( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0) n_params = 0;
  if (n_params > 1)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  if (n_params > 0)
    CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_TRIANGULAR;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_triangular_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_triangular;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_triangular; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_triangular;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.H = 0.5;   /* default is symmetric triangular distribution */
  
  /* copy parameters */
  switch (n_params) {
  case 1:
    DISTR.H = H;
    n_params = 1;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.H < 0. || DISTR.H > 1.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"H < 0 || H > 1");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.H;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0.;        /* left boundary  */
  DISTR.domain[1] = 1.;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_triangular() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef phi  
/*---------------------------------------------------------------------------*/
