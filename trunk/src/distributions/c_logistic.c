/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_logistic.c                                                 *
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
 *  Logistic distribution [3; ch.23, p.115]                                  *
 *                                                                           *
 *  pdf:       f(x) = exp(-(x-alpha)/beta) * (1 + exp(-(x-alpha)/beta))^(-2) *
 *  cdf:       F(x) = (1 + exp(-(x-alpha)/beta))^(-1)                        *
 *  domain:    infinity < x < infinity                                       *
 *  constant:  1 / beta                                                      *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  beta  > 0   ... scale                                             *
 *     1:  alpha       ... location                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp(-x) * (1 + exp(-x))^(-2)                           *
 *  cdf:       F(x) = (1 + exp(-x))^(-1)                                     *
 *  domain:    infinity < x < infinity                                       *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *     0:  beta  = 1                                                         *
 *     1:  alpha = 0                                                         *
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
static const char distr_name[] = "logistic";

/* parameters */
#define beta   params[0]
#define alpha  params[1]

/* function prototypes                                                       */
static double _unur_pdf_logistic(double x, double *params, int n_params);
static double _unur_dpdf_logistic(double x, double *params, int n_params);
static double _unur_cdf_logistic(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_logistic( double x, double *params, int n_params )
{ 
  register double ex;

  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - alpha) / beta;
  case 0: default: /* standard */
    ex = exp(-x);
    return (NORMCONSTANT * ex / ((1. + ex) * (1. + ex)));
  }
} /* end of _unur_pdf_logistic() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_logistic( double x, double *params, int n_params )
{ 
  register double factor = 1.;
  register double ex;

  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    factor = 1. / beta;
    x = (x - alpha) / beta;
  case 0: default: /* standard */
    ex = exp(x);
    return (factor * NORMCONSTANT * ex * (1. - ex) / ((1.+ex)*(1.+ex)*(1.+ex)));
  }
} /* end of unur_dpdf_logistic() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_logistic( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - alpha) / beta;
  case 0: default: /* standard */
    return ( 1. / (1. + exp(-x)) );
  }
} /* end of _unur_cdf_logistic() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_logistic( double *params, int n_params )
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
  distr->id = UNUR_DISTR_LOGISTIC;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_logistic_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_logistic;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_logistic; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_logistic;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.beta  = 1.;
  DISTR.alpha = 0.;
  
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.alpha = alpha;
  case 1:
    DISTR.beta = beta;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.beta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"beta <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant */
  DISTR.NORMCONSTANT = 1. / DISTR.beta;

  /* mode and area below p.d.f. */
/*    DISTR.mode = 0.; */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
/*   		 UNUR_DISTR_SET_MODE   | */
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_logistic() */

/*---------------------------------------------------------------------------*/
#undef alpha
#undef beta 
/*---------------------------------------------------------------------------*/
