/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_exponential.c                                              *
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
 *  distr: Exponential distribution [2; ch.19, p.494]                        *
 *                                                                           *
 *  pdf:       f(x) = exp( - (x-theta)/sigma )                               *
 *  domain:    x >= theta                                                    *
 *  constant:  sigma                                                         *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma > 0  (1)  ... scale                                         *
 *     1:  theta      (0)  ... location                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp(-x)                                                  *
 *  domain:  x >= 0                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  sigma = 1                                                         *
 *     1:  theta = 0                                                         *
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
static const char distr_name[] = "exponential";

/* parameters */
#define sigma  params[0]
#define theta  params[1]

#define DISTR distr->data.cont

/* function prototypes                                                       */
static double _unur_pdf_exponential(double x, UNUR_DISTR *distr);
static double _unur_dpdf_exponential(double x, UNUR_DISTR *distr);
static double _unur_cdf_exponential(double x, UNUR_DISTR *distr);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_exponential( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:                  /* non standard */
    x = (x-theta) / sigma; /* -> standardize */
  case 0: default:         /* standard */
    return ( (x<0.) ? 0. : exp(-x) / sigma );
  }
} /* end of _unur_pdf_exponential() */

/*---------------------------------------------------------------------------*/
  
double
_unur_dpdf_exponential( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:                  /* non standard */
    return ( (x<theta) ? 0. : -exp( -(x-theta)/sigma ) / (sigma*sigma));
  case 0: default:         /* standard */
    return ( (x<0.) ? 0. : -exp(-x) );
  }
} /* end of _unur_dpdf_exponential() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_exponential( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:                  /* non standard */
    x = (x-theta) / sigma; /* -> standardize */
  case 0: default:         /* standard */
    return ( (x<0.) ? 0. : 1.-exp(-x) );
  }
} /* end of _unur_cdf_exponential() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_exponential( double *params, int n_params )
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
  distr->id = UNUR_DISTR_EXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_exponential_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_exponential;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_exponential; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_exponential;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.sigma = 1.;
  DISTR.theta = 0.;
  
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
  case 1:
    DISTR.sigma = sigma;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.sigma <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.theta;   /* theta */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.theta;     /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_exponential() */

/*---------------------------------------------------------------------------*/
#undef sigma 
#undef theta 
#undef DISTR
/*---------------------------------------------------------------------------*/
