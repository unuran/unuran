/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_weibull.c                                                  *
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
 *  distr: Weibull distribution [2; ch.21, p.628]                            *
 *                                                                           *
 *  pdf:       f(x) = ((x-zeta)/alpha)^(c-1) * exp(-((x-zeta)/alpha)^c)      *
 *  domain:    zeta < x < infinity                                           *
 *  constant:  c / alpha                                                     *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  c     > 0         ... shape                                       *
 *     1:  alpha > 0   (1)   ... scale                                       *
 *     2:  zeta        (0)   ... location                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = x^(c-1) * exp(-x^c)                                    *
 *  domain:    0 < x < infinity                                              *
 *  constant:  c                                                             *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  c     >0    ... shape                                             *
 *                                                                           *
 *     1:  alpha = 1                                                         *
 *     2:  zeta  = 0                                                         *
 *                                                                           *
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
static const char distr_name[] = "weibull";

/* parameters */
#define c      params[0]
#define alpha  params[1]
#define zeta   params[2]

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_weibull(double x, UNUR_DISTR *distr);
static double _unur_dpdf_weibull(double x, UNUR_DISTR *distr);
static double _unur_cdf_weibull(double x, UNUR_DISTR *distr);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_weibull( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x - zeta) / alpha;
  case 1: default: /* standard */
    if (x < 0.)
      return 0.;
    if (x==0. && c==1.)
      return NORMCONSTANT;

    /* else */
    return (exp (-pow (x, c) + (c-1.) * log (x)) * NORMCONSTANT);
  }
} /* end of _unur_pdf_weibull() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_weibull( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double factor = 1.;
  register double xc;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
    /* standardize */
    factor = 1. / alpha;
    x = (x - zeta) / alpha;
  case 1: default: /* standard */
    if (x < 0.)
      return 0.;
    if (x==0. && c==1.)
      return 0.; 

    /* else */
    xc = -pow (x, c);
    return (exp (xc + (c-2.) * log (x)) * (-1. - c * (xc-1.)) * NORMCONSTANT);
  }
} /* end of unur_dpdf_weibull() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_weibull( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 3:  /* non standard */
    /* standardize */
    x = (x - zeta) / alpha;
  case 1: default: /* standard */
    if (x <= 0.)
      return 0.;

    /* else */
    return (1. - exp(-pow (x, c)));
  }
} /* end of _unur_cdf_weibull() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_weibull( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_WEIBULL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_weibull_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_weibull;  /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_weibull; /* pointer to derivative of PDF */
  DISTR.cdf  = _unur_cdf_weibull;  /* pointer to CDF               */

  /* default parameters */
  DISTR.alpha = 1.;
  DISTR.zeta  = 0.;
  
  /* copy parameters */
  DISTR.c = c;
  switch (n_params) {
  case 3:
    DISTR.zeta = zeta;
  case 2:
    DISTR.alpha = alpha;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.c <= 0. || DISTR.alpha <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"c <= 0 || alpha <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant */
  NORMCONSTANT = DISTR.c / DISTR.alpha;

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.c<=1.) ? 0. : DISTR.alpha * pow((DISTR.c - 1.)/DISTR.c, 1./DISTR.c) + DISTR.zeta;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.zeta;      /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_weibull() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
/*---------------------------------------------------------------------------*/
