/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_powerexponential.c                                         *
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
 *  distr: power-exponential distribution                                    *
 *    (also called Subbotin distribution)                                    *
 *                                                                           *
 *  (see also [3; ch.24, p.195] for a different definition;                  *
 *  also called Subbotin distribution)                                       *
 *                                                                           *
 *  pdf:       exp(-abs(x)^tau)                                              *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  2 * Gamma(1+1/tau)                                            *
 *                                                                           *
 *  parameters: 1                                                            *
 *     0:  tau > 0 ... shape                                                 *
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
static const char distr_name[] = "powerexponential";

/* parameters */
#define tau  params[0]

#define DISTR distr->data.cont

/* function prototypes                                                       */
static double _unur_pdf_powerexponential(double x, UNUR_DISTR *distr);
static double _unur_dpdf_powerexponential(double x, UNUR_DISTR *distr);
static double _unur_cdf_powerexponential(double x, UNUR_DISTR *distr);
inline static double _unur_lognormconstant_powerexponential(double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_powerexponential( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  return exp( - pow(fabs(x), tau) - LOGNORMCONSTANT);
} /* end of _unur_pdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_powerexponential( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  register double tmp;

  if (x == 0.)    /* derivative is not defined, but ...        */
    return 0.;    /* a tangent parallel to x-axis is possible. */

  tmp = exp( -pow(fabs(x),tau) - LOGNORMCONSTANT + (tau-1.)*log(fabs(x)) ) * tau;

  /* sign ! */
  return ( (x<0.) ? tmp : -tmp );
} /* end of _unur_dpdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_powerexponential( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;
  register double cdf;

  /* compute cdf(abs(x)) - cdf(0) */
  cdf = _unur_incgamma(pow(fabs(x),tau),1./tau) / 2.;
  return ((x<0.) ? 0.5 - cdf : 0.5 + cdf);

} /* end of _unur_cdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_powerexponential(double *params, int n_params)
{ 
  return  _unur_gammaln(1+1/tau) + M_LN2;
} /* end of _unur_lognormconstant_powerexponential() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_powerexponential( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = _unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_POWEREXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_powerexponential_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_powerexponential;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_powerexponential; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_powerexponential;  /* pointer to c.d.f.               */

  /* copy parameter */
  DISTR.tau = tau;

  /* check parameter sigma */
  if (DISTR.tau <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"tau <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT = _unur_lognormconstant_powerexponential(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  DISTR.mode = 0;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA ); 

  /* return pointer to object */
  return distr;

} /* end of unur_distr_powerexponential() */

/*---------------------------------------------------------------------------*/
#undef tau
#undef DISTR
/*---------------------------------------------------------------------------*/
