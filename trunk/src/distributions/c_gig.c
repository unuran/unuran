/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gig.c                                                      *
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
 *  distr: Generalized Inverse Gaussian (GIG) distribution [2; ch.15, p.284] *
 *                                                                           *
 *  pdf:       f(x) = x^(theta-1) * exp( -1/2 * omega * (x/eta + eta/x))     *
 *  domain:    0 < x < infinity                                              *
 *  constant:  2 * eta^theta K_theta(omega)                                  *
 *             [K_theta(.) ... modified Bessel function of third kind]       *
 * K_theta(x) = 1/2 * int_-inf^inf  cosh(theta*u) * exp(-x*cosh(u)) du       *
 *                                                                           *
 *                              inf                                          *
 *                               -                                           *
 *                          1   |                                            *
 *         K_theta(omega) = -   | cosh(theta*u) * exp(-omega*cosh(u)) du     *
 *                          2   |                                            *
 *                             -                                             *
 *                            -inf                                           *
 *                                                                           *
 *                                                                           *
 *  parameters: 3                                                            *
 *     0:  theta            ... shape                                        *
 *     1:  omega > 0        ... scale                                        *
 *     2:  eta   > 0   (1)  ... shape                                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Standard form:                                                           *
 *                                                                           *
 *  pdf:       f(x) = x^(theta-1) * exp( -1/2 * omega * (x + 1/x))           *
 *  domain:    0 < x < infinity                                              *
 *  constant:  2 * K_theta(omega)                                            *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  theta       ... shape                                             *
 *     1:  omega > 0   ... scale                                             *
 *                                                                           *
 *     2:  eta   = 1                                                         *
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

static const char distr_name[] = "gig";

/* parameters */
#define theta  params[0]    /* shape */
#define omega  params[1]    /* scale */
#define eta    params[2]    /* shape */

#define DISTR distr->data.cont
/* #define NORMCONSTANT (distr->data.cont.norm_constant) */

/* function prototypes                                                       */
static double _unur_pdf_gig(double x, UNUR_DISTR *distr);
static double _unur_dpdf_gig(double x, UNUR_DISTR *distr);
/*  static double _unur_cdf_gig(double x, UNUR_DISTR *distr); */

/*---------------------------------------------------------------------------*/

double
_unur_pdf_gig(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return (exp( (theta-1.) * log(x) - 0.5 * omega * (x/eta + eta/x) ));

} /* end of _unur_pdf_gig() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_gig(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( exp( (theta-3.) * log(x) - 0.5 * omega * (x/eta + eta/x) )
	   * (eta*eta*omega + 2.*eta*(theta-1.)*x - omega*x*x) / (2*eta) );

} /* end of _unur_dpdf_gig() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gig( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_GIG;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_gig_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_gig;   /* pointer to PDF                */
  DISTR.dpdf = _unur_dpdf_gig;  /* pointer to derivative of PDF  */
  DISTR.cdf  = NULL;            /* _unur_cdf_gig; pointer to CDF */

  /* default parameters */
  DISTR.eta  = 1.;

  /* copy parameters */
  DISTR.theta = theta;
  DISTR.omega = omega;
  switch (n_params) {
  case 3:
    DISTR.eta = eta;
    n_params = 3;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameters omega and eta */
  if (DISTR.omega <= 0. || DISTR.eta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"omega <= 0 || eta <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  /*    DISTR.LOGNORMCONSTANT = ? */

  /* mode and area below p.d.f. */
  /*    DISTR.mode = ? */
  /*    DISTR.area = ? */

  /* domain */
  DISTR.domain[0] = 0.;          /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN );
		 /* UNUR_DISTR_SET_MODE   |  */
		 /* UNUR_DISTR_SET_PDFAREA ); */
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_gig() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
