/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_chi.c                                                      *
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
 *  distr: Chi distribution [2; ch.18, p.417]                                *
 *  (Do not confuse with Chi^2 distribution!)                                *
 *                                                                           *
 *  pdf:       f(x) = x^(nu-1) * exp( -x^2/2 )                               *
 *  domain:    0 <= x < infinity                                             *
 *  constant:  2^((nu/2)-1) * Gamma(nu/2)                                    *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  nu > 0   ... shape (degrees of freedom)                           *
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

static const char distr_name[] = "chi";

/*---------------------------------------------------------------------------*/
/* parameters */
#define nu  params[0]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#ifdef HAVE_UNUR_SF_INCOMPLETE_GAMMA
#  define HAVE_CDF
#else
#  undef  HAVE_CDF
#endif

/* can we compute the area below the pdf ? */
#ifdef HAVE_UNUR_SF_LN_GAMMA
#  define HAVE_AREA
#else
#  undef  HAVE_AREA
#endif

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_chi(double x, UNUR_DISTR *distr);
static double _unur_dpdf_chi(double x, UNUR_DISTR *distr);
#ifdef HAVE_CDF
static double _unur_cdf_chi(double x, UNUR_DISTR *distr);
#endif

/*---------------------------------------------------------------------------*/

double
_unur_pdf_chi(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return (pow(x,nu - 1.) * exp(-x*x/2. - LOGNORMCONSTANT));

} /* end of _unur_pdf_chi() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_chi(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  return ( pow(x,nu - 2.) * exp(-x*x/2. - LOGNORMCONSTANT) * (nu - 1. - x*x) );
} /* end of _unur_dpdf_chi() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_chi(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support of p.d.f. */
    return 0.;

  return _unur_sf_incomplete_gamma(x*x/2.,nu/2.);
} /* end of _unur_cdf_chi() */

#endif

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_chi( double *params, int n_params )
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
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CHI;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_chi_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_chi;   /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_chi;  /* pointer to derivative of p.d.f. */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_chi;   /* pointer to c.d.f.            */
#endif

  /* copy parameters */
  DISTR.nu = nu;

  /* check parameter nu */
  if (DISTR.nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
#ifdef HAVE_AREA
  LOGNORMCONSTANT = _unur_sf_ln_gamma(nu/2.) - M_LN2 * (nu/2. - 1.);
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and area below p.d.f. */
  DISTR.mode = (nu >= 1.) ? sqrt(nu - 1.) : 0.;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0.;          /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_chi() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
