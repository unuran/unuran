/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_negativebinomial.c                                         *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] N.L. Johnson, S. Kotz and A.W. Kemp                                 *
 *       Univariate Discrete Distributions,                                  *
 *       2nd edition                                                         *
 *       John Wiley & Sons, Inc., New York, 1992                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  distr: Negative Binomial distribution  [1; ch.5.1, p.200]                *
 *                                                                           *
 *  pmf:       p(k) = (k+r-1 \choose r-1) * p^k * (1-p)^r                    *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  0 < p < 1                                                         *
 *     1:      r > 0                                                         *
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

static const char distr_name[] = "negativebinomial";

/*---------------------------------------------------------------------------*/
/* parameters */
#define p  params[0]
#define r  params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
/* #define NORMCONSTANT (distr->data.discr.norm_constant) */

/*---------------------------------------------------------------------------*/
/* do we have the pmf of the distribution ? */
#ifdef HAVE_UNUR_SF_LN_GAMMA
#  define HAVE_PMF
#else
#  undef  HAVE_PMF
#endif

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_negativebinomial(int k, UNUR_DISTR *distr);
#endif
/*  static double _unur_cdf_negativebinomial(int k, UNUR_DISTR *distr);  */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_negativebinomial(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (k<0) return 0.;

  else
    return (pow( p, (double)k ) * pow( 1.-p, r ) 
	    * exp( _unur_sf_ln_gamma(k+r) - _unur_sf_ln_gamma(r) - _unur_sf_ln_gamma(k+1.) ) );
} /* end of _unur_pmf_negativebinomial() */

#endif

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_negativebinomial( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 2) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_NEGATIVEBINOMIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_negativebinomial_init;
   
  /* functions */
#ifdef HAVE_PMF
  DISTR.pmf  = _unur_pmf_negativebinomial;   /* pointer to p.m.f.            */
#endif
  /* DISTR.cdf  = _unur_cdf_negativebinomial;   pointer to c.d.f.            */

  /* copy parameters */
  DISTR.p = p;
  DISTR.r = r;

  /* check parameters */
  if (DISTR.p <= 0. || DISTR.p >= 1. || DISTR.r <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || r <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */

  /* mode and area below p.d.f. */
  /*    DISTR.mode = 0.; */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0.;          /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
                 /* UNUR_DISTR_SET_MODE   |  */
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_negativebinomial() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
