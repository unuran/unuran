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
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the PMF of the distribution ? */
#ifdef HAVE_UNUR_SF_LN_GAMMA
#  define HAVE_PMF
#  define HAVE_SUM
#else
#  undef  HAVE_PMF
#  undef  HAVE_SUM
#endif

/* no CDF */
#undef HAVE_CDF

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_negativebinomial(int k, UNUR_DISTR *distr);
#endif
#ifdef HAVE_CDF
static double _unur_cdf_negativebinomial(int k, UNUR_DISTR *distr); 
#endif

static int _unur_upd_mode_negativebinomial( UNUR_DISTR *distr );
#ifdef HAVE_SUM
static int _unur_upd_sum_negativebinomial( UNUR_DISTR *distr );
#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_negativebinomial(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (k<0) return 0.;

  else
    return (pow( p, (double)k ) * pow( 1.-p, r ) 
	    * exp( _unur_sf_ln_gamma(k+r) - _unur_sf_ln_gamma(k+1.) - LOGNORMCONSTANT ) );
} /* end of _unur_pmf_negativebinomial() */

#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_negativebinomial(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (k<0) return 0.;

  else
    return 1.;  /** TODO **/

} /* end of _unur_cdf_negativebinomial() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_negativebinomial( UNUR_DISTR *distr )
{
  double m;

  m = (DISTR.r * (1. - DISTR.p) - 1.) / DISTR.p;

  DISTR.mode = (m<0) ? 0 : (int) (m+1);

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_negativebinomial() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_SUM

int
_unur_upd_sum_negativebinomial( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.r);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.sum = ( _unur_cdf_negativebinomial( DISTR.domain[1],distr) 
		 - _unur_cdf_negativebinomial( DISTR.domain[0]-1,distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_sum_negativebinomial() */

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
  DISTR.pmf  = _unur_pmf_negativebinomial;   /* pointer to PMF */
#endif
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_negativebinomial;   /* pointer to CDF */
#endif

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

  /* domain: [0, infinity] */
  DISTR.domain[0] = 0;           /* left boundary  */
  DISTR.domain[1] = INT_MAX;     /* right boundary */

  /* log of normalization constant */
#ifdef HAVE_SUM
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.r);
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and sum over PMF */
  _unur_upd_mode_negativebinomial(distr);
  DISTR.sum = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_negativebinomial; /* funct for computing mode */
#ifdef HAVE_SUM
  DISTR.upd_sum  = _unur_upd_sum_negativebinomial;  /* funct for computing area */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_SUM
		 UNUR_DISTR_SET_PMFSUM |
#endif
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_negativebinomial() */

/*---------------------------------------------------------------------------*/
#undef p
#undef r
#undef DISTR
/*---------------------------------------------------------------------------*/
