/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_hypergeometric.c                                           *
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
 *  distr: Hypergeometric distribution  [1; ch.6, p.237]                     *
 *                                                                           *
 *  pmf:       p(k) = (M \choose k) * (N-M \choose n-k) / (N \choose n)      *
 *  domain:    max(0,n-N+M) <= k <= min(n,M)                                 *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  N   >= M                                                          *
 *     1:  M   >= 1                                                          *
 *     2:  n   n >= 1                                                        *
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

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "hypergeometric";

/*---------------------------------------------------------------------------*/
/* parameters */
#define N  params[0]
#define M  params[1]
#define n  params[2]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the PMF of the distribution ? */
#ifdef HAVE_UNUR_SF_LN_FACTORIAL
#  define HAVE_PMF
#  define HAVE_SUM
#else
#  undef  HAVE_PMF
#  undef  HAVE_SUM
#endif

/** In Cephes there is no CDF for the hypergeometric distribution**/
#undef HAVE_CDF

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_hypergeometric( int k, const UNUR_DISTR *distr );
#endif
#ifdef HAVE_CDF
static double _unur_cdf_hypergeometric( int k, const UNUR_DISTR *distr ); 
#endif

static int _unur_upd_mode_hypergeometric( UNUR_DISTR *distr );
#ifdef HAVE_SUM
static int _unur_upd_sum_hypergeometric( UNUR_DISTR *distr );
#endif
static int _unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params );

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if ( k<max(0,(n-N+M-0.5)) || k > min(n,M)+0.5 ) 
    return 0.;

  else
    return exp( LOGNORMCONSTANT - _unur_sf_ln_factorial(k) - _unur_sf_ln_factorial(M-k) -
                _unur_sf_ln_factorial(n-k) - _unur_sf_ln_factorial(N-M-n+k) );

} /* end of _unur_pmf_hypergeometric() */

#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_hypergeometric(int k, const UNUR_DISTR *distr)
{ 

  /* Not included in CEPHES-library !!*/

} /* end of _unur_cdf_hypergeometric() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_hypergeometric( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ( (DISTR.n + 1) * (DISTR.M + 1.) / (DISTR.N + 2.) );

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_hypergeometric() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_SUM

int
_unur_upd_sum_hypergeometric( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;

  /* log of normalization constant: none */
  LOGNORMCONSTANT = _unur_sf_ln_factorial(M) + _unur_sf_ln_factorial(N-M) + _unur_sf_ln_factorial(n) +
    _unur_sf_ln_factorial(N-n) - _unur_sf_ln_factorial(N);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.sum = ( _unur_cdf_hypergeometric( DISTR.domain[1],distr) 
		 - _unur_cdf_hypergeometric( DISTR.domain[0]-1,distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_sum_hypergeometric() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_set_params_hypergeometric( UNUR_DISTR *distr, const double *params, int n_params )
{
  int nh;

  /* check number of parameters for distribution */
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,0);

  /* check parameters */
  if (M <= 0. || N <=0. || n <= 0. || n >= N || M >= N ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"M, N, n must be > 0 and n<N M<N");
    return 0;
  }

  /* copy parameters for standard form */

  nh = (int)(N+0.5);
  if(fabs(nh-N)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.N = nh; 

  nh = (int)(M+0.5);
  if(fabs(nh-M)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.M = nh; 

  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closest integer value");
  DISTR.n = nh; 

  /* default parameters: none */
  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = (int) (max(0,(DISTR.n - DISTR.N + DISTR.M + 0.5)));  /* left boundary  */
    DISTR.domain[1] = (int) (min(DISTR.n, DISTR.M) + 0.5);                 /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_hypergeometric() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_hypergeometric( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_HYPERGEOMETRIC;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_hypergeometric_init;
   
  /* functions */
#ifdef HAVE_PMF
  DISTR.pmf  = _unur_pmf_hypergeometric;   /* pointer to PMF */
#endif
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_hypergeometric;   /* pointer to CDF */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_SUM
		 UNUR_DISTR_SET_PMFSUM |
#endif
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (!_unur_set_params_hypergeometric(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
#ifdef HAVE_SUM
  _unur_upd_sum_hypergeometric( distr );
#else
  LOGNORMCONSTANT = 0.;
#endif
  /* log of normalization constant: none */

  /* mode and sum over PMF */
  _unur_upd_mode_hypergeometric(distr);
  DISTR.sum = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_hypergeometric;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_hypergeometric; /* funct for computing mode */
#ifdef HAVE_SUM
  DISTR.upd_sum  = _unur_upd_sum_hypergeometric;  /* funct for computing area */
#endif

  /* return pointer to object */
  return distr;

} /* end of unur_distr_hypergeometric() */

/*---------------------------------------------------------------------------*/
#undef N
#undef M
#undef n
#undef DISTR
/*---------------------------------------------------------------------------*/









