/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_binomial.c                                                 *
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
 *  distr: Binomial distribution  [1; ch.3, p.105]                           *
 *                                                                           *
 *  pmf:       p(k) = (n \choose k) * p^k * (1-p)^(n-k)                      *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  n >= 1                                                            *
 *     1:  0 < p < 1                                                         *
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

static const char distr_name[] = "binomial";

/*---------------------------------------------------------------------------*/
/* parameters */
#define n  params[0]
#define p  params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr
#define LOGNORMCONSTANT (distr->data.discr.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the PMF of the distribution ? */
#ifdef HAVE_UNUR_SF_LN_FACTORIAL
#  define HAVE_PMF
#else
#  undef  HAVE_PMF
#endif

/** TODO: haben wir eine CDF (kannst du da in cephes nachschauen **/
#ifdef HAVE_UNUR_SF_INCOMPLETE_BETA
# define HAVE_CDF
#else
#  undef HAVE_CDF
#endif
/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_binomial(int k, UNUR_DISTR *distr);
#endif
#ifdef HAVE_CDF
static double _unur_cdf_binomial(int k, UNUR_DISTR *distr); 
#endif

static int _unur_upd_mode_binomial( UNUR_DISTR *distr );
static int _unur_upd_sum_binomial( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_binomial(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if ( k<0 || k>n ) return 0.;

  else
    return exp( k * log(p) + (n-k) * log(1.-p) +
                 LOGNORMCONSTANT - _unur_sf_ln_factorial(k) - _unur_sf_ln_factorial(n-k) ) ;

} /* end of _unur_pmf_binomial() */

#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_binomial(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (k<0) return 0.;

  else if (k==0) return exp(n*(log(1.-p)));
  else if(k>=n) return(1.);
  else return(_unur_sf_incomplete_beta(1.-p, n-k, k+1.));


} /* end of _unur_cdf_binomial() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_binomial( UNUR_DISTR *distr )
{
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_binomial() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_binomial( UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  /* log of normalization constant: none */
  LOGNORMCONSTANT = _unur_sf_ln_factorial(n);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.sum = ( _unur_cdf_binomial( DISTR.domain[1],distr) 
		 - _unur_cdf_binomial( DISTR.domain[0]-1,distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_sum_binomial() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_binomial( double *params, int n_params )
{
  register struct unur_distr *distr;
  int nh;

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
  distr->id = UNUR_DISTR_BINOMIAL;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_binomial_init;
   
  /* functions */
#ifdef HAVE_PMF
  DISTR.pmf  = _unur_pmf_binomial;   /* pointer to PMF */
#endif
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_binomial;   /* pointer to CDF */
#endif

  /* copy parameters */
  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closets integer value!!,");
  DISTR.n = nh;  /** TODO: hier eventuel runden (falls nur int erlaubt sind) **/
  DISTR.p = p;

  /* check parameters */
  if (DISTR.p <= 0. || DISTR.p >= 1. || DISTR.n <= 0.) { 
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"p <= 0 || p >= 1 || n <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* domain: [0, infinity] */
  DISTR.domain[0] = 0;           /* left boundary  */
  DISTR.domain[1] = n;          /* right boundary */

  /* log of normalization constant: none */
  LOGNORMCONSTANT = _unur_sf_ln_factorial(n);

  /* mode and sum over PMF */
  DISTR.mode = (int) ((DISTR.n + 1) * DISTR.p);
  DISTR.sum = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_binomial; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_binomial;  /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_binomial() */

/*---------------------------------------------------------------------------*/
#undef p
#undef r
#undef DISTR
/*---------------------------------------------------------------------------*/











