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

#include <source_distributions.h>

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
#else
#  undef  HAVE_PMF
#endif

/** In Cephes there is no CDF for the hypergeometric distribution**/
#undef HAVE_CDF

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_hypergeometric(int k, UNUR_DISTR *distr);
#endif
#ifdef HAVE_CDF
static double _unur_cdf_hypergeometric(int k, UNUR_DISTR *distr); 
#endif

static int _unur_upd_mode_hypergeometric( UNUR_DISTR *distr );
static int _unur_upd_sum_hypergeometric( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_hypergeometric(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  double x;

  if ( k<max(0,n-N+M) || k > min(n,M) ) return 0.;


  else
    {

    return exp( LOGNORMCONSTANT - _unur_sf_ln_factorial(k) - _unur_sf_ln_factorial(M-k) -
                _unur_sf_ln_factorial(n-k) - _unur_sf_ln_factorial(N-M-n+k) );

    }

      
      


} /* end of _unur_pmf_hypergeometric() */

#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_hypergeometric(int k, UNUR_DISTR *distr)
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

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_hypergeometric( double *params, int n_params )
{
  register struct unur_distr *distr;
  int nh;

  /* check new parameter for generator */
  if (n_params < 3) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 3; }
  CHECK_NULL(params,NULL);

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

  /* copy parameters */
  nh = (int)(N+0.5);
  if(fabs(nh-N)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"N was rounded to the closets integer value!!,");
  DISTR.N = nh; 
  nh = (int)(M+0.5);
  if(fabs(nh-M)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"M was rounded to the closets integer value!!,");
  DISTR.M = nh; 
  nh = (int)(n+0.5);
  if(fabs(nh-n)>0.001)
    _unur_warning(distr_name,UNUR_ERR_DISTR_DOMAIN,"n was rounded to the closets integer value!!,");
  DISTR.n = nh; 

  /* check parameters */
  if (DISTR.M <= 0. || DISTR.N <=0. || DISTR.n <= 0. || DISTR.n >= DISTR.N || DISTR.M >= DISTR.N ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"M, N, n must be > 0 and n<N M<N");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* domain: [0, infinity] */
  DISTR.domain[0] = max(0,n-N+M);       /* left boundary  */
  DISTR.domain[1] = min(n,M);           /* right boundary */

  /** TODO: wie soll man das mit dem domain machen, wenn die parameter
      der verteilung geaendert werden?? 
      (das ist die erste verteilung wo das problem auftritt.
      einfach auf [0,INT_MAX] setzen?
      oder den benutzer sagen, dass er bei dieser verteilung den domain nicht 
      veraendern darf??
  **/

  /* log of normalization constant: none */
  LOGNORMCONSTANT = _unur_sf_ln_factorial(M) + _unur_sf_ln_factorial(N-M) + _unur_sf_ln_factorial(n) +
    _unur_sf_ln_factorial(N-n) - _unur_sf_ln_factorial(N);

  /* mode and sum over PMF */
  _unur_upd_mode_hypergeometric(distr);
  DISTR.sum = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_hypergeometric; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_hypergeometric;  /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_hypergeometric() */

/*---------------------------------------------------------------------------*/
#undef N
#undef M
#undef n
#undef DISTR
/*---------------------------------------------------------------------------*/









