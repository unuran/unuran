/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_poisson.c                                                  *
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
 *  distr: Poisson distribution  [1; ch.4, p.151]                            *
 *                                                                           *
 *  pmf:       p(k) = theta^k / k!                                           *
 *  domain:    0 <= k < infinity                                             *
 *  constant:  exp( -theta )                                                 *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  theta > 0  ... shape                                              *
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

static const char distr_name[] = "poisson";

/* parameters */
#define theta  params[0]

#define DISTR distr->data.discr
/*  #define LOGNORMCONSTANT (distr->data.discr.norm_constant) */

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#if defined(HAVE_UNUR_SF_GAMMA) && defined(HAVE_UNUR_SF_INCOMPLETE_GAMMA)
#  define HAVE_CDF
#else
#  undef  HAVE_CDF
#endif

/* can we compute the area below the pdf ? */
#ifdef HAVE_UNUR_SF_LN_FACTORIAL
#  define HAVE_PMF
#else
#  undef  HAVE_PMF
#endif

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
#ifdef HAVE_PMF
static double _unur_pmf_poisson(int k, UNUR_DISTR *distr);
#endif
#ifdef HAVE_CDF
static double _unur_cdf_poisson(int k, UNUR_DISTR *distr);      
#endif

static int _unur_upd_mode_poisson( UNUR_DISTR *distr );
static int _unur_upd_sum_poisson( UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

#ifdef HAVE_PMF

double
_unur_pmf_poisson(int k, UNUR_DISTR *distr)
{ 
  if (k>=0)
    return exp( -DISTR.theta + k * log(DISTR.theta) - _unur_sf_ln_factorial(k) );
  else
    return 0.;
} /* end of _unur_pmf_poisson() */

#endif

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_poisson(int k, UNUR_DISTR *distr)
{ 
  if (k>=0)
    return ( (k+1.) * _unur_sf_incomplete_gamma(k,DISTR.theta) / _unur_sf_ln_gamma(k+2.) );
  else
    return 0.;
} /* end of _unur_cdf_poisson() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_poisson( UNUR_DISTR *distr )
{
  DISTR.mode = (int) DISTR.theta;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  /* o.k. */
  return 1;
} /* end of _unur_upd_mode_poisson() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_sum_poisson( UNUR_DISTR *distr )
{
  /* log normalization constant: none */

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.sum = 1.;
    return 1;
  }
  
#ifdef HAVE_CDF
  /* else */
  DISTR.sum = ( _unur_cdf_poisson( DISTR.domain[1],distr) 
		 - _unur_cdf_poisson( DISTR.domain[0]-1,distr) );
  return 1;
#else
  return 0;
#endif

} /* end of _unur_upd_sum_poisson() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_poisson( double *params, int n_params )
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
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_POISSON;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_poisson_init;
   
  /* functions */
#ifdef HAVE_PMF
  DISTR.pmf  = _unur_pmf_poisson;   /* pointer to PMF */
#endif
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_poisson;   /* pointer to CDF */
#endif

  /* copy parameters */
  DISTR.theta = theta;

  /* check parameter */
  if (DISTR.theta <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* domain: [0, inifinty] */
  DISTR.domain[0] = 0;           /* left boundary  */
  DISTR.domain[1] = INT_MAX;     /* right boundary */

  /* log normalization constant: none */

  /* mode and sum over PMF */
  DISTR.mode = (int) DISTR.theta;
  DISTR.sum = 1.;

  /* function for updating derived parameters */
  DISTR.upd_mode = _unur_upd_mode_poisson; /* funct for computing mode */
  DISTR.upd_sum  = _unur_upd_sum_poisson;  /* funct for computing area */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PMFSUM |
		 UNUR_DISTR_SET_MODE );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_poisson() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef DISTR
/*---------------------------------------------------------------------------*/
