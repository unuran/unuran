/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_chisquare.c                                                *
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
 *  distr: Chisquare distribution [2; ch.18, p.416]                          *
 *                                                                           *
 *  pdf:       f(x) = x^((nu/2)-1) * exp( -x/2 )                             *
 *  domain:    0 <= x < infinity                                             *
 *  constant:  2^(nu/2) * Gamma(nu/2)                                        *
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

static const char distr_name[] = "chisquare";

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
static double _unur_pdf_chisquare( double x, UNUR_DISTR *distr );
static double _unur_dpdf_chisquare( double x, UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_chisquare( double x, UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_chisquare( UNUR_DISTR *distr );
#ifdef HAVE_AREA
static int _unur_upd_area_chisquare( UNUR_DISTR *distr );
#endif
static int _unur_set_params_chisquare( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_chisquare(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  if (nu == 2.)
    return exp(-x/2. - LOGNORMCONSTANT);

  return (pow(x,nu/2. - 1.) * exp(-x/2. - LOGNORMCONSTANT));

} /* end of _unur_pdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_chisquare(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support */
    return 0.;

  if (nu == 2.)
    return ( -exp(-x/2. - LOGNORMCONSTANT) / 2. );

  return ( pow(x,nu/2. - 2.) * exp(-x/2. - LOGNORMCONSTANT) * (nu - 2. - x)/2. );
} /* end of _unur_dpdf_chisquare() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_chisquare(double x, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;

  if (x <= 0.)
    /* out of support of p.d.f. */
    return 0.;

  return _unur_sf_incomplete_gamma(x/2.,nu/2.);
} /* end of _unur_cdf_chisquare() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_chisquare( UNUR_DISTR *distr )
{
  DISTR.mode = (DISTR.nu >= 2.) ? (DISTR.nu - 2.) : 0.;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_chisquare() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_AREA

int
_unur_upd_area_chisquare( UNUR_DISTR *distr )
{
  /* normalization constant */
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2.);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }

#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_chisquare( DISTR.domain[1],distr) 
		 - _unur_cdf_chisquare( DISTR.domain[0],distr) );
  return 1;
#else
  return 0;
#endif
  
} /* end of _unur_upd_area_chisquare() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_set_params_chisquare( UNUR_DISTR *distr, double *params, int n_params )
{

  /* check number of parameters for distribution */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return 0; }
  if (n_params > 1) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 1; }
  CHECK_NULL(params,0);

  /* check parameter nu */
  if (nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"nu <= 0");
    return 0;
  }

  /* copy parameters for standard form */
  DISTR.nu = nu;

  /* copy optional parameters: none */

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = 0.;          /* left boundary  */
    DISTR.domain[1] = INFINITY;    /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_chisquare() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_chisquare( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CHISQUARE;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = NULL;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_chisquare;   /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_chisquare;  /* pointer to derivative of PDF */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_chisquare;   /* pointer to CDF               */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
#ifdef HAVE_AREA
		 UNUR_DISTR_SET_PDFAREA |
#endif
		 UNUR_DISTR_SET_MODE );
                
  /* set parameters for distribution */
  if (!_unur_set_params_chisquare(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
#ifdef HAVE_AREA
  LOGNORMCONSTANT = _unur_sf_ln_gamma(DISTR.nu/2.) + M_LN2 * (DISTR.nu/2.);
#else
  LOGNORMCONSTANT = 0.;
#endif

  /* mode and area below p.d.f. */
  DISTR.mode = (DISTR.nu >= 2.) ? (DISTR.nu - 2.) : 0.;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_chisquare;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_chisquare; /* funct for computing mode */
#ifdef HAVE_AREA
  DISTR.upd_area  = _unur_upd_area_chisquare; /* funct for computing area */
#endif

  /* return pointer to object */
  return distr;

} /* end of unur_distr_chisquare() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
