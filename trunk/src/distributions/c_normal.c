/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_normal.c                                                   *
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
 *  distr: Normal (Gaussian) distribution [2; ch.13, p.80]                   *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * ((x-mu)/sigma)^2 )                         *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  1 / (sigma * sqrt(2 pi))                                      *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mu          (0.)  ... location                                    *
 *     1:  sigma > 0   (1.)  ... scale                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp( - x^2 / 2)                                        *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  sqrt(2 pi)                                                    *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  mu    = 0.                                                        *
 *     1:  sigma = 1.                                                        *
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

static const char distr_name[] = "normal";

/*---------------------------------------------------------------------------*/
/* parameters */
#define mu    params[0]
#define sigma params[1]

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont
#define LOGNORMCONSTANT (distr->data.cont.norm_constant)

/*---------------------------------------------------------------------------*/
/* do we have the cdf of the distribution ? */
#ifdef HAVE_UNUR_SF_CDFNORMAL
#  define HAVE_CDF
#else
#  undef  HAVE_CDF
#endif

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */
static double _unur_pdf_normal( double x, UNUR_DISTR *distr );
static double _unur_dpdf_normal( double x, UNUR_DISTR *distr );
#ifdef HAVE_CDF
static double _unur_cdf_normal( double x, UNUR_DISTR *distr );
#endif

static int _unur_upd_mode_normal( UNUR_DISTR *distr );
static int _unur_upd_area_normal( UNUR_DISTR *distr );
static int _unur_set_params_normal( UNUR_DISTR *distr, double *params, int n_params );

/*---------------------------------------------------------------------------*/

double
_unur_pdf_normal( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - mu) / sigma;

  /* standard form */

  return exp(-x*x/2. - LOGNORMCONSTANT); 

} /* end of _unur_pdf_normal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_normal( double x, UNUR_DISTR *distr )
{
  register double *params = DISTR.params;
  register double factor = 1.;

  if (DISTR.n_params > 0) {
    /* standardize */
    factor = 1./sigma;
    x = (x - mu) / sigma;
  }

  /* standard form */

  return ( -x * exp(-x*x/2. - LOGNORMCONSTANT) * factor );

} /* end of _unur_dpdf_normal() */

/*---------------------------------------------------------------------------*/

#ifdef HAVE_CDF

double
_unur_cdf_normal( double x, UNUR_DISTR *distr ) 
{
  register double *params = DISTR.params;

  if (DISTR.n_params > 0)
    /* standardize */
    x = (x - mu) / sigma;

  /* standard form */

  return _unur_sf_cdfnormal(x);

} /* end of _unur_cdf_normal() */

#endif

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_normal( UNUR_DISTR *distr )
{
  DISTR.mode = DISTR.mu;

  /* mode must be in domain */
  if (DISTR.mode < DISTR.domain[0]) 
    DISTR.mode = DISTR.domain[0];
  else if (DISTR.mode > DISTR.domain[1]) 
    DISTR.mode = DISTR.domain[1];

  return 1;
} /* end of _unur_upd_mode_normal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_normal( UNUR_DISTR *distr )
{
  /* log of normalization constant */
  LOGNORMCONSTANT = log(M_SQRTPI * M_SQRT2 * DISTR.sigma);

  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.area = 1.;
    return 1;
  }

#ifdef HAVE_CDF
  /* else */
  DISTR.area = ( _unur_cdf_normal( DISTR.domain[1],distr) 
		 - _unur_cdf_normal( DISTR.domain[0],distr) );
  return 1;
#else
  return 0;
#endif
} /* end of _unur_upd_area_normal() */

/*---------------------------------------------------------------------------*/

int
_unur_set_params_normal( UNUR_DISTR *distr, double *params, int n_params )
{
  /* check number of parameters for distribution */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,0);

  /* check parameter sigma */
  if (n_params > 1 && sigma <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR_DOMAIN,"sigma <= 0");
    return 0;
  }

  /* copy parameters for standard form: none */

  /* default parameters */
  DISTR.mu    = 0.;
  DISTR.sigma = 1.;

  /* copy optional parameters */
  switch (n_params) {
  case 2:
    DISTR.sigma = sigma;
  case 1:
    DISTR.mu = mu;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
    break;
  }

  /* store number of parameters */
  DISTR.n_params = n_params;

  /* set (standard) domain */
  if (distr->set & UNUR_DISTR_SET_STDDOMAIN) {
    DISTR.domain[0] = -INFINITY;       /* left boundary  */
    DISTR.domain[1] = INFINITY;        /* right boundary */
  }

  return 1;
} /* end of _unur_set_params_normal() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_normal( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_NORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = _unur_stdgen_normal_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_normal;   /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_normal;  /* pointer to derivative of PDF */
#ifdef HAVE_CDF
  DISTR.cdf  = _unur_cdf_normal;   /* pointer to CDF               */
#endif

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* set parameters for distribution */
  if (!_unur_set_params_normal(distr,params,n_params)) {
    free(distr);
    return NULL;
  }

  /* log of normalization constant */
  LOGNORMCONSTANT = log(M_SQRTPI * M_SQRT2 * DISTR.sigma);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.mu;
  DISTR.area = 1.;

  /* function for setting parameters and updating domain */
  DISTR.set_params = _unur_set_params_normal;

  /* function for updating derived parameters */
  DISTR.upd_mode  = _unur_upd_mode_normal; /* funct for computing mode */
  DISTR.upd_area  = _unur_upd_area_normal; /* funct for computing area */
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_normal() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef sigma
#undef DISTR
/*---------------------------------------------------------------------------*/
