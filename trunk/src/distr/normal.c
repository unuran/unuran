/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      normal.c                                                     *
 *                                                                           *
 *   Normalization constants for pdf OMITTED!                                *
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
 *  Normal (Gaussian) distribution [2; ch.13, p.80]                          *
 *                                                                           *
 *  pdf:     f(x) = exp( -1/2 * ((x-mu)/sigma)^2 )                           *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mu          ... location                                          *
 *     1:  sigma > 0   ... scale                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp( - x^2 / 2)                                          *
 *  domain:  -infinity < x < infinity                                        *
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

#include <float.h>
#include <stdlib.h>

#include <unur_distr.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "normal";

#define mu    (params[0])
#define sigma (params[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_normal( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - mu) / sigma;
  case 0:  /* standard */
    return exp(-x*x/2.); 
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_normal( double x, double *params, int n_params )
{
  register double factor = 1.;

  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    factor = 1./sigma;
    x = (x - mu) / sigma;
  case 0:  /* standard */
    return ( -x * exp(-x*x/2.) * factor );
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_normal( double x, double *params, int n_params ) 
{
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - mu) / sigma;
  case 0:  /* standard */
    return _unur_cdf_normal_ext(x);
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_normal() */

/*---------------------------------------------------------------------------*/

double
unur_mode_normal( double *params, int n_params )
{
  switch (n_params) {
  case 2:  /* non standard */
    CHECK_NULL(params,RETURN_NULL);
    return mu;
  case 0:  /* standard */
    return 0.;
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }
} /* end of unur_mode_normal() */

/*---------------------------------------------------------------------------*/

double
unur_area_normal( double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    return M_SQRTPI * M_SQRT2 * sigma;
  case 0:  /* standard */
    return M_SQRTPI * M_SQRT2;
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_normal() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

#define SET_DATA(i,routinename) \
  COOKIE_SET(specialgen+(i),CK_SPECIALGEN_CONT); \
  specialgen[(i)].data.cont.sample = (routinename); \
  specialgen[(i)].data.cont.is_inversion = FALSE; \
  specialgen[(i)].routine_name = #routinename

/*---------------------------------------------------------------------------*/

#define NORMAL_N_VARIANT  7

inline static struct unur_specialgen *
_unur_distr_normal_specialgen( void )
{
  static struct unur_specialgen *specialgen = NULL;

  if (specialgen == NULL) {
    /* allocate memory */
    specialgen = _unur_malloc( NORMAL_N_VARIANT * sizeof(struct unur_specialgen));
    if (!specialgen) return NULL;   /* error */

    /* set data. [0] = DEFAULT */
    SET_DATA(0,unur_cstd_sample_normal_bm);    /* Box-Muller method */
    SET_DATA(1,unur_cstd_sample_normal_pol);   /* Polarmethod with rejection */
    SET_DATA(2,unur_cstd_sample_normal_nquo);  /* "Naive" ratio-of-uniforms */
    SET_DATA(3,unur_cstd_sample_normal_quo);   /* Ratio-of-uniforms with squeeze */
    SET_DATA(4,unur_cstd_sample_normal_leva);  /* Ratio-of-uniforms with quadratic bounding curves */
    SET_DATA(5,unur_cstd_sample_normal_kr);    /* Kindermann-Ramage method */
    SET_DATA(6,unur_cstd_sample_normal_acr);   /* Acceptance-complement ratio */
  }

  return specialgen;

} /* end of _unur_distr_normal_specialgen() */ 

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_normal( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0 || n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params > 0)
    CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_NORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* add list of special generators */
  distr->n_specialgen = NORMAL_N_VARIANT;
  distr->specialgen   = _unur_distr_normal_specialgen();
                
  /* functions */
  DISTR.pdf  = unur_pdf_normal;   /* pointer to p.d.f.            */
  DISTR.dpdf = unur_dpdf_normal;  /* pointer to derivative of p.d.f. */
  DISTR.cdf  = unur_cdf_normal;   /* pointer to c.d.f.            */

  /* default parameters */
  DISTR.params[0] = 0.;        /* default for mu */
  DISTR.params[1] = 1.;        /* default for sigma */

  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.params[1] = sigma;
  case 1:
    DISTR.params[0] = mu;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.params[1] <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = unur_mode_normal(DISTR.params,DISTR.n_params);
  DISTR.area = unur_area_normal(DISTR.params,DISTR.n_params);

  /* domain */
  DISTR.domain[0] = -INFINITY;   /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_normal() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef sigma
/*---------------------------------------------------------------------------*/
