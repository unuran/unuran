/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      chi2.c                                                       *
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
 *  Chisquare distribution [2; ch.18, p.416]                                 *
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

#include <unur_defs.h>

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "chisquare";

/* parameters */
#define nu  params[0]

/* function prototypes                                                       */
static double _unur_pdf_chisquare(double x, double *params, int n_params);
static double _unur_dpdf_chisquare(double x, double *params, int n_params);
static double _unur_mode_chisquare(double *params, int n_params);
static double _unur_lognormconstant_chisquare(double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_chisquare(double x, double *params, int n_params)
{ 
  if (x <= 0.)
    return 0.;

  if (nu == 2.)
    return exp(-x/2. - LOGNORMCONSTANT);

  return (pow(x,nu/2. - 1.) * exp(-x/2. - LOGNORMCONSTANT));

} /* end of _unur_pdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_chisquare(double x, double *params, int n_params)
{ 
  if (x <= 0.)
    return 0.;

  if (nu == 2.)
    return ( -exp(-x/2. - LOGNORMCONSTANT) / 2. );

  return ( pow(x,nu/2. - 2.) * exp(-x/2. - LOGNORMCONSTANT) * (nu - 2. - x)/2. );
} /* end of _unur_dpdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_chisquare(double x, double *params, int n_params)
{ 
  if (x <= 0.)
    return 0.;

  return _unur_cdf_chisquare_ext(x,nu);
} /* end of _unur_cdf_chisquare() */

/*---------------------------------------------------------------------------*/

double
_unur_mode_chisquare( double *params, int n_params )
{
  return (nu >= 2.) ? (nu/4. - 0.5) : 0.;
} /* end of _unur_mode_chisquare() */

/*---------------------------------------------------------------------------*/

double
_unur_lognormconstant_chisquare( double *params, int n_params )
{
  return ( _unur_gammaln(nu/2.) - M_LN2 * (nu/2.));
} /* end of _unur_lognormconstant_chisquare() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_chisquare( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  CHECK_NULL(params,RETURN_NULL);
  if (n_params != 1) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_CHISQUARE;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = NULL;            /* _unur_stdgen_chisquare_init; */
   
  /* functions */
  DISTR.pdf  = _unur_pdf_chisquare;   /* pointer to p.d.f.            */
  DISTR.dpdf = _unur_dpdf_chisquare;  /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_chisquare;   /* pointer to c.d.f.            */

  /* copy parameters */
  DISTR.nu = nu;

  /* check parameter lambda */
  if (DISTR.nu <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"shape parameter nu <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.LOGNORMCONSTANT = _unur_lognormconstant_chisquare(DISTR.params,DISTR.n_params);

  /* mode and area below p.d.f. */
  DISTR.mode = _unur_mode_chisquare(DISTR.params,DISTR.n_params);
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 0        ;   /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_chi2square() */

/*---------------------------------------------------------------------------*/
#undef nu
/*---------------------------------------------------------------------------*/

