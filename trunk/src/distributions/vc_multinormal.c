/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multinormal.c                                              *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [5] S. Kotz, N. Balakrishnan, and N.L. Johnson                          *
 *       Continuous Multivariate Distributions,                              *
 *       Volume 1: Models and Applications                                   *
 *       John Wiley & Sons, Inc., New York, 2000                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  distr: Multinormal distribution [5; ch.45, p.107]                        *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * (x-mu)^t . Sigma^(-1) . (x-mu) )           * 
 *  domain:    Reals^2                                                       *
 *  constant:  (2 pi)^(dim/2) * sqrt(det(Sigma))                             *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  mean  ... mu      (0-vector)                                      *
 *     1:  covar ... Sigma   (identity matrix)                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:       f(x) = exp( -1/2 * x^t . x )                                  *
 *  domain:    Reals^2                                                       *
 *  constant:  (2 pi)^(dim/2)                                                *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     mean  = (0,...,0)  ... 0-vector                                       *
 *     covar = identity matrix                                               *
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
#include <distr/distr.h>
#include <distr/cvec.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "multinormal";

/*---------------------------------------------------------------------------*/
/* parameters */

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*---------------------------------------------------------------------------*/
/* function prototypes                                                       */

/** TODO:
    static double _unur_pdf_multinormal( const double *x, const UNUR_DISTR *distr );
    static double _unur_dpdf_multinormal( double *result, const double *x, const UNUR_DISTR *distr );
    static int _unur_upd_mode_multinormal( UNUR_DISTR *distr );
    static int _unur_upd_area_multinormal( UNUR_DISTR *distr );
**/

/*---------------------------------------------------------------------------*/
#if 0
/*---------------------------------------------------------------------------*/

double
_unur_pdf_multinormal( double *x, const UNUR_DISTR *distr )
{ 
  ;
} /* end of _unur_pdf_multinormal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_multinormal( double *result, const double *x, const UNUR_DISTR *distr )
{
  ;
} /* end of _unur_dpdf_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_multinormal( UNUR_DISTR *distr )
{
  return UNUR_SUCCESS;
} /* end of _unur_upd_mode_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_multinormal( UNUR_DISTR *distr )
{
  return UNUR_SUCCESS;
} /* end of _unur_upd_area_multinormal() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Make distribution object                                               **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal( int dim, const double *mean, const double *covar )
{
  struct unur_distr *distr;
  struct unur_distr *marginal;

  /* get new (empty) distribution object */
  distr = unur_distr_cvec_new(dim);

  /* check new parameter for generator */
  if (distr == NULL) {
    /* error: dim < 2 */
    return NULL;
  }

  /* set distribution id */
  distr->id = UNUR_DISTR_MNORMAL;

  /* name of distribution */
  distr->name = distr_name;

  /* how to get special generators */
  DISTR.init = NULL;

  /* functions */
  /* DISTR.pdf  = _unur_pdf_multinormal;      pointer to p.d.f.            */
  /* DISTR.dpdf = _unur_dpdf_multinormal;     pointer to derivative of p.d.f. */

  /* copy (and check) parameters */
  if ((unur_distr_cvec_set_mean(distr,mean)!=UNUR_SUCCESS) ||
      (unur_distr_cvec_set_covar(distr,covar)!=UNUR_SUCCESS) ) {
    unur_distr_free( distr );
    return NULL;
  }

  /* set marginal distributions */
  marginal = unur_distr_normal(NULL,0);
  unur_distr_cvec_set_marginals(distr,marginal);
  unur_distr_free(marginal);

  /* copy other parameters of distribution */
  /* none */

  /* number of other parameters */
  /* DISTR.n_params = 0;  ... default */

  /* domain */

  /* log of normalization constant */
  /* LOGNORMCONSTANT = 1; */

  /* mode */
  DISTR.mode = _unur_malloc( distr->dim * sizeof(double) );
  memcpy( DISTR.mode, DISTR.mean, distr->dim * sizeof(double) );

  /* volume below p.d.f. */
  /** TODO **/
  /*    DISTR.volume = 1.; */

  /* indicate which parameters are set (additional to mean and covariance) */
  distr->set |= ( UNUR_DISTR_SET_STDDOMAIN |
		  /* UNUR_DISTR_SET_PDFVOLUME | */
		  UNUR_DISTR_SET_MODE );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
