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
 *     mean  ... mu      (0-vector)                                          *
 *     covar ... Sigma   (identity matrix)                                   *
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

#include <source_distributions.h>

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
    static double _unur_pdf_multinormal( double *x, UNUR_DISTR *distr );
    static double _unur_dpdf_multinormal( double *result, double *x, UNUR_DISTR *distr );
    static int _unur_upd_mode_multinormal( UNUR_DISTR *distr );
    static int _unur_upd_area_multinormal( UNUR_DISTR *distr );
**/

/*---------------------------------------------------------------------------*/
#if 0
/*---------------------------------------------------------------------------*/

double
_unur_pdf_multinormal( double *x, UNUR_DISTR *distr )
{ 
  ;
} /* end of _unur_pdf_multinormal() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_multinormal( double *result, double *x, UNUR_DISTR *distr )
{
  ;
} /* end of _unur_dpdf_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_mode_multinormal( UNUR_DISTR *distr )
{
  ;
} /* end of _unur_upd_mode_multinormal() */

/*---------------------------------------------------------------------------*/

int
_unur_upd_area_multinormal( UNUR_DISTR *distr )
{
  ;
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
unur_distr_multinormal( int dim, double *mean, double *covar )
{
  register struct unur_distr *distr;

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
  if (! (unur_distr_cvec_set_mean( distr, mean ) &&
	 unur_distr_cvec_set_covar( distr, covar ) ) ) {
    unur_distr_free( distr );
    return NULL;
  }

  /* copy other parameters of distribution */
  /* none */

  /* number of other parameters */
  /* DISTR.n_params = 0;  ... default */

  /* domain */

  /* log of normalization constant */
  /* LOGNORMCONSTANT = 1; */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.mean;
  /** TODO **/
  /*    DISTR.volume = 1.; */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_STDDOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFVOLUME );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
