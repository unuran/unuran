/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_logarithmic.c                                              *
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
 *  distr: Logarithmic distribution  [1; ch.7, p.285]                        *
 *                                                                           *
 *  pmf:       p(k) = theta^k / k                                            *
 *  domain:    1 <= k < infinity                                             *
 *  constant:  -1 / log(1 - theta)                                           *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  0 < theta < 1  ... shape                                          *
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

static const char distr_name[] = "logarithmic";

/* parameters */
#define theta  params[0]

#define DISTR distr->data.discr

/* function prototypes                                                       */
static double _unur_pmf_logarithmic(int k, UNUR_DISTR *distr);
/*  static double _unur_cdf_logarithmic(int k, UNUR_DISTR *distr);       */

/*---------------------------------------------------------------------------*/

double
_unur_pmf_logarithmic(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;
  return ((k<1) ? 0. : pow( theta, k ) / k * NORMCONSTANT);
} /* end of _unur_pmf_logarithmic() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_logarithmic( double *params, int n_params )
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
  distr = _unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_LOGARITHMIC;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_logarithmic_init;
   
  /* functions */
  DISTR.pmf  = _unur_pmf_logarithmic;   /* pointer to p.m.f.            */
  /* DISTR.cdf  = _unur_cdf_logarithmic;   pointer to c.d.f.            */

  /* copy parameters */
  DISTR.theta = theta;

  /* check parameter */
  if (DISTR.theta <= 0. || DISTR.theta >= 1.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0 || theta >= 1");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.NORMCONSTANT = -1. / log( 1.-DISTR.theta);

  /* mode and area below p.d.f. */
  /*    DISTR.mode = 0.; */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 1.;          /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
                 /* UNUR_DISTR_SET_MODE   |  */
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_logarithmic() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
