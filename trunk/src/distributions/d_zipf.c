/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_zipf.c                                                     *
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
 *  distr: Zipf distribution (or Zeta distribution)  [1; ch.11.20, p.465]    *
 *                                                                           *
 *  pmf:       p(k) = (k + tau)^(-(rho+1))                                   *
 *  domain:    1 <= k < infinity                                             *
 *  constant:  1 / sum_k=1^infinity  (k+tau)^(rho+1)                       *
 *             [ zeta(a) = sum_1^infinity (x^a) ... Riemann zeta function ]  *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  rho  > 0   ... shape                                              *
 *     1:  tau >= 0                                                          *
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

static const char distr_name[] = "zipf";

/* parameters */
#define rho  params[0]
#define tau  params[1]

#define DISTR distr->data.discr

/* function prototypes                                                       */
static double _unur_pmf_zipf(int k, UNUR_DISTR *distr);
/*  static double _unur_cdf_zipf(int k, UNUR_DISTR *distr);       */

/*---------------------------------------------------------------------------*/

double
_unur_pmf_zipf(int k, UNUR_DISTR *distr)
{ 
  register double *params = DISTR.params;
  return ((k<1) ? 0. : exp( log(k + tau) * (-rho - 1.) ) );
} /* end of _unur_pmf_zipf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_zipf( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1) {
    _unur_error(distr_name,UNUR_ERR_DISTR_NPARAMS,"too few"); return NULL; }
  if (n_params > 2)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
  CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_discr_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_ZIPF;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_zipf_init;
   
  /* functions */
  DISTR.pmf  = _unur_pmf_zipf;   /* pointer to p.m.f.            */
  /* DISTR.cdf  = _unur_cdf_zipf;   pointer to c.d.f.            */

  /* default parameters */
  DISTR.tau = 0.;
  
  /* copy parameters */
  DISTR.rho = rho;
  switch (n_params) {
  case 2:
    DISTR.tau = tau;
  default:
  }

  /* check parameters */
  if (DISTR.rho <= 0. || DISTR.tau < 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"rho <= 0 || tau < 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */

  /* mode and area below p.d.f. */
  /*    DISTR.mode = 0.; */
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = 1.;          /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN );
                 /* UNUR_DISTR_SET_MODE   |  */
                 /* UNUR_DISTR_SET_PDFAREA ); */
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_zipf() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
