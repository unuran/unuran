/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_slash.c                                                    *
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
 *  distr: Slash distribution [2; ch.12, p.63]                               *
 *                                                                           *
 *  pdf:       f(x) =  (1 - exp(-x^2/2)) / x^2    for x != 0                 *
 *  pdf:       f(x) =  1 / 2                      for x == 0                 *
 *  domain:    -infinity < x < infinity                                      *
 *  constant:  sqrt(2 * pi)                                                  *
 *                                                                           *
 *  parameters: none                                                         *
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

static const char distr_name[] = "slash";

/* parameters */
/* none */

#define DISTR distr->data.cont
#define NORMCONSTANT (distr->data.cont.norm_constant)

/* function prototypes                                                       */
static double _unur_pdf_slash(double x, UNUR_DISTR *distr);
static double _unur_dpdf_slash(double x, UNUR_DISTR *distr);
/*  static double _unur_cdf_slash(double x, UNUR_DISTR *distr); */

/*---------------------------------------------------------------------------*/

double
_unur_pdf_slash(double x, UNUR_DISTR *distr)
{
  if (x == 0.)
    return (0.5 * NORMCONSTANT);
  else
    return ((1. - exp(-x*x/2.)) / (x*x) * NORMCONSTANT);
} /* end of _unur_pdf_slash() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_slash(double x, UNUR_DISTR *distr)
{ 
  register double xsq = x * x;

  if (x == 0.)
    return 0.;
  else
    return ((-2. + exp(-xsq/2.) * (2. + xsq)) / (xsq * x));

  /** TODO: NORMCONSTANT ?? **/

} /* end of _unur_dpdf_slash() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_slash( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params > 0)
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_SLASH;

  /* name of distribution */
  distr->name = distr_name;
             
  /* how to get special generators */
  DISTR.init = _unur_stdgen_slash_init;
   
  /* functions */
  DISTR.pdf  = _unur_pdf_slash;   /* pointer to PDF               */
  DISTR.dpdf = _unur_dpdf_slash;  /* pointer to derivative of PDF */
  /* DISTR.cdf  = _unur_cdf_slash;   pointer to CDF               */

  /* copy parameters: none */

  /* check parameters: none */

  /* number of arguments */
  DISTR.n_params = 0;

  /* normalization constant */
  NORMCONSTANT = 1. / (M_SQRT2 * M_SQRTPI);

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;   /* left boundary  */
  DISTR.domain[1] = INFINITY;    /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   | 
  		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_slash() */

/*---------------------------------------------------------------------------*/
#undef nu
#undef DISTR
/*---------------------------------------------------------------------------*/
