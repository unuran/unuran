/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pareto.c                                                     *
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
 *  Pareto distribution (of first kind) [2; ch.20, p.574]                    *
 *                                                                           *
 *  pdf:       f(x) = x^(-(a+1))                                              *
 *  domain:    x >= k                                                         *
 *  constant:  a * k^a                                                       *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  k > 0   ... location, shape                                       *
 *     1:  a > 0   ... shape                                                 *
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
static const char distr_name[] = "pareto";

/* parameters */
#define k  params[0]
#define a  params[1]

/* function prototypes                                                       */
static double _unur_pdf_pareto(double x, double *params, int n_params);
static double _unur_dpdf_pareto(double x, double *params, int n_params);
static double _unur_cdf_pareto(double x, double *params, int n_params);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_pareto( double x, double *params, int n_params )
{ 
  return ( (x<k) ? 0. : pow(x,-(a+1.))/NORMCONSTANT );
} /* end of _unur_pdf_pareto() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_pareto( double x, double *params, int n_params )
{ 
  return ( (x<k) ? 0. : (1.-a) * pow(x,-(a+2.))/NORMCONSTANT );
} /* end of _unur_dpdf_pareto() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_pareto( double x, double *params, int n_params )
{ 
  return ( (x<k) ? 0. : (1. - pow(k/x,a)) );
} /* end of _unur_cdf_pareto() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_pareto( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params != 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_PARETO;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = NULL;            /* _unur_stdgen_pareto_init; */

  /* functions */
  DISTR.pdf  = _unur_pdf_pareto;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_pareto; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_pareto;  /* pointer to c.d.f.               */

  /* copy parameters */
  DISTR.k = k;
  DISTR.a = a;

  /* check parameters k and a */
  if (DISTR.k <= 0. || DISTR.a <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"k <= 0 or a <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* log of normalization constant */
  DISTR.NORMCONSTANT = DISTR.a * pow(DISTR.k,DISTR.a);

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.k;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = DISTR.k;         /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
  		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_pareto() */

/*---------------------------------------------------------------------------*/
#undef k
#undef a
/*---------------------------------------------------------------------------*/
