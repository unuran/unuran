/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pareto.c                                                     *
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
 *  Pareto distribution (of first kind) [2; ch.20, p.574]                    *
 *                                                                           *
 *  pdf:     f(x) = x^(-(a+1))                                               *
 *  domain:  x >= k                                                          *
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

#include <float.h>
#include <stdlib.h>

#include <unur_distr.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
static char distr_name[] = "pareto";

#define k (params[0])
#define a (params[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_pareto( double x, double *params, int n_params )
{ 
  return ( (x<k) ? 0. : pow(x,-(a+1.)) );
} /* end of unur_pdf_pareto() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_pareto( double x, double *params, int n_params )
{ 
  return ( (x<k) ? 0. : (1.-a) * pow(x,-(a+2.)) );
} /* end of unur_dpdf_pareto() */

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

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set distribution id */
  distr->id = UNUR_DISTR_PARETO;

  /* name of distribution */
  distr->name = distr_name;
                
  /* functions */
  DISTR.pdf  = unur_pdf_pareto;  /* pointer to p.d.f.               */
  DISTR.dpdf = unur_dpdf_pareto; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = NULL;             /* pointer to c.d.f.               */

  /* copy parameters */
  DISTR.params[0] = k;
  DISTR.params[1] = a;

  /* check parameters k and a */
  if (DISTR.params[0] <= 0. || DISTR.params[1] <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"k <= 0 or a <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = 0.;     /* unur_mode_pareto(DISTR.params,DISTR.n_params); */
  DISTR.area = 1.;     /* unur_area_pareto(DISTR.params,DISTR.n_params); */

  /* domain */
  DISTR.domain[0] = DISTR.params[0]; /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN );

/*  		 UNUR_DISTR_SET_MODE   | */
/*  		 UNUR_DISTR_SET_PDFAREA ); */

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_pareto() */

/*---------------------------------------------------------------------------*/
#undef k
#undef a
/*---------------------------------------------------------------------------*/
