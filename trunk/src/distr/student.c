/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      student.c                                                    *
 *                                                                           *
 *   Normalization constants for pdf OMITTED!                                *
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [3] N.L. Johnson, S. Kotz and N. Balakrishnan                           *
 *       Continuous Univariate Distributions,                                *
 *       Volume 2, 2nd edition                                               *
 *       John Wiley & Sons, Inc., New York, 1995                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Student distribution or t-distribution [3; ch. 28; p. 362]               *
 *                                                                           *
 *  pdf:     f(x) = ( 1 + (x^2)/nu )^(-(nu+1)/2)                             *
 *  domain:  -infinity < x < infintiy                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0: a >= 1  ... shape                                                  *
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

#include <unur_methods.h>
#include <unur_distr.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_umalloc.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
static const char distr_name[] = "student";

#define nu (params [0])
/*---------------------------------------------------------------------------*/

double
unur_pdf_student( double x, double *params, int n_params )
{
  return pow( (1. + x*x/nu), (-nu-1.)*0.5 );
}  /* end of unur_pdf_student() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_student( double x, double *params, int n_params )
{
  return ( (-nu-1.)*x/nu * pow( (1. + x*x/nu), (-nu-3.)*0.5 ) );
} /* end of unur_dpdf_student() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_student( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params != 1) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params > 0)
    CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_STUDENT;

  /* name of distribution */
  distr->name = distr_name;
                
  /* functions */
  DISTR.pdf  = unur_pdf_student;  /* pointer to p.d.f.               */
  DISTR.dpdf = unur_dpdf_student; /* pointer to derivative of p.d.f. */
  /* DISTR.cdf = unur_cdf_student;   pointer to c.d.f.               */

  /* copy parameters */
  DISTR.params[0] = nu;

  /* check parameter sigma */
  if (DISTR.params[0] <= 0.) {
    _unur_error(distr_name , UNUR_ERR_DISTR,"scale parameter nu <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  /* DISTR.mode = unur_mode_exponential(DISTR.params,DISTR.n_params); */
  /* DISTR.area = unur_area_exponential(DISTR.params,DISTR.n_params); */

  /* domain */
  DISTR.domain[0] = -INFINITY;        /* left boundary  */
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
} /* end of unur_distr_student() */

/*---------------------------------------------------------------------------*/
#undef nu
/*---------------------------------------------------------------------------*/

