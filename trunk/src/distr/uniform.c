/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      uniform.c                                                    *
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
 *  Uniform distribution [3; ch.26, p.276]                                   *
 *                                                                           *
 *  pdf:     f(x) = 1 / (b-a)                                                *
 *  domain:  a <= x <= b                                                     *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  a      ... location                                               *
 *     1:  b (>a) ... location                                               *
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

static char distr_name[] = "Uniform distribution";

#define a (params[0])
#define b (params[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_uniform( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */
    return ((x < 0. || x > 1.) ? 0. : 1.);
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_uniform( double x, double *params, int n_params )
{ 
  return 0.;
} /* end of unur_dpdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_uniform( double x, double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x-a) / (b-a);

  case 0:  /* standard */
    if (x<=0.) return 0.;
    if (x>=1.) return 1.;
    return x;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_mode_uniform( double *params, int n_params )
{ 
  switch (n_params) {
  case 2:  /* non standard */
    return (a+b)/2.;

  case 0:  /* standard */
    return 0.5;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_area_uniform(double *params, int n_params)
{ 
  switch (n_params) {
  case 2:  /* non standard */
    return b-a;

  case 0:  /* standard */
    return 1.;
    
  default:
    _unur_error(distr_name ,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_uniform() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_uniform( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params != 0 && n_params != 2) {
    _unur_warning(NULL,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params>0)
    CHECK_NULL(params,RETURN_NULL);

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );

  /* set magiv cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* functions */
  distr->data.cont.pdf  = unur_pdf_uniform;  /* pointer to p.d.f.            */
  distr->data.cont.dpdf = unur_dpdf_uniform; /* pointer to derivative of p.d.f. */
  distr->data.cont.cdf  = unur_cdf_uniform;  /* pointer to c.d.f.            */

  /* copy parameters */
  switch (n_params) {
  case 0:
    distr->data.cont.params[0] = 0.;         /* default for a */
    distr->data.cont.params[1] = 1.;         /* default for b */
    break;
  case 2:
    distr->data.cont.params[0] = params[0];  /* a */
    distr->data.cont.params[1] = params[1];  /* b */
    break;
  }

  /* check parameters a and b */
  if (distr->data.cont.params[0] >= distr->data.cont.params[1]) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"invalid domain: a >= b!");
    free( distr ); return NULL;
  }

  /* number of arguments */
  distr->data.cont.n_params = n_params;

  /* mode and area below p.d.f. */
  distr->data.cont.mode = unur_mode_uniform(distr->data.cont.params,distr->data.cont.n_params);
  distr->data.cont.area = unur_area_uniform(distr->data.cont.params,distr->data.cont.n_params);

  /* domain */
  distr->data.cont.domain[0] = distr->data.cont.params[0]; /* left boundary  */
  distr->data.cont.domain[1] = distr->data.cont.params[1]; /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );
                
  /* return pointer to object */
  return distr;

} /* end of unur_distr_uniform() */

/*---------------------------------------------------------------------------*/
#undef a
#undef b
/*---------------------------------------------------------------------------*/

