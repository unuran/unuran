/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      exponential.c                                                *
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
 *  Exponential distribution [2; ch.19, p.494]                               *
 *                                                                           *
 *  pdf:     f(x) = exp( - (x-theta)/sigma )                                 *
 *  domain:  x >= theta                                                      *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  sigma > 0  ... scale                                              *
 *     1:  theta      ... location                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  standard form                                                            *
 *                                                                           *
 *  pdf:     f(x) = exp(-x)                                                  *
 *  domain:  x >= 0                                                          *
 *                                                                           *
 *  parameters:                                                              *
 *     none                                                                  *
 *                                                                           *
 *     0:  sigma = 1                                                         *
 *     1:  theta = 0                                                         *
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
static char distr_name[] = "exponential";

#define sigma (params[0])
#define theta (params[1])
/*---------------------------------------------------------------------------*/

double
unur_pdf_exponential( double x, double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    /* standardize */
    x = (x-theta) / sigma;
  case 0:  /* standard */
    return ( (x<0.) ? 0. : exp(-x) );

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_pdf_exponential() */

/*---------------------------------------------------------------------------*/
  
double
unur_dpdf_exponential( double x, double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    return ( (x<theta) ? 0. : -exp( -(x-theta)/sigma ) / sigma);
  case 0:  /* standard */
    return ( (x<0.) ? 0. : -exp(-x) );

  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_dpdf_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_cdf_exponential( double x, double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    /* standardize */
    x = (x-theta) / sigma;
  case 0:  /* standard */
    return ( (x<0.) ? 0. : 1.-exp(-x) );
    
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_cdf_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_area_exponential( double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    return sigma;
  case 0:  /* standard */
    return 1.;
    
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_area_exponential() */

/*---------------------------------------------------------------------------*/

double
unur_mode_exponential( double *params, int n_params )
{
  switch (n_params) {

  case 2:  /* non standard */
    return theta;
  case 0:  /* standard */
    return 0.;
    
  default:
    _unur_error(distr_name,UNUR_ERR_NPARAM,"");
    return 0.;
  }

} /* end of unur_mode_exponential() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_exponential( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0 || n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  if (n_params > 0)
    CHECK_NULL(params,RETURN_NULL);

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set distribution id */
  distr->id = UNUR_DISTR_EXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* functions */
  DISTR.pdf  = unur_pdf_exponential;  /* pointer to p.d.f.               */
  DISTR.dpdf = unur_dpdf_exponential; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = unur_cdf_exponential;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.params[0] = 1.;        /* default for sigma */
  DISTR.params[1] = 0.;        /* default for theta */
  
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.params[1] = theta;
  case 1:
    DISTR.params[0] = sigma;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameter sigma */
  if (DISTR.params[0] <= 0.) {
    _unur_error(distr_name,UNUR_ERR_DISTR,"scale parameter sigma <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  DISTR.mode = unur_mode_exponential(DISTR.params,DISTR.n_params);
  DISTR.area = unur_area_exponential(DISTR.params,DISTR.n_params);

  /* domain */
  DISTR.domain[0] = DISTR.params[1]; /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_MODE   |
		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_exponential() */

/*---------------------------------------------------------------------------*/
#undef sigma 
#undef theta 
/*---------------------------------------------------------------------------*/
