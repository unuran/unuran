/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      powerexpon.c                                                 *
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
 *  Power-exponential (Subbotin) distribution [3; ch.24, p.195]              *
 *                                                                           *
 *  pdf:     exp(-1/2 * abs((x-theta)/phi) ^ (2/delta) )                     *
 *  domain:  -infinity < x < infinity                                        *
 *                                                                           *
 *  parameters:                                                              *
 *     0:  delta > 0 ... shape                                               *
 *     1:  theta     ... location                                            *
 *     2:  phi > 0   ... scale                                               *
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
static char distr_name[] = "power-exponential";

#define delta (params[0])
#define theta (params[1])
#define phi   (params[2])
/*---------------------------------------------------------------------------*/

double
unur_pdf_powerexponential( double x, double *params, int n_params )
{ 
  register double z;
  z = (x - theta) / phi;
  return exp( - pow( abs(z), 2./delta ) * 0.5 );
} /* end of unur_pdf_powerexponential() */

/*---------------------------------------------------------------------------*/

double
unur_dpdf_powerexponential( double x, double *params, int n_params )
{
  register double z, tmp;

  z = (x - theta) / phi;

  if (z == 0.)    /* derivative is not defined, but ...        */
    return 0.;    /* a tangent parallel to x-axis is possible. */

  tmp = exp( - pow( abs(z), 2./delta ) * 0.5 ) * pow(abs(z),2./delta-1.) / (delta*phi);

  /* sign ! */
  return ( (z<0.) ? tmp : -tmp );
} /* end of unur_dpdf_powerexponential() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_powerexponential( double *params, int n_params )
{
#define DISTR distr->data.cont
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 1 || n_params > 3) {
    _unur_warning(distr_name,UNUR_ERR_GENERIC,"invalid number parameter");
    return NULL;
  }
  CHECK_NULL(params,RETURN_NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_POWEREXPONENTIAL;

  /* name of distribution */
  distr->name = distr_name;
                
  /* functions */
  DISTR.pdf  = unur_pdf_powerexponential;  /* pointer to p.d.f.               */
  DISTR.dpdf = unur_dpdf_powerexponential; /* pointer to derivative of p.d.f. */
  /* DISTR.cdf = unur_cdf_powerexponential;    pointer to c.d.f.               */

  /* default parameters */
  DISTR.params[1] = 0.;        /* default for theta */
  DISTR.params[2] = 1.;        /* default for phi   */

  /* copy parameters */
  DISTR.params[0] = delta;
  switch (n_params) {
  case 3:
    DISTR.params[2] = phi;
  case 1:
    DISTR.params[1] = theta;
  default:
    n_params = 3;
  }

  /* check parameter sigma */
  if (DISTR.params[0] <= 0. || DISTR.params[2] <= 0.) {
    _unur_error(distr_name ,UNUR_ERR_DISTR,"delta <= 0 or phi <= 0.");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* mode and area below p.d.f. */
  /* DISTR.mode = unur_mode_powerexponential(DISTR.params,DISTR.n_params); */
  /* DISTR.area = unur_area_powerexponential(DISTR.params,DISTR.n_params); */

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
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
} /* end of unur_distr_powerexponential() */

/*---------------------------------------------------------------------------*/
#undef theta
#undef delta
#undef phi  
/*---------------------------------------------------------------------------*/



