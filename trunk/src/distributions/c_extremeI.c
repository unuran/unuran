/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_extremeI.c                                                 *
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
 *  distr: Extreme value type I distribution  [3; ch.22, p.2]                *
 *  (also Gumbel-type distribution)                                          *
 *                                                                           *
 *                                                                           *
 *  Type I (also Gumbel-type distribution)                                   *
 *                                                                           *
 *  cdf:       F(x) = exp( -exp( -(x-zeta)/theta ) )                         *
 *  pdf:       f(x) = exp( -exp( -(x-zeta)/theta ) - (x-zeta)/theta )        *               
 *  domain:    -infinity < x <infinity                                       *
 *  constant:  1/theta                                                       *
 *                                                                           *
 *  parameters: 2                                                            *
 *     0:  zeta      (0)  ... location                                       *
 *     1:  theta >0  (1)  ... scale                                          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Standard Form:                                                            *
 *                                                                           *
 *  cdf:       F(x) = exp( -exp(-x) )                                        *
 *  pdf:       f(x) = exp( -exp(-x) - x )                                    *
 *  domain:    -infinity < x <infinity                                       *
 *  constant:  1                                                             *
 *                                                                           *
 *  parameters: none                                                         *
 *                                                                           *
 *     0:  zeta  = 0                                                         *
 *     1:  theta = 1                                                         *
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

static const char distr_name[] = "extremeI";

/* parameters */
#define zeta   params[0]    /* location */
#define theta  params[1]    /* scale */

#define DISTR distr->data.cont

/* function prototypes                                                       */
static double _unur_pdf_extremeI(double x, UNUR_DISTR *distr);
static double _unur_dpdf_extremeI(double x, UNUR_DISTR *distr);
static double _unur_cdf_extremeI(double x, UNUR_DISTR *distr);

/*---------------------------------------------------------------------------*/

double
_unur_pdf_extremeI( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - zeta) / theta;
  case 0: default: /* standard */
    return ( exp( -exp(-x) - x ) / theta );
  }
} /* end of _unur_pdf_extremeI() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_extremeI( double x, UNUR_DISTR *distr )
{ 
  register double factor = 1.;
  register double expx;
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:  /* non standard */
    /* standardize */
    factor = 1. / (theta * theta);
    x = (x - zeta) / theta;
  case 0: default: /* standard */
    expx = exp(-x);
    return ( exp( -expx + x ) * (expx - 1.) * factor );
  }
} /* end of unur_dpdf_extremeI() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_extremeI( double x, UNUR_DISTR *distr )
{ 
  register double *params = DISTR.params;

  switch (DISTR.n_params) {
  case 2:  /* non standard */
    /* standardize */
    x = (x - zeta) / theta;
  case 0: default: /* standard */
    return exp( -exp( -x) );
  }
} /* end of _unur_cdf_extremeI() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_extremeI( double *params, int n_params )
{
  register struct unur_distr *distr;

  /* check new parameter for generator */
  if (n_params < 0) n_params = 0;
  if (n_params > 2) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_NPARAMS,"too many");
    n_params = 2; }
  if (n_params > 0)
    CHECK_NULL(params,NULL);

  /* get new (empty) distribution object */
  distr = unur_distr_cont_new();

  /* set distribution id */
  distr->id = UNUR_DISTR_EXTREME_I;

  /* name of distribution */
  distr->name = distr_name;
                
  /* how to get special generators */
  DISTR.init = _unur_stdgen_extremeI_init;

  /* functions */
  DISTR.pdf  = _unur_pdf_extremeI;  /* pointer to p.d.f.               */
  DISTR.dpdf = _unur_dpdf_extremeI; /* pointer to derivative of p.d.f. */
  DISTR.cdf  = _unur_cdf_extremeI;  /* pointer to c.d.f.               */

  /* default parameters */
  DISTR.zeta  = 0.;
  DISTR.theta = 1.;
  
  /* copy parameters */
  switch (n_params) {
  case 2:
    DISTR.theta = theta;
  case 1:
    DISTR.zeta = zeta;
    n_params = 2;           /* number of parameters for non-standard form */
  default:
  }

  /* check parameters */
  if (DISTR.theta <= 0. ) {
    _unur_error(distr_name,UNUR_ERR_DISTR_DOMAIN,"theta <= 0");
    free( distr ); return NULL;
  }

  /* number of arguments */
  DISTR.n_params = n_params;

  /* normalization constant: not required */

  /* mode and area below p.d.f. */
  DISTR.mode = DISTR.zeta;
  DISTR.area = 1.;

  /* domain */
  DISTR.domain[0] = -INFINITY;       /* left boundary  */
  DISTR.domain[1] = INFINITY;        /* right boundary */

  /* indicate which parameters are set */
  distr->set = ( UNUR_DISTR_SET_PARAMS | 
		 UNUR_DISTR_SET_DOMAIN |
		 UNUR_DISTR_SET_STDDOMAIN |
 		 UNUR_DISTR_SET_MODE   |
  		 UNUR_DISTR_SET_PDFAREA );

  /* return pointer to object */
  return distr;

} /* end of unur_distr_extremeI() */

/*---------------------------------------------------------------------------*/
#undef c    
#undef alpha
#undef zeta 
#undef DISTR
/*---------------------------------------------------------------------------*/
