/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_misc.c                                                          *
 *                                                                           *
 *   miscelleanous routines                                                  *
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

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

/* global variable for default debugging flags */
unsigned _unur_default_debugflag = UNUR_DEBUGFLAG_DEFAULT;

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Set debuging flags                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_debug( struct unur_par *par, unsigned debug )
     /*----------------------------------------------------------------------*/
     /* set debugging flag for generator                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   debug ... debugging flag                                           */
     /*----------------------------------------------------------------------*/
{
  CHECK_NULL(par,0);

#ifdef UNUR_ENABLE_LOGGING
  par->debug = debug;
  return 1;
#else
  _unur_warning("DEBUG",_UNUR_ERR_COMPILE,"debugging, #define UNUR_ENABLE_LOGGING");
  return 0;
#endif

} /* end of unur_set_debug() */
  
/*---------------------------------------------------------------------------*/

int
unur_set_default_debug( struct unur_par *par, unsigned debug )
     /*----------------------------------------------------------------------*/
     /* set default debugging flag for generator                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   debug ... debugging flag                                           */
     /*----------------------------------------------------------------------*/
{
  CHECK_NULL(par,0);

  _unur_default_debugflag = debug;
  return 1;
} /* end of unur_set_default_debug() */
  
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**                                                                         **/
/**  Get data about generator object                                        **/
/**                                                                         **/
/*****************************************************************************/

int
unur_get_dimension( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get dimension of generator for multivariate distribution             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   dimension of distribution                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  switch (gen->method) {
  case UNUR_METH_RECT:
    COOKIE_CHECK(gen,CK_RECT_GEN,0);
    return gen->data.rect.dim;
    break;
  default:
    /* method unknown assume dim = 1 */
    /** TODO: distinguish between univariate --> dim = 1 and
	unknown multivariate **/
    return 1.;
  }

} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  misc                                                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

#define ARCMEAN_HARMONIC 1.e5  /* use harmonic mean when abs larger than this value */

double
_unur_arcmean( double x0, double x1 )
     /*----------------------------------------------------------------------*/
     /* compute "arctan mean" of two numbers.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x0, x1 ... two numbers                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   mean                                                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   "arctan mean" = tan(0.5*(arctan(x0)+arctan(x1)))                   */
     /*                                                                      */
     /*   a combination of arithmetical mean (for x0 and x1 close to 0)      */
     /*   and the harmonic mean (for |x0| and |x1| large).                   */
     /*----------------------------------------------------------------------*/
{
  /** TODO: possible over/underflow (?) **/

  /* we need x0 < x1 */
  if (x0>x1) {double tmp = x0; x0=x1; x1=tmp;}

  if (x1 < -ARCMEAN_HARMONIC || x0 > ARCMEAN_HARMONIC)
    /* use harmonic mean */
    return 2./(1./x0 + 1./x1);
  
  return tan( (((x0<=-INFINITY) ? -M_PI/2. : atan(x0)) + ((x1>=INFINITY) ? M_PI/2. : atan(x1))) / 2. );

} /* end of _unur_arcmean() */

/*---------------------------------------------------------------------------*/




