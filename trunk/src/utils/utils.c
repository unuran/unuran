/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      utils.c                                                      *
 *                                                                           *
 *   auxilliary subroutines                                                  *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Mon Oct 18 20:46:18 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_utils.h>
#include <unur_math.h>

/*---------------------------------------------------------------------------*/

#define ARCMEAN_HARMONIC 1.e5  /* use harmonic mean when abs larger than this value */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  misc                                                                   **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

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

