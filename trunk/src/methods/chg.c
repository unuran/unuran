/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      chg.h                                                        *
 *                                                                           *
 *   change parameters for (existing) generators                             *
 *                                                                           *
 *   PARAMETER: struct unur_gen *                                            *
 *                                                                           *
 *   return:                                                                 *
 *     old parameter                                                         *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Fri Dec 17 09:55:33 CET 1999                         *
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

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for the distribution and its p.d.f.                         **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of discrete distributions                    **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of continuous distributions                  **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for all generators                                          **/
/**                                                                         **/
/*****************************************************************************/

UNUR_URNG_TYPE
unur_chg_urng( struct unur_gen *gen, UNUR_URNG_TYPE urng )
/*---------------------------------------------------------------------------*/
/* set uniform random number generator                                       */
/*                                                                           */
/* parameters:                                                               */
/*   gen     ... pointer to generator object                                 */
/*   urng    ... pointer to uniform random number generator                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
{
  UNUR_URNG_TYPE urng_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng,NULL);

  urng_old = gen->urng;

  gen->urng = urng;

  return urng_old;
} /* end of unur_chg_urng() */

/*---------------------------------------------------------------------------*/

