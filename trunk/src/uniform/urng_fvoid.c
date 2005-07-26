/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.c                                                      *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type FVOID:    *
 *   double(*urng)(void), NULL) and global state variable.                   *
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
#include <unur_source.h>
#include "unur_uniform.h"
#include "urng.h"
#include "urng_fvoid.h"
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_fvoid_new( double (*random)(void), int (*reset)(void) )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type FVOID                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   random  ... pointer to uniform random number generator             */
     /*   reset   ... pointer to reset function for URNG                     */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_new( (double(*)(void*)) random, NULL );
  unur_urng_set_reset( urng, (void(*)(void*)) reset );
  return urng;
} /* end of unur_urng_fvoid_new() */

/*---------------------------------------------------------------------------*/
#endif   /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/

