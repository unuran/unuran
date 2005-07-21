/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_prng.c                                                       *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type PRNG.     *
 *   (Lendl's prng package, see http://statistik.wu-wien.ac.at/prng/ or      *
 *   http://random.mat.sbg.ac.at/.                                           *
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
#include "urng.h"

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
#ifdef UNURAN_HAS_PRNG
/*---------------------------------------------------------------------------*/
#ifndef HAVE_LIBPRNG
# error
# error +-----------------------------------------------------+
# error ! You have defined UNURAN_HAS_PRNG in unuran_config.h +
# error ! but Otmar Lendl`s PRNG library is not installed.    +
# error +-----------------------------------------------------+
# error
#endif
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_prngptr_new( struct prng *prng )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type PRNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prng ... pointer to generator structure                            */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_new( (double(*)(void*)) prng->get_next, prng );
  unur_urng_set_reset(urng, (void(*)(void*)) prng->reset);
  unur_urng_set_delete(urng, (void(*)(void*)) prng->destroy);
  return urng;
} /* end of unur_urng_prngptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_prng_new( const char *prngstr )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type PRNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prngstr ... string that describes generator                        */
     /*----------------------------------------------------------------------*/
{
  return unur_urng_prngptr_new(prng_new(prngstr));
} /* end of unur_urng_prng_new() */

/*---------------------------------------------------------------------------*/
#endif   /* #ifdef UNURAN_HAS_PRNG */
/*---------------------------------------------------------------------------*/
#endif   /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
