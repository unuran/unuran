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
#include "unur_uniform.h"
#include "urng.h"

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
#ifdef UNURAN_HAS_RNGSTREAMS
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_rngstreamptr_new( RngStream rngstream )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type RNGSTREAMS.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urngstr ... pointer to generator structure                         */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = unur_urng_new( (_unur_urng_doublevoidptr) RngStream_RandU01, rngstream );
  unur_urng_set_reset(urng, (_unur_urng_intvoidptr) RngStream_ResetStartStream);
  unur_urng_set_delete(urng, (_unur_urng_voidvoidptr) RngStream_DeleteStream);
  return urng;
} /* end of unur_urng_rngstreamptr_new() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_rngstream_new( const char *urngstr )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type RNGSTREAMS.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   prngstr ... string that describes generator                        */
     /*----------------------------------------------------------------------*/
{
  return unur_urng_rngstreamptr_new(RngStream_CreateStream(urngstr));
} /* end of unur_urng_prng_new() */

/*---------------------------------------------------------------------------*/
#endif   /* #ifdef UNURAN_HAS_RNGSTREAMS */
/*---------------------------------------------------------------------------*/
#endif   /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
