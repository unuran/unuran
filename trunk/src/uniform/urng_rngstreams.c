/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_rngstreams.c                                                 *
 *                                                                           *
 *   routines to get new URNG object with sampling routine of type           *
 *   RNGSTREAMSPRNG (Pierre L'Ecuyer's RNGSTREAMS package).                  *
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
#include "urng_rngstreams.h"
/*---------------------------------------------------------------------------*/
#if defined(UNURAN_HAS_RNGSTREAMS) && UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
#ifndef HAVE_LIBRNGSTREAMS
# error
# error +------------------------------------------------------------+
# error ! You have defined UNURAN_HAS_RNGSTREAMS in unuran_config.h  +
# error ! but Pierre L`Ecuyer`s RNGSTREAMS library is not installed. +
# error +------------------------------------------------------------+
# error
#endif
/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_urng_rngstreamptr_new( RngStream rngstream )
     /*----------------------------------------------------------------------*/
     /* get new URNG object of type RNGSTREAMS.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rngstream ... pointer to generator structure                       */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;

  /* check argument */
  if (rngstream == NULL) {
    _unur_error("URNG",UNUR_ERR_NULL,"Cannot create RNGSTREAM object");
    return NULL;
  }

  urng = unur_urng_new( (double(*)(void*)) RngStream_RandU01, rngstream );
  unur_urng_set_reset    (urng, (void(*)(void*)) RngStream_ResetStartStream);
  unur_urng_set_delete   (urng, (void(*)(void*)) RngStream_DeleteStream);
  unur_urng_set_anti     (urng, (void(*)(void*,int)) RngStream_SetAntithetic);
  unur_urng_set_nextsub  (urng, (void(*)(void*)) RngStream_ResetNextSubstream);
  unur_urng_set_resetsub (urng, (void(*)(void*)) RngStream_ResetStartSubstream);

  /* There only a function for seeding the RngStreams package, but no  */
  /* function for seeding an individual random stream.                 */

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
#endif /* defined(UNURAN_HAS_RNGSTREAMS) && UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
