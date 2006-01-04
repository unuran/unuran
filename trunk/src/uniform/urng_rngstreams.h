/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using uniform of type RNGSTREAMSPRNG:         *
 *     Pierre L'Ecuyer's RNGSTREAMS package.                                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
#ifndef URNG_RNGSTREAMS_H_SEEN
#define URNG_RNGSTREAMS_H_SEEN
/*---------------------------------------------------------------------------*/
#if defined(UNURAN_HAS_RNGSTREAMS) && UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/
#include <RngStream.h>
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-RNGSTREAMS  Interface to L'Ecuyer's RNGSTREAMS random number generators

   =UP URNG [50]

   =DESCRIPTION
      URNGs from Pierre L'Ecuyer's @code{RngStream} library for multiple 
      independent streams of pseudo-random numbers. 
      It allows to split a randpm stream into many substreams.
      It is available from
      A GNU-style package is available from
      @uref{http://statistik.wu-wien.ac.at/software/RngStreams/}.

   =HOWTOUSE
      This library has to be installed before compiling UNURAN and
      UNURAN_HAS_RNGSTREAMS has to be defined in @file{src/unuran_config.h}.
      Do not forget to link your executables against this library.

      The following routines are supported for URNG objects of
      type RANDOMSTREAMS:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_reset() 
      @item unur_urng_nextsub() 
      @item unur_urng_resetsub() 
      @item unur_urng_anti() 
      @item unur_urng_free()
      @end itemize

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_rngstream_new( const char *urngstr );
/*
   Make object for URNGs from Pierre L'Ecuyer's @file{RngStream}
   library. This library provides multiple independent streams of
   pseudo-random numbers and is available from
   @uref{http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/}.
   @var{urngstr} is an arbitrary string to label a stream. It need not
   be unique.
*/

UNUR_URNG *unur_urng_rngstreamptr_new( RngStream rngstream );
/*
   Similar to unur_urng_rngstream_new() but it uses a pointer to a 
   generator object as returned by @code{RngStream_CreateStream()}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* defined(UNURAN_HAS_RNGSTREAMS) && UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
#endif  /* URNG_RNGSTREAMS_H_SEEN */
/*---------------------------------------------------------------------------*/
