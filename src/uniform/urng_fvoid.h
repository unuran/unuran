/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_fvoid.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using uniform of type FVOID                   *
 *     (i.e. routine without an argment:   double uniform(void)              *
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
#ifndef URNG_FVOID_H_SEEN
#define URNG_FVOID_H_SEEN
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-FVOID  Simple interface for uniform random number generators

   =UP URNG [10]

   =DESCRIPTION
      Simple interface for URNGs of type @code{FVOID}, i.e.,
      routines without an argment:  @code{double uniform(void)}.

      If independent versions of the same URNG should be used, a copy of
      the subroutine has to be implement in the program code (with
      different names, of course).
      UNURAN contains some build-in URNGs of this type in directory
      @file{src/uniform/}.
      
   =HOWTOUSE
      Create an URNG object using unur_urng_fvoid_new(). 
      By this call a pointer to the sampling routine and (optional) a
      pointer to a reset routine are copied into the URNG object.
      Other functions, like seeding the URNG, switching to antithetic
      random number, or jumping to next substream, can be added to the
      URNG object by the respective calls, e.g. by
      unur_urng_set_seed().

      The following routines are supported for URNG objects of
      type FVOID:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_seed()   [optional]
      @item unur_urng_reset()   [optional]
      @item unur_urng_free()
      @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

UNUR_URNG *unur_urng_fvoid_new( double (*random)(void), int (*reset)(void) );
/*
   Make a URNG object for a genertor that consists of a single
   function call with a global state variable.

   @emph{Notice:} If independent versions of the same URNG should be
   used, copies of the subroutine with different names has to be
   implement in the program code.

   If there is no reset function use NULL for the second argument.
   
   UNURAN contains some build-in URNGs of this type in directory
   @file{src/uniform/}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* #if UNUR_URNG_TYPE == UNUR_URNG_GENERIC */
/*---------------------------------------------------------------------------*/
#endif  /* URNG_FVOID_H_SEEN */
/*---------------------------------------------------------------------------*/
