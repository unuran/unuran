/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_prng.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *     Function prototypes for using uniform of type PRNG:                   *
 *     Otmar Lendl's prng package,                                           *
 *     see http://statistik.wu-wien.ac.at/prng/ or                           *
 *     http://random.mat.sbg.ac.at/.                                         *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
#ifndef URNG_PRNG_H_SEEN
#define URNG_PRNG_H_SEEN
/*---------------------------------------------------------------------------*/
#if defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN)
/*---------------------------------------------------------------------------*/
#include <prng.h>
/*---------------------------------------------------------------------------*/

/* 
   =NODE  URNG-PRNG  Interface to Otmar Lendl's pseudo-random number generators

   =UP URNG [40]

   =DESCRIPTION
      URNGs from the @code{prng} library. It provides a very
      flexible way to sample form arbitrary URNGs by means of an object
      oriented programing paradigma. Similarly to the UNURAN library
      independent generator objects can be build and used.
      
      This library has been developed by the pLab group at the university
      of Salzburg (Austria, EU) and implemented by Otmar Lendl.
      It is available from
      @uref{http://statistik.wu-wien.ac.at/prng/}
      or from the pLab site at
      @uref{http://random.mat.sbg.ac.at/}.

   =HOWTOUSE
      This library has to be installed before compiling UNURAN and
      UNURAN_HAS_PRNG has to be defined in @file{src/unuran_config.h}.
      Do not forget to link your executables against @file{libprng}.

      The following routines are supported for URNG objects of
      type PRNG:

      @itemize @minus
      @item unur_urng_sample()
      @item unur_urng_sample_array()
      @item unur_urng_reset() 
      @item unur_urng_free()
      @end itemize

   =END

*/

/*---------------------------------------------------------------------------*/

/* =ROUTINES */

/*---------------------------------------------------------------------------*/

UNUR_URNG *unur_urng_prng_new( const char *prngstr );
/*
   Make object for URNGs from Otmar Lendl's @file{prng} package. 
   @var{prngstr} is a string that contains the necessary information
   to create a uniform random number generator. For the format of this
   string see the @file{prng} user manual.

   The @file{prng} library provides a very flexible way to sample form
   arbitrary URNGs by means of an object oriented programing
   paradigma. Similarly to the UNURAN library independent generator
   objects can be build and used. The library has been developed
   and implemented by Otmar Lendl as member of the pLab group at the
   university of Salzburg (Austria, EU). 

   It is available via anonymous ftp from
   @uref{http://statistik.wu-wien.ac.at/prng/}
   or from the pLab site at
   @uref{http://random.mat.sbg.ac.at/}.
*/

UNUR_URNG *unur_urng_prngptr_new( struct prng *urng );
/*
   Similar to unur_urng_prng_new() but it uses a pointer to a
   generator object as returned by @code{prng_new(prngstr)};
   see @file{prng} manual for details.
*/

/* =END */

/*---------------------------------------------------------------------------*/
#endif  /* defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN) */
/*---------------------------------------------------------------------------*/
#endif  /* URNG_PRNG_H_SEEN */
/*---------------------------------------------------------------------------*/
