/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_urng.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros and function prototypes for using uniform         *
 *         random number generators inside UNURAN.                           *
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
#ifndef __X_URNG_H_SEEN
#define __X_URNG_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/*
  Using uniform random number generators.

  Each generator has a pointer to a uniform (pseudo-) random number
  generator (URNG). It can be set via the unur_set_urng() call. It is
  also possible change this pointer via unur_get_urng() or change the
  URNG for an existing generator object by means of unur_get_urng();

  By this very flexible concept it is possible that each generator has
  its own (independent) URNG or several generators can share the same
  URNG. 

  If no URNG is provided for a parameter or generator object a default
  generator is used which is the same for all generators. This URNG is
  defined in unuran_config.h at compile time. A pointer to
  this default URNG can be obtained via
  unur_get_default_urng(). Nevertheless it is also possible to
  override this default URNG by another one by means of the
  unur_set_default_urng() call. However this only takes effect for new
  parameter objects.

  The pointer to a URNG is of type @code{UNUR_URNG*}. Its definition 
  depends on the compilation switch @code{UNUR_URNG_TYPE} in unuran_config.h. 
  Currently we have two possible switches (other values would result
  in a compilation error):

  1. UNUR_URNG_TYPE == UNUR_URNG_POINTER

  This uses URNGs of type @code{double uniform(void)}.
  If independent versions of the same URNG should be used, a copy of
  the subroutine has to be implement in the program code (with
  different names, of course).

  2. UNUR_URNG_TYPE == UNUR_URNG_PRNG

  This uses the URNGs from the prng library. It provides a very
  flexible way to sample form arbitrary URNGs by means of an object
  oriented programing paradigma. Similarly to the UNURAN library
  independent generator objects can be build and used.
  Here @code{UNUR_URNG*} is simply a pointer to such a uniform
  generator object.

  This library has been developed by the pLab group at the university
  of Salzburg (Austria, EU) and implemented by Otmar Lendl.
  It is available via anonymous ftp from
  http://random.mat.sbg.ac.at/ftp/pub/software/gen/.

  It is possible to use other interfaces to URNGs without much
  troubles. If you need such a new interface please email the authors
  of the UNURAN library.

*/

/*---------------------------------------------------------------------------*/
/* uniform random number generator                                           */

/* We have to define the following macros:

   UNUR_URNG_DEFAULT
      ... name|pointer of default urng (depends on UNUR_URNG_TYPE)
          to be set in unuran_config.h

   _unur_call_urng(urng)
      ... function call to urng 
*/
/* Remark: UNUR_URNG_DEFAULT and _unur_call_urng() should be defined in      */
/*         source_urng.h. However this would be a confusing splitting of     */
/*         this few lines of code.                                           */

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER
/*---------------------------------------------------------------------------*/

/* prototype for uniform rng  */
double UNUR_URNG_DEFAULT(void);

/* type of uniform random number generator                                   */
typedef double (UNUR_URNG)(void);

/* function call to uniform rng */
#define _unur_call_urng(urng)        ((*(urng))())

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* header file from prng library */
#include <prng.h>

/* type of uniform random number generator                                   */
typedef struct prng UNUR_URNG;

/* function call to uniform rng */
#define _unur_call_urng(urng)        (prng_get_next(urng))

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* set, get or change uniform RNG for generator                              */

/*
  URNG for generator objects.
*/

int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
/*
  Use the URNG @code{urng} for the new generator. This overrides the
  default URNG.
*/

UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
/*
  Change the URNG for the given generator. It returns the pointer to
  the old URNG that has been used by the generator.
*/

UNUR_URNG *unur_get_urng( UNUR_GEN *generator );
/*
  Get the pointer to the URNG that is used by the generator.
  This is usefull if two generators should share the same URNG.
*/

/*---------------------------------------------------------------------------*/
/* get and set default uniform RNG                                           */
/* (defined in src/utils/urng.c)                                             */

/*
  Default URNG.
*/

UNUR_URNG *unur_get_default_urng( void );
/*
  Get the pointer to the default URNG. The default URNG is used by all
  generators where no URNG was set explicitly by a unur_set_urng()
  call.
*/

UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );
/*
  Change the default URNG for new generator objects. 
*/

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __X_URNG_H_SEEN */
/*---------------------------------------------------------------------------*/
