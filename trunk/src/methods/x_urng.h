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

  typedef double (UNUR_URNG)(void);
  

  2. UNUR_URNG_TYPE == UNUR_URNG_PRNG

*/

/*---------------------------------------------------------------------------*/
/* uniform random number generator                                           */

/* We have to define the following macros:

   UNUR_URNG_DEFAULT
      ... name|pointer of default urng (depends on UNUR_URNG_TYPE)
          to be set in unuran_config.h

   _unur_call_urng(gen)
      ... function call to urng (via UNUR_GEN)
*/

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_POINTER
/*---------------------------------------------------------------------------*/

/* prototype for uniform rng  */
double UNUR_URNG_DEFAULT(void);

/* type of uniform random number generator                                   */
typedef double (UNUR_URNG)(void);

/* function call to uniform rng */
#define _unur_call_urng(gen)        ((*(gen->urng))())

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* header file from prng library */
#include <prng.h>

/* type of uniform random number generator                                   */
typedef struct prng UNUR_URNG;

/* function call to uniform rng */
#define _unur_call_urng(gen)        (prng_get_next(gen->urng))

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* get and set default uniform RNG                                           */
/* (defined in src/utils/urng.c)                                             */

UNUR_URNG *unur_get_default_urng( void );
UNUR_URNG *unur_set_default_urng( UNUR_URNG *urng_new );

/*---------------------------------------------------------------------------*/
/* set, get or change uniform RNG for generator                              */

int unur_set_urng( UNUR_PAR *parameters, UNUR_URNG *urng );
UNUR_URNG *unur_chg_urng( UNUR_GEN *generator, UNUR_URNG *urng );
UNUR_URNG *unur_get_urng( UNUR_GEN *generator );

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __X_URNG_H_SEEN */
/*---------------------------------------------------------------------------*/
