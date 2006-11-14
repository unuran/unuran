/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_default_source.h                                             *
 *                                                                           *
 *   macros for default uniform random number generators                     *
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
#ifndef URNG_DEFAULT_SOURCE_H_SEEN
#define URNG_DEFAULT_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

#if UNUR_URNG_DEFAULT_TYPE == UNUR_URNG_PRNG
/* use type PRNG */
#  ifndef UNURAN_HAS_PRNG
#    error Choosen default URNG type requires PRNG library
#  else 
#    define UNUR_URNG_DEFAULT      (unur_urng_prng_new(UNUR_URNG_DEFAULT_PRNG))
#    define UNUR_URNG_AUX_DEFAULT  (unur_urng_prng_new(UNUR_URNG_AUX_DEFAULT_PRNG))
#  endif

#elif UNUR_URNG_DEFAULT_TYPE == UNUR_URNG_RNGSTREAMS
/* use type RNGSTREAMS */
#  ifndef UNURAN_HAS_RNGSTREAMS
#    error Choosen default URNG type requires RNGSTREAMS library
#  else 
#    define UNUR_URNG_DEFAULT      (unur_urng_rngstream_new("URNG_main"))
#    define UNUR_URNG_AUX_DEFAULT  (unur_urng_rngstream_new("URNG_aux"))
#  endif

#elif UNUR_URNG_DEFAULT_TYPE == UNUR_URNG_GSL
/* use type GSL */
#  ifndef UNURAN_HAS_GSL
#    error Choosen default URNG type requires GSL library
#  else 
#    define UNUR_URNG_DEFAULT      (unur_urng_gsl_new(UNUR_URNG_DEFAULT_GSL))
#    define UNUR_URNG_AUX_DEFAULT  (unur_urng_gsl_new(UNUR_URNG_AUX_DEFAULT_GSL))
#  endif

#elif UNUR_URNG_DEFAULT_TYPE == UNUR_URNG_FVOID
/* use type FVOID */
#  define UNUR_URNG_DEFAULT      (unur_urng_fvoid_new(UNUR_URNG_DEFAULT_FVOID, \
					              UNUR_URNG_DEFAULT_RESET_FVOID))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_fvoid_new(UNUR_URNG_AUX_DEFAULT_FVOID, \
					              UNUR_URNG_AUX_DEFAULT_RESET_FVOID))
#else
#error UNUR_URNG_DEFAULT_TYPE not valid !!
#endif

/*---------------------------------------------------------------------------*/
/* The following types are for backward compatibility.                       */
/* Its usage ist strongly deprecated!                                        */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_FVOID
/*---------------------------------------------------------------------------*/

#define UNUR_URNG_DEFAULT        UNUR_URNG_DEFAULT_FVOID
#define UNUR_URNG_AUX_DEFAULT    UNUR_URNG_AUX_DEFAULT_FVOID

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

#define UNUR_URNG_DEFAULT        (prng_new(UNUR_URNG_DEFAULT_PRNG))
#define UNUR_URNG_AUX_DEFAULT    (prng_new(UNUR_URNG_AUX_DEFAULT_PRNG))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAMS
/*---------------------------------------------------------------------------*/

#define UNUR_URNG_DEFAULT        (RngStream_CreateStream("URNG_main"))
#define UNUR_URNG_AUX_DEFAULT    (RngStream_CreateStream("URNG_aux"))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GSL
/*---------------------------------------------------------------------------*/

#define UNUR_URNG_DEFAULT        (gsl_rng_alloc(UNUR_URNG_DEFAULT_GSL))
#define UNUR_URNG_AUX_DEFAULT    (gsl_rng_alloc(UNUR_URNG_AUX_DEFAULT_GSL))

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/
#endif  /* #ifndef URNG_DEFAULT_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/

