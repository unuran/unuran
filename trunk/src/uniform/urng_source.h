/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_source.h                                                     *
 *                                                                           *
 *   macros and prototypes for included uniform random number generators     *
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
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    ((urng)->sampleunif((urng)->params))

/* reset uniform RNG */
/* #define _unur_call_reset(urng) */
/* is defined in unur_API.c */

/* default generators */
#if defined(UNURAN_HAS_RNGSTREAMS)
/* use type RNGSTREAMS */
#  define UNUR_URNG_DEFAULT      (unur_urng_rngstream_new("URNG_main"))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_rngstream_new("URNG_aux"))

#elif defined(UNURAN_HAS_GSL)
/* use type GSL */
#  define UNUR_URNG_DEFAULT      (unur_urng_gsl_new(UNUR_URNG_DEFAULT_GSL))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_gsl_new(UNUR_URNG_AUX_DEFAULT_GSL))

#elif defined(UNURAN_HAS_PRNG)
/* use type PRNG */
#  define UNUR_URNG_DEFAULT      (unur_urng_prng_new(UNUR_URNG_DEFAULT_PRNG))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_prng_new(UNUR_URNG_AUX_DEFAULT_PRNG))

#else
/* use type FVOID */
#  define UNUR_URNG_DEFAULT      (unur_urng_fvoid_new(UNUR_URNG_DEFAULT_FVOID, \
					              UNUR_URNG_DEFAULT_RESET_FVOID))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_fvoid_new(UNUR_URNG_AUX_DEFAULT_FVOID, \
					              UNUR_URNG_AUX_DEFAULT_RESET_FVOID))
#endif

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_FVOID
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    ((*(urng))())

/* reset uniform RNG */
#define _unur_call_reset(urng) \
 ( (urng)==(unur_urng_MRG31k3p) ? (unur_urng_MRG31k3p_reset()) \
    : ( (urng)==(unur_urng_fish) ? (unur_urng_fish_reset()) \
      : ( (urng)==(unur_urng_mstd) ? (unur_urng_mstd_reset()) \
	  : (UNUR_FAILURE) )))

/* default generators */
#define UNUR_URNG_DEFAULT        UNUR_URNG_DEFAULT_FVOID
#define UNUR_URNG_AUX_DEFAULT    UNUR_URNG_AUX_DEFAULT_FVOID

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (prng_get_next(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (prng_reset(urng),UNUR_SUCCESS)

/* default generators */
#define UNUR_URNG_DEFAULT        (prng_new(UNUR_URNG_DEFAULT_PRNG))
#define UNUR_URNG_AUX_DEFAULT    (prng_new(UNUR_URNG_AUX_DEFAULT_PRNG))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAMS
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (RngStream_RandU01(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (RngStream_ResetStartStream(urng),UNUR_SUCCESS)

/* default generators */
#define UNUR_URNG_DEFAULT        (RngStream_CreateStream("URNG_main"))
#define UNUR_URNG_AUX_DEFAULT    (RngStream_CreateStream("URNG_aux"))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GSL
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (gsl_rng_uniform_pos(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (UNUR_FAILURE)

/* default generators */
#define UNUR_URNG_DEFAULT        (gsl_rng_alloc(UNUR_URNG_DEFAULT_GSL))
#define UNUR_URNG_AUX_DEFAULT    (gsl_rng_alloc(UNUR_URNG_AUX_DEFAULT_GSL))

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/

