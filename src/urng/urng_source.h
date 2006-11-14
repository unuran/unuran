/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_source.h                                                     *
 *                                                                           *
 *   macros for calling and resetting uniform random number generators       *
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
#ifndef URNG_SOURCE_H_SEEN
#define URNG_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    ((urng)->sampleunif((urng)->state))

/* reset uniform RNG */
/* #define _unur_call_reset(urng) */
/* is defined in unur_API.c */

/*---------------------------------------------------------------------------*/
/* The following types are for backward compatibility.                       */
/* Its usage ist strongly deprecated!                                        */
/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (prng_get_next(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (prng_reset(urng),UNUR_SUCCESS)

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAMS
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (RngStream_RandU01(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (RngStream_ResetStartStream(urng),UNUR_SUCCESS)

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GSL
/*---------------------------------------------------------------------------*/

/* function call to uniform RNG */
#define _unur_call_urng(urng)    (gsl_rng_uniform_pos(urng))

/* reset uniform RNG */
#define _unur_call_reset(urng)   (UNUR_FAILURE)

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/
#endif  /* #ifndef URNG_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/

