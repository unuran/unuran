/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_typedefs.h                                                   *
 *                                                                           *
 *   type UNUR_URNG for uniform random number generators                     *
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
#ifndef URNG_TYPEDEFS_H_SEEN
#define URNG_TYPEDEFS_H_SEEN
/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

typedef struct unur_urng_generic UNUR_URNG;

#ifdef UNURAN_HAS_PRNG
#  include <prng.h>
#endif

#ifdef UNURAN_HAS_RNGSTREAMS
#  include <RngStreams.h>
#endif

#ifdef UNURAN_HAS_GSL
#  include <gsl/gsl_rng.h>
#endif

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_FVOID
/*---------------------------------------------------------------------------*/

typedef double (UNUR_URNG)(void);

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

#include <prng.h>
typedef struct prng UNUR_URNG;

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAMS
/*---------------------------------------------------------------------------*/

#include <RngStreams.h>
typedef struct RngStream_InfoState UNUR_URNG;

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GSL
/*---------------------------------------------------------------------------*/

#include <gsl/gsl_rng.h>
typedef gsl_rng UNUR_URNG;

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/
#endif  /* #ifndef URNG_TYPEDEFS_H_SEEN */
/*---------------------------------------------------------------------------*/

