/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_urng.h                                                       *
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
#ifndef __UNUR_URNG_H_SEEN
#define __UNUR_URNG_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/*---------------------------------------------------------------------------*/
/* uniform random number generator                                           */

/* We have to define the following macros:
   (see also unuran.h)

   UNUR_URNG_DEFAULT
      ... name|pointer of default urng (depends on UNUR_URNG_INVOKE)
          to be set in unuran_config.h

   _unur_call_urng(gen)
      ... function call to urng (via struct unur_gen)
*/

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER
/*---------------------------------------------------------------------------*/

/* prototype for uniform rng  */
double UNUR_URNG_DEFAULT(void);

/* type of uniform random number generator                                   */
typedef double (*UNUR_URNG_TYPE)(void);

/* function call to uniform rng */
#define _unur_call_urng(gen)        ((*(gen->urng))())

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* header file from prng library */
#include <prng.h>

/* type of uniform random number generator                                   */
typedef struct prng *UNUR_URNG_TYPE;

/* function call to uniform rng */
#define _unur_call_urng(gen)        (prng_get_next(gen->urng))

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_INVOKE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_INVOKE */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_URNG_H_SEEN */
/*---------------------------------------------------------------------------*/






