/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      unuran_config.h                                              *
 *                                                                           *
 *   compiler switches, compile time options and default values              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   author: Wolfgang.Hoermann @ statistik.wu-wien.ac.at                     *
 *           Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Wed Jul 28 14:25:52 1999                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifndef __UNURAN_CONFIG_H_SEEN
#define __UNURAN_CONFIG_H_SEEN
/*---------------------------------------------------------------------------*/

#include <config.h>

/*---------------------------------------------------------------------------*/
/* Debugging.                                                                 */
/*                                                                           */
/* Debugging flags are                                                       */
/*                                                                           */
/*   UNUR_DB_STDERR    ...    write warnings and errors on stderr            */
/*   UNUR_DB_LOG       ...    write warnings and infos into logfile          */
/*   UNUR_DB_COOKIES   ...    use magic cookies                              */
/*   UNUR_DB_CHECKNULL ...    check for NULL pointer                         */
/*   UNUR_DB_CHECKARGS ...    check arguments                                */
/*   UNUR_DB_INFO      ...    write info about generator into logfile        */

#if 1
#define UNUR_DEBUG (  \
          UNUR_DB_STDERR    | \
          UNUR_DB_LOG       | \
          UNUR_DB_COOKIES   | \
          UNUR_DB_CHECKNULL | \
          UNUR_DB_CHECKARGS | \
          UNUR_DB_INFO      | \
0 )
#endif

//#define UNUR_DEBUG 0

/*            UNUR_DB_INFO      | \ */

/*---------------------------------------------------------------------------*/
/* Set default flag for debugging of generators:                             */
/*    0  ... do not write any infos about generator.                         */
/*    1  ... write only pameters and resulting structure of generator.       */
/*    ~0 ... show all available information.                                 */
/* (Detailed discription of possible flags (besides 1) in generator files.)  */

#define UNUR_DEBUGFLAG_DEFAULT   (~0UL)

/*---------------------------------------------------------------------------*/
/* (Default) name of log file.                                               */
/* use "stdout" to write all infos to stdout.                                */

/*  #define UNUR_LOG_FILE "unuran.log" */
#define UNUR_LOG_FILE "stdout"

/*---------------------------------------------------------------------------*/
/* Invoke uniform random number generator.                                   */
/*                                                                           */
/* There are several ways to use a uniform pseudo random number generator:   */
/*                                                                           */
/* UNUR_URNG_LINKED:                                                         */
/*     Link a routine without an argument with a given name                  */
/*     (e.g.: "rand(void)")                                                  */
/*                                                                           */
/* UNUR_URNG_POINTER:                                                        */
/*     Use a pointer to the routine without an argment                       */
/*     (e.g.: "double (*urng)(void)")                                        */
/*                                                                           */
/* UNUR_URNG_PRNG:                                                           */
/*     Use a pointer to a routine from the plab library                      */
/*     (see http://random.mat.sbg.ac.at/ftp/pub/software/gen/)               */
/*                                                                           */

/*  #define UNUR_URNG_INVOKE UNUR_URNG_LINKED */
#define UNUR_URNG_INVOKE UNUR_URNG_POINTER
/*  #define UNUR_URNG_INVOKE UNUR_URNG_PRNG */

/*---------------------------------------------------------------------------*/
/* Default name of uniform random number generator.                          */

#if UNUR_URNG_INVOKE != UNUR_URNG_PRNG
/* valid name of a C routine (case LINKED and POINTER)                       */

#define UNUR_URNG_DEFAULT uniform

#else
/* valid parameter (char) string for prng-2.2 (case PRNG)                    */

#define UNUR_URNG_DEFAULT "LCG(2147483647,950706376,0,1)"

#endif

/*---------------------------------------------------------------------------*/

/** TODO: we have to find a better place for these macros                    */
#define UNUR_MALLOC_SIZE   10
#define HAVE_GETTIMEOFDAY

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/



