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
#ifndef __UNURAN_CONFIG_H_SEEN
#define __UNURAN_CONFIG_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

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

#else

#define UNUR_DEBUG 0

#endif

/*            UNUR_DB_INFO      | \ */

/*---------------------------------------------------------------------------*/
/* Set default flag for debugging of generators:                             */
/*    0  ... do not write any infos about generator.                         */
/*    1  ... write only pameters and resulting structure of generator.       */
/*    ~0 ... show all available information.                                 */
/* (Detailed discription of possible flags (besides 1) in generator files.)  */

#define UNUR_DEBUGFLAG_DEFAULT   (~0u)

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
/* UNUR_URNG_POINTER:                                                        */
/*     Use a pointer to the routine without an argment                       */
/*     (e.g.: "double (*urng)(void)")                                        */
/*                                                                           */
/* UNUR_URNG_PRNG:                                                           */
/*     Use a pointer to a routine from the plab library                      */
/*     (see http://random.mat.sbg.ac.at/ftp/pub/software/gen/)               */
/*                                                                           */

#define UNUR_URNG_INVOKE UNUR_URNG_POINTER
/*  #define UNUR_URNG_INVOKE UNUR_URNG_PRNG */

/*---------------------------------------------------------------------------*/
/* Default name of uniform random number generator.                          */

#if UNUR_URNG_INVOKE == UNUR_URNG_POINTER

/* valid name of a C routine                                                 */
#define UNUR_URNG_DEFAULT uniform

#elif UNUR_URNG_INVOKE == UNUR_URNG_PRNG

/* valid parameter (char) string for prng-2.2                                */
#define UNUR_URNG_DEFAULT "LCG(2147483647,950706376,0,1)"

#else

#error UNUR_URNG_INVOKE not valid !!

#endif  /* UNUR_URNG_INVOKE */

/*---------------------------------------------------------------------------*/

/** TODO: we have to find a better place for these macros                    */
#define UNUR_MALLOC_SIZE   10
#define HAVE_GETTIMEOFDAY

#define UNUR_DISTR_MAXPARAMS  5 /* maximal number of parameters for distribution */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/



