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
/* Invoke uniform random number generator.                                   */
/*                                                                           */
/* There are several ways to use a uniform pseudo random number generator:   */
/*                                                                           */
/* UNUR_URNG_POINTER:                                                        */
/*     Use a pointer to the routine without an argment, i.e.                 */
/*        double (*urng)(void);                                              */
/*                                                                           */
/* UNUR_URNG_PRNG:                                                           */
/*     Use a pointer to a routine from the plab library                      */
/*     (see http://random.mat.sbg.ac.at/ftp/pub/software/gen/)               */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/* define the possible compiler switches                                     */
#define UNUR_URNG_POINTER  2     /* use a pointer to a routine               */
#define UNUR_URNG_PRNG     3     /* use a pointer to gen. from prng-3.0      */
/*---------------------------------------------------------------------------*/

/* set type of uniform generator                                             */
/*  #define UNUR_URNG_TYPE UNUR_URNG_POINTER */
#define UNUR_URNG_TYPE UNUR_URNG_PRNG

/*---------------------------------------------------------------------------*/
/* Default name of uniform random number generator.                          */

#if UNUR_URNG_TYPE == UNUR_URNG_POINTER

/* valid name of a C routine                                                 */
#define UNUR_URNG_DEFAULT uniform
#define UNUR_URNG_AUX_DEFAULT uniform    /** TODO **/

#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG

/* valid parameter (char) string for prng-2.2                                */
#define UNUR_URNG_DEFAULT "LCG(2147483647,950706376,0,1)"
#define UNUR_URNG_AUX_DEFAULT "LCG(2147483647,16807,0,1)"

#else
#error UNUR_URNG_TYPE not valid !!
#endif  /* UNUR_URNG_TYPE */


/*---------------------------------------------------------------------------*/
/* Default name of log file.                                                 */
/*                                                                           */
/* Use "stdout" to write all infos to stdout.                                */
/* If no log file should be used, #undef the macro UNUR_ENABLE_LOGGING below.*/

/*  #define UNUR_LOG_FILE "unuran.log" */
#define UNUR_LOG_FILE "stdout"

/*---------------------------------------------------------------------------*/
/* Define this macro to switch on writing information about the              */
/* generator into log file.                                                  */
/*                                                                           */
/* If no log file should be used at all, #undef this macro.                  */

#define UNUR_ENABLE_LOGGING 1

/*---------------------------------------------------------------------------*/
/* Every message about a generator that is written into the log file or      */
/* or to the stderr starts with an identifier. Every generator has its own   */
/* unique identifier (well if there are not more than 999 generators).       */
/* It is composed by the generator type, followed by a dot and three digits. */
/* Building such a generator id can be disabled by undefining the following  */
/* macro. If it is disabled the generator type (a constant string) is used   */
/* for the identifier (which is not unique any more if more than one         */
/* generator of the same type is created).                                   */

#define UNUR_ENABLE_GENID  1

/*---------------------------------------------------------------------------*/
/* Set default flag for debugging of generators:                             */
/*                                                                           */
/*   UNUR_DEBUG_OFF    ... switch off debugging information                  */
/*   UNUR_DEBUG_INIT   ... pameters and structure of generator only          */
/*   UNUR_DEBUG_SETUP  ... information for setup step                        */
/*   UNUR_DEBUG_ADAPT  ... trace adaptive steps                              */ 
/*   UNUR_DEBUG_SAMPLE ... trace sampling                                    */
/*   UNUR_DEBUG_ALL    ... write all available debugging information         */
/*                                                                           */
/* Detailed discription of possible flags in file `./methods/x_debug.h'      */
/*                                                                           */
/* Debugging information is written into the log file.                       */
/* It only works if additionally UNUR_ENABLE_LOGGING is defined (see above). */

#define UNUR_DEBUGFLAG_DEFAULT   UNUR_DEBUG_ALL


/*---------------------------------------------------------------------------*/
/* Warnings and error messages.                                              */
/*                                                                           */
/* UNURAN produces a lot of (more or less useful) warnings and error         */
/* messages. These three compiler switches controll their output.            */

/* Enable warnings and error messages.                                       */
/* #undef this macro to suppress warnings and error messages.                */

#define UNUR_WARNINGS_ON  1

/* Write warnings and error messages into log file.                          */
/* It only works if additionally UNUR_ENABLE_LOGGING is defined (see above). */
/* #undef this macro to suppress output into log file.                       */

#define UNUR_ENABLE_LOGFILE  1

/* Write warnings and error messages to stderr.                              */
/* #undef this macro to suppress output on stderr.                           */
/*  #define UNUR_ENABLE_STDERR  1 */

/* Notice that if neither UNUR_ENABLE_LOGFILE nor UNUR_ENABLE_STDERR is      */
/* defined, then there are warnigs and error messages at all.                */
/* However then it is recommend not to define UNUR_WARNINGS_ON, either.      */

/*---------------------------------------------------------------------------*/
/* Check for invalide NULL pointer.                                          */
/*                                                                           */
/* UNURAN expects that every generator object is not NULL. Thus the user     */
/* is responsible to check the result of each unur_init() call!!             */
/*                                                                           */
/* However checking against an invalid NULL pointer can be switched on for   */
/* each pointer that occurs by defining  UNUR_ENABLE_CHECKNULL.              */

#define UNUR_ENABLE_CHECKNULL 1

/*---------------------------------------------------------------------------*/
/* Debugging tools.                                                          */
/* (for development only. there is no need to set these flags unless         */
/* changes are made in the library.)                                         */

/* use magic cookies to validate type of pointer */
#define UNUR_COOKIES  1


/*---------------------------------------------------------------------------*/
/* Precision for floating point arithmetic                                   */

/* Epsilon for comparision of two doubles:                                   */
/* Two doubles are considered equal if there relative difference is          */
/* less than UNUR_EPSILON (see file `./methods/source_fp.h' for details).    */

#define UNUR_EPSILON  1e-14    /* should be about 100 * DBL_EPSILON          */

/* Square root of machine epsilon. It is used to compare two doubles         */
/* when round-off errors have to be considered                               */
/* (see file `./methods/source_fp.h' for details).                           */

#define UNUR_SQRT_DBL_EPSILON  FLT_EPSILON  /* should be about sqrt(DBL_EPS) */


/*---------------------------------------------------------------------------*/
/* Distribution objects.                                                     */

/* Maximal number of parameters for the PDF of a distribution.               */
/* (It must be at least 5!)                                                  */
#define UNUR_DISTR_MAXPARAMS  5

/* Maximal size of automatically created probability vectors.                */
#define UNUR_MAX_AUTO_PV    100000


/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/


