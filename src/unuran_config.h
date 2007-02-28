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
#ifndef UNURAN_CONFIG_H_SEEN
#define UNURAN_CONFIG_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *  Logging and debugging.                                                   *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Default name of log file.                                                 */
/*                                                                           */
/* If this macro is not defined all infos are written on 'stdout'.           */
/*                                                                           */
/* Remark: If no log file should be used,                                    */
/*         #undef the macro UNUR_ENABLE_LOGGING below.                       */

#define UNUR_LOG_FILE "unuran.log"

/*---------------------------------------------------------------------------*/
/* Define this macro to switch on writing information about the              */
/* generator into log file.                                                  */
/*                                                                           */
/* If no log file should be used at all, #undef this macro.                  */

#define UNUR_ENABLE_LOGGING 1

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
/* Detailed discription of possible flags in file `./utils/debug.h'          */
/*                                                                           */
/* Debugging information is written into the log file.                       */
/* It only works if additionally UNUR_ENABLE_LOGGING is defined (see above). */

#define UNUR_DEBUGFLAG_DEFAULT   UNUR_DEBUG_INIT

/*---------------------------------------------------------------------------*/
/* Warnings and error messages.                                              */
/*                                                                           */
/* UNURAN produces a lot of (more or less useful) warnings and error         */
/* messages. The following three compiler switches controll their output.    */

/* Enable warnings and error messages.                                       */
/* #undef this macro to suppress warnings and error messages.                */

#define UNUR_WARNINGS_ON  1

/* Write warnings and error messages into log file.                          */
/* It only works if additionally UNUR_ENABLE_LOGGING is defined (see above). */
/* #undef this macro to suppress output into log file.                       */

#define UNUR_ENABLE_LOGFILE  1

/* Notice that if UNUR_ENABLE_LOGFILE is not defined, then there are         */
/* no warnings and error messages.                                           */
/* However, then it is recommend _not_ to define UNUR_WARNINGS_ON, either.   */

/*---------------------------------------------------------------------------*/
/* Check for invalid NULL pointer.                                           */
/*                                                                           */
/* UNURAN expects that every generator object is not a NULL pointer.         */
/* Thus the user is responsible to check the result of each unur_init()      */
/* call!!                                                                    */
/*                                                                           */
/* However, checking against an invalid NULL pointer can be switched on      */
/* for each pointer that occurs by defining  UNUR_ENABLE_CHECKNULL.          */

#define UNUR_ENABLE_CHECKNULL 1

/*---------------------------------------------------------------------------*/
/* Debugging tools.                                                          */
/* (for development only. there is no need to set these flags unless         */
/* changes are made in the library.)                                         */

/* use magic cookies to validate type of pointer */
#define UNUR_COOKIES  1


/*****************************************************************************
 *  Compile time parameters for generators.                                  *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Distribution objects.                                                     */

/* Maximal number of parameters for the PDF of a distribution.               */
/* (It must be at least 5!)                                                  */
#define UNUR_DISTR_MAXPARAMS  5

/* Maximal size of automatically created probability vectors.                */
#define UNUR_MAX_AUTO_PV    100000


/*****************************************************************************
 *  Interface for uniform random number generators.                          *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Default generators                                                        */

/* Macros 'UNUR_URNG_DEFAULT' and 'UNUR_URNG_AUX_DEFAULT' must expand to     */
/* (a function returning) an object of type UNUR_URNG.                       */
/*                                                                           */
/* IMPORTANT!                                                                */
/*                                                                           */
/*   When a particular default is chosen then the corresponding library must */
/*   be installed and the wrapper functions must be enabled by defining the  */
/*   above macros. Moreover, the library must linked when creating an        */
/*   executable.                                                             */
/*                                                                           */
/*   If you have none of these libraries use type FVOID.                     */
/*                                                                           */

/* The following is switched on by the configure flag                        */
/*   --with-urng-default=rngstream                                           */

#ifdef UNUR_URNG_DEFAULT_RNGSTREAM
/* use type RNGSTREAM (recommended!) */
#  define UNUR_URNG_DEFAULT      (unur_urng_rngstream_new("URNG_main"))
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_rngstream_new("URNG_aux"))

/* else we assume that we have                                               */
/*   --with-urng-default=builtin                                             */

#else
/* use type FVOID (built-in) */
#  define UNUR_URNG_DEFAULT      (unur_urng_builtin())
#  define UNUR_URNG_AUX_DEFAULT  (unur_urng_builtin_aux())
#endif

/* However, it is also possible to use any other URNG.                       */
/* Then one must comment out the above #definitions (and thus disables the   */
/* the corresponding configure flags) and #defines new ones.                 */
/* Here are two examples:                                                    */

/* use type PRNG */
/* #define UNUR_URNG_DEFAULT      (unur_urng_prng_new("mt19937(19863)")) */
/* #define UNUR_URNG_AUX_DEFAULT  (unur_urng_prng_new("LCG(2147483647,16807,0,1)")) */

/* use type GSL */
/* #define UNUR_URNG_DEFAULT      (unur_urng_gsl_new(gsl_rng_mt19937)) */
/* #define UNUR_URNG_AUX_DEFAULT  (unur_urng_gsl_new(gsl_rng_cmrg)) */

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/
