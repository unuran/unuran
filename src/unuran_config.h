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
#ifndef UNURAN_CONFIG_H_SEEN
#define UNURAN_CONFIG_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *  Logging and debugging.                                                   *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Default name of log file.                                                 */
/*                                                                           */
/* Use "stdout" to write all infos to stdout.                                */
/* If no log file should be used, #undef the macro UNUR_ENABLE_LOGGING below.*/

#define UNUR_LOG_FILE "unuran.log"
/*  #define UNUR_LOG_FILE "stdout" */

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

#define UNUR_DEBUGFLAG_DEFAULT   UNUR_DEBUG_INIT

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
/* defined, then there are no warnings and error messages at all.            */
/* However then it is recommend _not_ to define UNUR_WARNINGS_ON, either.    */

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
/* Invoke uniform random number generator (URNG).                            */
/*                                                                           */
/* There are several ways to use a uniform pseudo random number generator:   */
/*                                                                           */
/* UNUR_URNG_GENERIC:                                                        */
/*     This the most flexible interface to URNGs.                            */
/*     All required data for a URNG are stored in an object.                 */
/*     Thus it is possible to use all the other generator types listed       */
/*     below without recompilation.                                          */
/*     The URNG objects are created and handled with appropriate calls.      */
/*     Another advantage of this approach is that a unified interface        */
/*     to different sources of random numbers is available.                  */
/*     In particular, one can use all sources of uniform random numbers      */
/*     for which the below compiler switches originally have been designed.  */
/*                                                                           */
/*...........................................................................*/
/*                                                                           */
/* The following compiler switches are for backward compatibility.           */
/* Their usage, however, is DEPRECIATED!                                     */
/*                                                                           */
/* These sources are still available under the new extended                  */
/* UNUR_URNG_GENERIC interface.                                              */
/*                                                                           */
/* UNUR_URNG_FVOID:                                                          */
/*     Use a pointer to the routine without an argment, i.e.                 */
/*        double uniform(void);                                              */
/*     E.g., the uniform generator included in UNURAN.                       */ 
/*                                                                           */
/* UNUR_URNG_PRNG:                                                           */
/*     Use a pointer to a uniform RNG object from the `prng' library.        */
/*     (see http://random.mat.sbg.ac.at/ftp/pub/software/gen/ or             */
/*     or http://statistik.wu-wien.ac.at/prng/).                             */
/*                                                                           */
/* UNUR_URNG_RNGSTREAMS:                                                     */
/*     Use a pointer to a uniform RNG object from Pierre L'Ecuyer's          */
/*     `RngStreams' library for multiple independent streams of              */
/*     pseudo-random numbers.                                                */
/*     (see http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/)         */
/*     A GNU-style package is available from                                 */
/*     http://statistik.wu-wien.ac.at/software/RngStreams/                   */
/*                                                                           */
/* UNUR_URNG_GSL:                                                            */
/*     Use a pointer to a uniform RNG object from the GNU Scientific         */
/*     Library (GSL).                                                        */
/*     (see http://www.gnu.org/software/gsl/)                                */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* IMPORTANT!                                                                */
/*   The generic interface has been changed and extended!                    */
/*   It is now the recommended interface to uniform random number generators */
/*   (URNG). For this reason there are calls to create and handle objects    */
/*   that contain all necessary data about a URNG (see ... ).                */
/*                                                                           */
/*   The structure cannot be used directly any more.                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Define the possible sources for uniform (pseudo-) random numbers.         */
#define UNUR_URNG_GENERIC    99   /* use a generic interface (recommended)   */
/*...........................................................................*/
#define UNUR_URNG_FVOID       1   /* use type `double urng(void)'            */
#define UNUR_URNG_PRNG        2   /* use prng-3.x                            */
#define UNUR_URNG_RNGSTREAMS  3   /* use RngStreams                          */
#define UNUR_URNG_GSL         4   /* use GNU Scientific Library              */
/*---------------------------------------------------------------------------*/
/* Set type of the API to uniform random number generators.                  */
/* (Define one of the following)                                             */

/* Recommended type */
#define UNUR_URNG_TYPE UNUR_URNG_GENERIC

/* For backward compatibility only! */
/* (These random number generators can be used more flexible via the         */
/*  generic interface!)                                                      */
/* #define UNUR_URNG_TYPE UNUR_URNG_FVOID */
/* #define UNUR_URNG_TYPE UNUR_URNG_PRNG */
/* #define UNUR_URNG_TYPE UNUR_URNG_RNGSTREAMS */
/* #define UNUR_URNG_TYPE UNUR_URNG_GSL */

/*---------------------------------------------------------------------------*/
/* Default generators                                                        */

/* Set a default generator type for the UNUR_URNG_GENERIC interface.         */
/* (Define one of the following)                                             */

/* IMPORTANT!                                                                */
/*   The corresponding library must be installed and linked when a           */
/*   particular default is chosen (with the exception of type                */
/*   UNUR_URNG_FVOID which is included in UNURAN).                           */
/*   If you have none of these libraries use UNUR_URNG_FVOID.                */

#define UNUR_URNG_DEFAULT_TYPE UNUR_URNG_FVOID
/* #define UNUR_URNG_DEFAULT_TYPE UNUR_URNG_PRNG */
/* #define UNUR_URNG_DEFAULT_TYPE UNUR_URNG_RNGSTREAMS */
/* #define UNUR_URNG_DEFAULT_TYPE UNUR_URNG_GSL */

/*...........................................................................*/
/* Default generators for different sources of uniform random numbers.       */
/* (These defaults are only used if the source is chosen via                 */
/*  UNUR_URNG_TYPE / UNUR_URNG_DEFAULT_TYPE.)                                */ 

/* type: FVOID */
#define UNUR_URNG_DEFAULT_FVOID           (unur_urng_MRG31k3p)
#define UNUR_URNG_DEFAULT_RESET_FVOID     (unur_urng_MRG31k3p_reset)
#define UNUR_URNG_AUX_DEFAULT_FVOID       (unur_urng_fish)
#define UNUR_URNG_AUX_DEFAULT_RESET_FVOID (unur_urng_fish_reset)

/* type: PRNG */
#define UNUR_URNG_DEFAULT_PRNG            ("mt19937(19863)")
#define UNUR_URNG_AUX_DEFAULT_PRNG        ("LCG(2147483647,16807,0,1)")

/* type: RNGSTREAMS */
/*   none required  */

/* type: GSL */
#define UNUR_URNG_DEFAULT_GSL             (gsl_rng_mt19937)
#define UNUR_URNG_AUX_DEFAULT_GSL         (gsl_rng_cmrg)

/*---------------------------------------------------------------------------*/
/* Enable interfaces to different sources.                                   */
/* Uncomment the following lines when the source of uniform random numbers   */
/* should be used. However, the corresponding libraries have to be linked    */
/* into each executable.                                                     */


/* Use Otmar Lendl's `prng' library                                          */
/*    http://statistik.wu-wien.ac.at/prng/                                   */
/* #define UNURAN_HAS_PRNG 1 */

/* Use Pierre L'Ecuyer's `RngStreams' library for multiple independent       */
/*    streams of pseudo-random numbers                                       */
/*    (see http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/).         */
/* UNURAN makes use of a GNU-style package which is available from           */
/*    http://statistik.wu-wien.ac.at/software/RngStreams/                    */
/* #define UNURAN_HAS_RNGSTREAMS 1 */

/* Use uniform RNG objects from the GNU Scientific Library (GSL)             */
/*    http://www.gnu.org/software/gsl/                                       */
/* #define UNURAN_HAS_GSL 1 */

/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/
