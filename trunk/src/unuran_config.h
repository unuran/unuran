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

/*  #define UNUR_ENABLE_CHECKNULL 1 */


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
/* Invoke uniform random number generator.                                   */
/*                                                                           */
/* There are several ways to use a uniform pseudo random number generator:   */
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
/* UNUR_URNG_RNGSTREAM:                                                      */
/*     Use a pointer to a uniform RNG object from Pierre L'Ecuyer's          */
/*     `RngStream' library for multiple multiple independent streams of      */
/*     pseudo-random numbers.                                                */
/*     (see http://www.iro.umontreal.ca/~lecuyer/myftp/streams00/c/)         */
/*                                                                           */
/* UNUR_URNG_GSL:                                                            */
/*     Use a pointer to a uniform RNG object from the GNU Scientific         */
/*     Library (GSL).                                                        */
/*     (see http://www.gnu.org/software/gsl/)                                */
/*                                                                           */
/* UNUR_URNG_GENERIC:                                                        */
/*     This a generic interface with limited support.                        */
/*     It uses a structure to store both a function call of type             */
/*     `double urng(void*)' and a void pointer to the parameter list.        */
/*     Both must be set directly using the structure (there are currently    */
/*     no calls that support this uniform RNG type).                         */
/*     To avoid a segmentation fault the UNURAN build-in uniform RNGs are    */
/*     as defaults (although these do not need a parameter list).            */
/*     See below for more details.                                           */
/*                                                                           */

/*---------------------------------------------------------------------------*/
/* Define the possible compiler switches.                                    */
#define UNUR_URNG_FVOID      1   /* use type `double urng(void)'             */
#define UNUR_URNG_PRNG       2   /* use prng-3.x                             */
#define UNUR_URNG_RNGSTREAM  3   /* use RngStream                            */
#define UNUR_URNG_GSL        4   /* use GNU Scientific Library               */
#define UNUR_URNG_GENERIC    99  /* use a generic interface (see below)      */
/*---------------------------------------------------------------------------*/

/* Set type of uniform random number generator.                              */
/*  #define UNUR_URNG_TYPE UNUR_URNG_FVOID */
/*  #define UNUR_URNG_TYPE UNUR_URNG_PRNG */
#define UNUR_URNG_TYPE UNUR_URNG_RNGSTREAM
/*  #define UNUR_URNG_TYPE UNUR_URNG_GSL */
/*  #define UNUR_URNG_TYPE UNUR_URNG_GENERIC */

/*---------------------------------------------------------------------------*/
#if UNUR_URNG_TYPE == UNUR_URNG_FVOID
/*---------------------------------------------------------------------------*/

/* Type of uniform random number generator.         (Don't touch this line!) */
typedef double (UNUR_URNG)(void);

/* Default name of uniform random number generator.                          */
/* Valid name of a C routine of type `double urng(void)'.                    */
#define UNUR_URNG_DEFAULT     unur_urng_MRG31k3p
#define UNUR_URNG_AUX_DEFAULT unur_urng_fish

/* Function call to uniform RNG.                    (Don't touch this line!) */
#define _unur_call_urng(urng)        ((*(urng))())

/* Function call to reset uniform RNG.                (For test suite only.) */
/* The routine must reset the the default uniform RNG (UNUR_URNG_DEFAULT).   */
/* Comment out the macro definition if a reset routine is not available.     */
#define unur_urng_reset(urng)        (unur_urng_MRG31k3p_reset());

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_PRNG
/*---------------------------------------------------------------------------*/

/* Header file from prng library.                  (Don't remove this line!) */
/* Make sure that this header file and the prng library is installed.        */
#include <prng.h>

/* Type of uniform random number generator.         (Don't touch this line!) */
typedef struct prng UNUR_URNG;

/* Default uniform random number generators.                                 */
/* These functions must return a pointer to a valid generator object of      */
/* type `struct prng *'  (see prng-3.x manual).                              */
#define UNUR_URNG_DEFAULT     (prng_new("mt19937(19863)"))
#define UNUR_URNG_AUX_DEFAULT (prng_new("LCG(2147483647,16807,0,1)"))

/* Function call to uniform RNG.                    (Don't touch this line!) */
#define _unur_call_urng(urng)        (prng_get_next(urng))

/* Function call to reset uniform RNG.                 (For test suite only) */
/*                                                  (Don't touch this line!) */
#define unur_urng_reset(urng)        (prng_reset(urng))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_RNGSTREAM
/*---------------------------------------------------------------------------*/

/* Header file from RngStream library.             (Don't remove this line!) */
/* Make sure that this header file and the RngStream library is installed.   */
#include <RngStream.h>

/* Type of uniform random number generator.         (Don't touch this line!) */
typedef struct RngStream_InfoState UNUR_URNG;

/* Default uniform random number generators.                                 */
/* These functions must return a pointer to a valid generator object of type */
/* `struct RngStream_InfoState *'  (see RngStream manual and RngStream.h).   */
#define UNUR_URNG_DEFAULT     (RngStream_CreateStream("URNG_main"))
#define UNUR_URNG_AUX_DEFAULT (RngStream_CreateStream("URNG_aux"))

/* Function call to uniform RNG.                    (Don't touch this line!) */
#define _unur_call_urng(urng)        (RngStream_RandU01(urng))

/* Function call to reset uniform RNG.                 (For test suite only) */
/*                                                  (Don't touch this line!) */
#define unur_urng_reset(urng)        (RngStream_ResetStartStream(urng))

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GSL
/*---------------------------------------------------------------------------*/

/* Header file from GNU Scientific Library.        (Don't remove this line!) */
/* Make sure that this header file and the GSL is installed.                 */
#include <gsl/gsl_rng.h>

/* Type of uniform random number generator.         (Don't touch this line!) */
typedef gsl_rng UNUR_URNG;

/* Default uniform random number generators.                                 */
/* These functions must return a pointer to a valid generator object of type */
/* `struct RngStream_InfoState *'  (see RngStream manual and RngStream.h).   */
#define UNUR_URNG_DEFAULT     (gsl_rng_alloc(gsl_rng_mt19937))
#define UNUR_URNG_AUX_DEFAULT (gsl_rng_alloc(gsl_rng_cmrg))

/* Function call to uniform RNG.                    (Don't touch this line!) */
#define _unur_call_urng(urng)        (gsl_rng_uniform_pos(urng))

/* Function call to reset uniform RNG.                 (For test suite only) */
/* No such routine for gsl uniform RNGs.                                     */
/*  #define unur_urng_reset(urng)        (exit(EXIT_FAILSURE)) */

/*---------------------------------------------------------------------------*/
#elif UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

/* This a generic interface with limited support.                            */
/*                                                                           */
/* IMPORTANT!                                                                */
/*   All functions and parameters should be set at run time:                 */
/*                                                                           */
/*     1. Allocate variable of type `struct unur_urng_generic'               */
/*        (or of type `UNUR_URNG', which is the same):                       */
/*          UNUR_URNG *urng;                                                 */
/*          urng = malloc(sizeof(UNUR_URNG));                                */
/*     2. Set function of type `double (*rand)(void *)' for sampling from    */
/*        uniform RNG:                                                       */
/*          urng->getrand = my_uniform_rng;                                  */
/*     3. Set pointer to parameters of for this function                     */
/*        (or NULL if no parameters are required):                           */
/*          urng->params = my_parameters;                                    */
/*     4. Use this uniform RNG:                                              */
/*          unur_urng_set_default(urng);             (set default generator) */
/*          unur_urng_set_default_aux(urng);    (set default aux. generator) */
/*        Notice that this must be done before UNURAN generator object       */
/*        are created. Of course urng can also used just for a particular    */
/*        generator object. Use the following and similar calls:             */ 
/*          unur_set_urng(par,urng);                                         */
/*          unur_chg_urng(gen,urng);                                         */
/*                                                                           */
/* To avoid a possible segmentation fault when this URNG type is set,        */
/* we use the build-in uniform RNGs as defaults.                             */
/* Be carefull when you change any macro definition!!                        */
/* Indead, there is no need to change these definitions (except              */
/* (un)commenting the below macro for reseting the uniform RNG), since all   */
/* pointers are set at run time!                                             */

/* Structure for storing function pointer and parameterlist.                 */
/*                                                  (Don't touch this line!) */
struct unur_urng_generic {
  double (*getrand)(void *params);  /* function for generating uniform RNG */
  void *params;                     /* list of parameters                  */
};

/* Type of uniform random number generator.         (Don't touch this line!) */
typedef struct unur_urng_generic UNUR_URNG;

/* Default uniform random number generators.                                 */
/* (Very ugly hack! Don't do this yourself unless you understand the code!)  */
typedef double (_unur_urng_voidptr)(void*);      /* cast function pointer */

#define UNUR_URNG_DEFAULT \
  malloc(sizeof(UNUR_URNG)); \
  urng_default->getrand = (_unur_urng_voidptr*)(unur_urng_MRG31k3p); \
  urng_default->params = NULL

#define UNUR_URNG_AUX_DEFAULT \
  malloc(sizeof(UNUR_URNG)); \
  urng_aux_default->getrand = (_unur_urng_voidptr*)(unur_urng_fish); \
  urng_aux_default->params = NULL

/* Function call to uniform RNG.                    (Don't touch this line!) */
#define _unur_call_urng(urng)      ((urng)->getrand((urng)->params))

/* Function call to reset uniform RNG.                 (For test suite only) */
/* The routine must reset the the default uniform RNG (UNUR_URNG_DEFAULT).   */
/* Comment out the macro definition if a reset routine is not available.     */
/*                                      (But don't change macro definition!) */
#define unur_urng_reset(urng)      (unur_urng_MRG31k3p_reset());

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
#error UNUR_URNG_TYPE not valid !!
/*---------------------------------------------------------------------------*/
#endif  /* UNUR_URNG_TYPE */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#endif  /* UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/


