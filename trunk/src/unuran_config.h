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
/* write info about generator and trace of generation into logfile.          */
#define UNUR_ENABLE_LOGGING

/*---------------------------------------------------------------------------*/
/* make a generator id (string stored in generator object)                   */
#define UNUR_ENABLE_GENID

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
/* Detailed discription of possible flags in method files.                   */
/*                                                                           */
#define UNUR_DEBUGFLAG_DEFAULT   UNUR_DEBUG_ALL


/*---------------------------------------------------------------------------*/
/* (Default) name of log file.                                               */
/* use "stdout" to write all infos to stdout.                                */

/*  #define UNUR_LOG_FILE "unuran.log" */
#define UNUR_LOG_FILE "stdout"


/*---------------------------------------------------------------------------*/
/* warnings and error messages.                                              */

/* enable warnings and error messages */
#define UNUR_WARNINGS_ON

/* write warnings and errors into logfile */
#define UNUR_ENABLE_LOGFILE

/* write warnings and errors to stderr */
/*  #define UNUR_ENABLE_STDERR */


/*---------------------------------------------------------------------------*/
/* check for invalide NULL pointer                                           */
/* (needs not be enabled unless the user is too lazy to verified pointers    */
/* to generator objects after initialization.)                               */
#define UNUR_ENABLE_CHECKNULL

/*---------------------------------------------------------------------------*/
/* Debugging tools                                                           */
/* (for development only. there is no need to set these flags unless         */
/* changes are made in the library.)                                         */

/* use magic cookies to validate type of pointer */
#define UNUR_COOKIES


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
/*---------------------------------------------------------------------------*/
/* define the possible compiler switches                                     */
#define UNUR_URNG_POINTER  2     /* use a pointer to a routine               */
#define UNUR_URNG_PRNG     3     /* use a pointer to gen. from prng-2.2      */
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
/* Precision for computations                                                */

/* epsilon */
#define UNUR_EPSILON  1e-14    /* should be about 100 * DBL_EPSILON          */

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/** TODO: we have to find a better place for these macros                    */
#define HAVE_GETTIMEOFDAY

#define UNUR_DISTR_MAXPARAMS  5 /* maximal number of parameters for distribution */

/* maximal size of automatically created probability vectors */
#define UNUR_MAX_AUTO_PV    100000

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __UNURAN_CONFIG_H_SEEN */
/*---------------------------------------------------------------------------*/


