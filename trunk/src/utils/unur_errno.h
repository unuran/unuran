/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_errno.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for error and message (to  *
 *         be written onto an outputstream) handling.                        *
 *                                                                           *
 *   USAGE:                                                                  *
 *         included in all source files that produce errors, warnings and    *
 *         other (debugging) messages.                                       *
 *         required for applications that need error handling.               *
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
#ifndef __UNUR_ERRNO_H_SEEN
#define __UNUR_ERRNO_H_SEEN
/*---------------------------------------------------------------------------*/

#include <errno.h>
#include <stdio.h>

#include <unur_defs.h>

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */

FILE *unur_set_stream( FILE *new_stream );
FILE *unur_get_stream( void );

void _unur_stream_printf( const char *genid, char *filename, int line, const char *format, ... );
char *_unur_make_genid( const char *gentype );


/*---------------------------------------------------------------------------*/
/* error types                                                               */

enum { 
  UNUR_SUCCESS           = 0, 

  UNUR_ERR_NULL          = 1,    /* invalid NULL pointer                      */
  UNUR_ERR_COOKIE        = 2,    /* invalid cookie                            */
  UNUR_ERR_ALLOC         = 3,    /* virtual memory exhausted                  */

  UNUR_ERR_NPARAM        = 11,   /* invalid number of parameters              */
  UNUR_ERR_PARAM         = 12,   /* invalid parameter                         */

  UNUR_ERR_SET           = 20,   /* SET failed (invalid parameter)            */
  UNUR_ERR_SET_INVALID   = 21,   /* SET failed (invalid parameter)            */
  UNUR_ERR_SET_NOTREQU   = 22,   /* SET failed (invalid parameter)            */
  UNUR_ERR_CHG           = 30,   /* CHG failed (invalid parameter)            */
  UNUR_ERR_CHG_INVALID   = 31,   /* CHG failed (invalid parameter)            */
  UNUR_ERR_CHG_NOTREQU   = 32,   /* CHG failed (invalid parameter)            */
  UNUR_ERR_GET           = 40,   /* GET failed (invalid parameter)            */
  UNUR_ERR_GET_INVALID   = 41,   /* GET failed (invalid parameter)            */
  UNUR_ERR_GET_NOTREQU   = 42,   /* GET failed (invalid parameter)            */

  UNUR_ERR_INIT          = 100,  /* INIT                                      */
  UNUR_ERR_INIT_FAILED   = 101,  /* INIT failed                               */
  UNUR_ERR_INIT_INVALID  = 102,  /* INIT failed (invalid parameter)           */
  UNUR_ERR_INIT_VIOLATE  = 103,  /* INIT failed (condition for method violated) */

  UNUR_ERR_SAMPLE        = 200,  /* SAMPLing error (condition for method violated) */
  UNUR_ERR_ADAPT         = 210,  /* ADAPTive step failed                      */
  UNUR_ERR_ADAPT_VIOLATE = 211,  /* ADAPTive step failed (condition for method violated) */

  UNUR_ERR_GENERIC       = 900,  /* generic error                             */
  UNUR_ERR_UNIMPLEMENTED = 999,  /* unimplemented feature                     */
  UNUR_ERR_UNKNOWN       = 1000, /* unknown error (report this!)              */

  UNUR_ERR_UNKNOWNDISTR  = 2000, /* cannot handle this distribution */

  /** TODO ??? **/
  UNUR_ERR_DISTR         = 1111,  /* invalid parameter for distribution        */


  /************************************************************/
  /** TODO **/


  UNUR_ERR_DISTR_INVALID = 500,   /* invalid distribution type */
  UNUR_ERR_PAR_INVALID   = 501,   /* invalid parameter object */
  UNUR_ERR_GEN_INVALID   = 502,   /* invalid generator object */

  UNUR_ERR_DISTR_REQUIRED = 503,   /* entry in distribution object required */

};


/*---------------------------------------------------------------------------*/
#if UNUR_DEBUG > 0
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/* warnings and error messages                                               */

const char *unur_get_strerror ( const int unur_errno );
/* return string that describes error                                        */

#define _unur_error(genid,errortype,str) \
   do { \
      _unur_stream_printf((genid),__FILE__,__LINE__,"error: %s %s", \
                          unur_get_strerror(errortype), (str) ); \
   } while (0)

#define _unur_warning(genid,errortype,str) \
   do { \
      _unur_stream_printf((genid),__FILE__,__LINE__,"warning: %s %s", \
                          unur_get_strerror(errortype), (str) ); \
   } while (0)

/*---------------------------------------------------------------------------*/
/* an identifier for the generator object                                    */

/* set generator id */
#define _unur_set_genid(par,gentype)  (par)->genid = _unur_make_genid(gentype)
#define _unur_free_genid(gen)         free((gen)->genid)

/*---------------------------------------------------------------------------*/
/* write infos into log file                                                 */

#if UNUR_DEBUG & UNUR_DB_INFO

/* debugging flags for generators */
#define _unur_set_debugflag_default(par)       (par)->debug = UNUR_DEBUGFLAG_DEFAULT
#define _unur_copy_debugflag(par,gen)          (gen)->debug = (par)->debug

#define _unur_print_if_default(par,flag)   if(!((par)->set & (flag))) fprintf(log,"  [default]")

#else

#define _unur_set_debugflag_default(par)
#define _unur_copy_debugflag(par,gen)

#define _unur_print_if_default(par,flag)

#endif   /* UNUR_DB_INFO */

/*---------------------------------------------------------------------------*/
#else    /* no debugging */
/*---------------------------------------------------------------------------*/

#define _unur_error(genid,errortype,str)      do { } while(0)
#define _unur_warning(genid,errortype,str)    do { } while(0)

#define _unur_set_genid(gen,gentype)
#define _unur_free_genid(gen)

#define _unur_set_debugflag_default(par)
#define _unur_copy_debugflag(par,gen)

#define _unur_print_if_default(par,flag)

/*---------------------------------------------------------------------------*/
#endif   /* UNUR_DEBUG */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_ERRNO_H_SEEN */
/*---------------------------------------------------------------------------*/






