/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unuran_errno.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines error codes.                                              *
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
#ifndef __X_ERRNO_H_SEEN
#define __X_ERRNO_H_SEEN
/*---------------------------------------------------------------------------*/

#include <stdio.h>

/*---------------------------------------------------------------------------*/

/*
  This chapter describes the way that UNURAN routines report
  errors. 

  @menu
  * Error reporting::             
  * Output streams::              
  @end menu

  @node Error reporting
  @section Error reporting

  UNURAN routines report an error whenever they cannot perform the
  task requested of them. 
  For example, apply transformed density rejection to a distribution
  that violates the T-concavity condition, or trying to set a
  parameter that is out of range.
  It might also happen that the setup fails for transformed density
  rejection for a T-concave distribution with some extreme density
  function simply because of round-off errors that makes the
  generation of a hat function numerically impossible.
  Situations like this may happen when using black box algorithms and
  you should check the return values of all routines.

  All ..._set_...(), and ..._chg_...() calls return @code{0} if it was
  not possible to set or change the desired parameters, e.g. because
  the given values are out of range, or simple because you have
  changed the method but not the corresponding set call and thus an
  invalid parameter or generator object is used.

  All routines that returns a pointer to the requested object will
  return a NULL pointer in case of error.
  (Thus you should always check the pointer to avoid possible
  segmentation faults. Sampling routines usually to not check the
  given pointer to the generator object. However you can switch on
  checking for NULL pointer defining the compiler switch 
  @code{UNUR_ENABLE_CHECKNULL} in unuran_config.h to avoid nasty
  segmentation faults.)

  The library distinguishes between two major classes of error:

  1. (fatal) errors: the library was not able to construct the
  requested object. 

  2. warnings: some problems encounters while construction a generator
  object. The routine has tried to solve the problem but the resulting
  object might not be what you want. For example, chosing a special
  variant of a method does not work and the initialization routine
  might switch to another variant. Then the generator produces random
  variates of the requested distribution but correlation induction
  is not possible. However it also might happen that 
  changing the domain of a distribution has failed. Then the generator
  produced random variates with too large/too small range, i.e. their
  distribution is different.

  It is obvious from the example that this distinction between errors
  and warning is rather crude and sometimes arbitrary. 

  UNURAN routines use the global variable @code{unuran_errno} to
  report errors, completely analogously to C library's
  @code{errno}. (However this approach is not thread-safe. There can 
  be only one instance of a global variable per program. Different
  threads of execution may overwrite @code{unuran_errno}
  simultaneously). 
  Thus when an error occurs the caller of the routine can examine the
  error code in @code{unuran_errno} to get more details about the
  reason why a routine failed. You get a short
  description of the error by a unur_get_strerror() call.
  All The error code numbers have prefix @code{UNUR_ERR_} and expand
  to non-zero constant unsigned integer values. 
  Error codes are divided into six main groups:

  @cindex error codes
*/
/*---------------------------------------------------------------------------*/
/* List of error codes:                                                      */
/*---------------------------------------------------------------------------*/

enum { 

  /** distribution object **/
  /*
    @code{UNUR_ERR_DISTR_...}
    Errors that occurred while handling distribution objects.
  */
  UNUR_ERR_DISTR_SET      = 0x11u,    /* set failed (invalid parameter)      */
  UNUR_ERR_DISTR_GET      = 0x12u,    /* get failed (parameter not set)      */
  UNUR_ERR_DISTR_NPARAMS  = 0x13u,    /* invalid number of parameters        */
  UNUR_ERR_DISTR_DOMAIN   = 0x14u,    /* parameter out of domain             */
  UNUR_ERR_DISTR_GEN      = 0x15u,    /* invalid variant for special generator */
  UNUR_ERR_DISTR_REQUIRED = 0x16u,    /* incomplete distribution object, entry missing */
  UNUR_ERR_DISTR_UNKNOWN  = 0x17u,    /* unknown distribution, cannot handle */
  UNUR_ERR_DISTR_INVALID  = 0x18u,    /* invalid distribution object         */
  UNUR_ERR_DISTR_DATA     = 0x19u,    /* data are missing                    */

  /** parameter object **/
  /*
    @code{UNUR_ERR_PAR_...}
    Errors that occurred while handling parameter objects.
  */
  UNUR_ERR_PAR_SET        = 0x21u,    /* set failed (invalid parameter)      */
  UNUR_ERR_PAR_VARIANT    = 0x22u,    /* invalid variant -> using default    */
  UNUR_ERR_PAR_INVALID    = 0x23u,    /* invalid parameter object            */

  /** generator object **/
  /*
    @code{UNUR_ERR_GEN_...}
    Errors that occurred while handling generator objects.
  */
  UNUR_ERR_GEN            = 0x31u,    /* bit for generator object            */
  UNUR_ERR_GEN_DATA       = 0x32u,    /* (possible) invalid data             */
  UNUR_ERR_GEN_CONDITION  = 0x33u,    /* condition for method violated       */
  UNUR_ERR_GEN_INVALID    = 0x34u,    /* invalid generator object            */

  /** misc **/
  /*
    @code{UNUR_ERR_...}
    Other errors.
  */
  UNUR_ERR_ROUNDOFF       = 0x02u,    /* (serious) round-off error           */
  UNUR_ERR_MALLOC         = 0x03u,    /* virtual memory exhausted            */
  UNUR_ERR_NULL           = 0x04u,    /* invalid NULL pointer                */ 
  UNUR_ERR_COOKIE         = 0x05u,    /* invalid cookie                      */
  UNUR_ERR_GENERIC        = 0x06u,    /* generic error                       */

  /** compilation switches **/
  /*
    @code{UNUR_ERR_COMPILE}
    Requested routine requires different compilation switches.
    Recompilation of library necessary.
  */
  UNUR_ERR_COMPILE        = 0x0eu,    /* not available, recompile library    */

  /** this should not happen **/
  /*
    @code{UNUR_ERR_SHOULD_NOT_HAPPEN}
    Internal error. This should not happen. 
    Please make a bug report.
  */
  UNUR_ERR_SHOULD_NOT_HAPPEN = 0x0fu, /* error should not happen, report this! */

};

/*---------------------------------------------------------------------------*/

/*
  @node Error handlers
  @section Error handlers
  @cindex Error handlers

  In addition to reporting error via the @code{unuran_errno} mechanism
  the library also provides an (optional) error handler. The error
  handler is called by the library functions when they are about to
  report an error. Then a short error diagnostics is written via two
  output streams. Both can be switched on/off by compiler flag
  @code{UNUR_WARNINGS_ON} in unuran_config.h.
  
  The first stream is the stderr. It can be enabled by defining 
  the macro @code{UNUR_ENABLE_STDERR} in unuran_config.h.

  The second stream can be set abritrarily by the unur_set_stream()
  call. If no such stream is given by the user a default stream is
  used by the library: all warnings and error messages are written
  into the file unuran.log in the current working directory.
  The name of this file defined by the macro @code{UNUR_LOG_FILE} in
  unuran_config.h. If the stdout should be used, define this macro by 
  "stdout".

  This output stream is also used to log descriptions of build generator
  objects and for writing debugging information.
  If you want to use this output stream for your own programs use 
  unur_get_stream() to get its file handler.
  This stream is enabled by the compiler switch
  @code{UNUR_ENABLE_LOGFILE} in unuran_config.h. 

  All warnings and error messages as well all debugging information
  are written onto the same output stream.
  To destinguish between the messages for different generators define
  the macro @code{UNUR_ENABLE_GENID} in unuran_config.h. Then every
  generator object has a unique identifier that is used for every
  message.

*/

/*---------------------------------------------------------------------------*/
/* warnings and error messages for given error number                        */

const char *unur_get_strerror ( const int unur_errno );
/*
  Get a short description for error code value.
*/

/*---------------------------------------------------------------------------*/
/* manipulate output stream                                                  */

FILE *unur_set_stream( FILE *new_stream );
/*
  Set new file handle for output stream; the old file handle is
  returned. The NULL pointer is not allowed. (If you want to disable
  logging of debugging information use 
  unur_set_default_debug(UNUR_DEBUG_OFF) instead.

  The output stream is used to report errors and warning, and
  debugging information. It is also used to log descriptions of
  build generator objects (when this feature is switched on; see also ?).
*/

FILE *unur_get_stream( void );
/*
  Get the file handle for the current output stream.
*/


/*---------------------------------------------------------------------------*/
/* global variable used to record errors                                     */

extern unsigned unur_errno;
/*
  This global variable is used to report diagnostics of error.
*/

/*---------------------------------------------------------------------------*/
#endif  /* __X_ERRNO_H_SEEN */
/*---------------------------------------------------------------------------*/
