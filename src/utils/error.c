/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: error.c                                                           *
 *                                                                           *
 *   routines for warnings and error messages                                *
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

#include <unur_source.h>

#include <time.h>
#include <stdarg.h>

/*---------------------------------------------------------------------------*/
/* prototypes for unuran error handlers                                      */

void _unur_error_handler_default( const char *objid, const char *file, int line, 
				  const char *errortype, int errorcode, const char *reason );

void _unur_error_handler_off( const char *objid, const char *file, int line, 
			      const char *errortype, int errorcode, const char *reason );

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

/* global variable used to record errors   */
int unur_errno = UNUR_SUCCESS;

/* error handler used to report errors in UNURAN */
static unur_error_handler_t *_unur_error_handler = _unur_error_handler_default;

/*---------------------------------------------------------------------------*/

void 
_unur_error_x( const char *objid, const char *file, int line, 
	       const char *errortype, int errorcode, const char *reason )
     /*----------------------------------------------------------------------*/
     /* error handler                                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... file name (inserted by __FILE__)                     */
     /*   line      ... line number in source file (inserted by __LINE__)    */ 
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNURAN error code                                    */
     /*   reason    ... (very) short description of reason for error         */
     /*----------------------------------------------------------------------*/
{
  unur_errno = errorcode;
  _unur_error_handler(objid, file, line, errortype, errorcode, reason);
} /* end of _unur_error_x() */

/*---------------------------------------------------------------------------*/

void
_unur_error_handler_default( const char *objid, const char *file, int line, 
			     const char *errortype, int errorcode, const char *reason )
     /*----------------------------------------------------------------------*/
     /* default error handler                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... file name (inserted by __FILE__)                     */
     /*   line      ... line number in source file (inserted by __LINE__)    */ 
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNURAN error code                                    */
     /*   reason    ... (very) short description of reason for error         */
     /*----------------------------------------------------------------------*/
{
  FILE *log = unur_get_stream();

  /* generator identifier known ? */
  if (!objid) objid = "UNURAN";

  fprintf(log,"%s: [%s] %s:%d - %s:\n", objid, errortype, file, line,
	  unur_get_strerror(errorcode));
  if (reason && strlen(reason))
    fprintf(log,"%s: ..>  %s\n", objid, reason);
  fflush(log);   /* in case of a segmentation fault */

} /* end of _unur_error_handler_default() */

/*---------------------------------------------------------------------------*/

void
_unur_error_handler_off( const char *objid, const char *file, int line, 
			 const char *errortype, int errorcode, const char *reason )
     /*----------------------------------------------------------------------*/
     /* disable error handler (except for logging)                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   objid     ... id/type of object                                    */
     /*   file      ... file name (inserted by __FILE__)                     */
     /*   line      ... line number in source file (inserted by __LINE__)    */ 
     /*   errortype ... "warning" or "error"                                 */
     /*   errorcode ... UNURAN error code                                    */
     /*   reason    ... (very) short description of reason for error         */
     /*----------------------------------------------------------------------*/
{
  return;
} /* end of _unur_error_handler_off() */

/*---------------------------------------------------------------------------*/
#ifdef UNUR_COOKIES

void
_unur_error_cookies( const char *file, int line, unsigned observed, unsigned expected )
     /*----------------------------------------------------------------------*/
     /* print error message: invalid cookie detected                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   file      ... file name (inserted by __FILE__)                     */
     /*   line      ... line number in source file (inserted by __LINE__)    */
     /*   observed  ... observed cookie                                      */
     /*   expected  ... expected cookie                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "observed = %#x, expected = %#x", observed, expected );
  _unur_error_x( "COOKIE", file, line, "error", UNUR_ERR_COOKIE, reason->text);
  _unur_string_free( reason );
} /* end of _unur_error_cookies() */

#endif
/*---------------------------------------------------------------------------*/

const char *
unur_get_strerror ( const int unur_errno )
     /*----------------------------------------------------------------------*/
     /* return string that describes error                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   unur_error ... error code                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to charater string                                         */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  switch (unur_errno) {

    /** procedure executed successfully **/
  case UNUR_SUCCESS:
    return "(no error)";

    /** distribution object **/
  case UNUR_ERR_DISTR_NPARAMS:
    return "(distribution) invalid number of parameters";
  case UNUR_ERR_DISTR_DOMAIN:
    return "(distribution) parameter out of domain";
  case UNUR_ERR_DISTR_GEN:
    return "(distribution) invalid variant for special generator";
  case UNUR_ERR_DISTR_INVALID:
    return "(distribution) invalid distribution object";
  case UNUR_ERR_DISTR_REQUIRED:
    return "(distribution) incomplete distribution object, entry missing";
  case UNUR_ERR_DISTR_UNKNOWN:
    return "(distribution) unknown distribution, cannot handle";
  case UNUR_ERR_DISTR_SET:
    return "(distribution) set failed (invalid parameter)";
  case UNUR_ERR_DISTR_GET:
    return "(distribution) get failed (parameter not set)";
  case UNUR_ERR_DISTR_DATA:
    return "(distribution) data are missing (cannot execute)";
  case UNUR_ERR_DISTR_PROP:
    return "(distribution) desired property does not exist";

    /** parameter object **/
  case UNUR_ERR_PAR_SET:
    return "(parameter) set failed, invalid parameter -> using default";
  case UNUR_ERR_PAR_VARIANT:
    return "(parameter) invalid variant -> using default";
  case UNUR_ERR_PAR_INVALID:
    return "(parameter) invalid parameter object";

    /** generator object **/
  case UNUR_ERR_GEN_DATA:
    return "(generator) (possible) invalid data";
  case UNUR_ERR_GEN_CONDITION:
    return "(generator) condition for method violated";
  case UNUR_ERR_GEN_INVALID:
    return "(generator) invalid generator object";
  case UNUR_ERR_GEN_SAMPLING:
    return "(generator) sampling error";
  case UNUR_ERR_NO_REINIT:
    return "(generator) reinit missing";
  case UNUR_ERR_GEN:
    return "(generator)";


    /** uniform random number generator (URNG) object **/
  case UNUR_ERR_URNG:
    return "(URNG)";
  case UNUR_ERR_URNG_MISS:
    return "(URNG) missing functionality";
    
    /** string parser **/
  case UNUR_ERR_STR:
    return "(parser) invalid string";
  case UNUR_ERR_STR_UNKNOWN:
    return "(parser) unknown keyword";
  case UNUR_ERR_STR_SYNTAX:
    return "(parser) syntax error";
  case UNUR_ERR_STR_INVALID:
    return "(parser) invalid parameter";
  case UNUR_ERR_FSTR_SYNTAX:
    return "(function parser) syntax error";
  case UNUR_ERR_FSTR_DERIV:
    return "(function parser) cannot derivate function";

    /** misc **/
  case UNUR_ERR_DOMAIN:
    return "argument out of domain";
  case UNUR_ERR_ROUNDOFF:
    return "(serious) round-off error";
   case UNUR_ERR_MALLOC:
    return "virtual memory exhausted";
  case UNUR_ERR_NULL:
    return "invalid NULL pointer";
  case UNUR_ERR_COOKIE:
    return "invalid cookie";
  case UNUR_ERR_SILENT:
    return "(silent error)";
  case UNUR_ERR_GENERIC:
    return "(generic error)";
  case UNUR_ERR_INF:
    return "invalid infinity occured";
  case UNUR_ERR_NAN:
    return "NaN occured";

    /** compilation switches **/
  case UNUR_ERR_COMPILE:
    return "not available, recompile library";

    /** this should not happen **/
  case UNUR_ERR_SHOULD_NOT_HAPPEN:
  default:
    return "error should not happen, report this!";

  }

} /* end of unur_get_strerror() */

/*---------------------------------------------------------------------------*/

