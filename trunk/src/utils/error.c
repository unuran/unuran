/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      debug.c                                                      *
 *                                                                           *
 *   routines for debugging                                                  *
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

#include <time.h>
#include <stdarg.h>

#include <unur_errno.h>
#include <unur_umalloc.h>

/*---------------------------------------------------------------------------*/

char * 
_unur_make_genid( const char *gentype )
     /*----------------------------------------------------------------------*/
     /* make a new generator identifier                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gentype ... type of generator                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer generator id (char string)                                 */
     /*----------------------------------------------------------------------*/
{
  static int count = 0;   /* counter for identifiers */
  char *genid;

  /* allocate memory for identifier */
  genid = _unur_malloc(sizeof(char)*(strlen(gentype) + 6));

  /* make new identifier */
  ++count; count %= 1000;      /* 1000 different generators should be enough */
  sprintf(genid,"%s.%03d",gentype,count);

  return genid;

} /* end of _unur_make_genid() */

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
  case UNUR_SUCCESS:
    return "success";
    case UNUR_ERR_NULL:
      return "invalid NULL pointer";
  case UNUR_ERR_COOKIE:
    return "invalid cookie";
  case UNUR_ERR_ALLOC:
    return "virtual memory exhausted";

  case UNUR_ERR_NPARAM:
    return "invalid number of parameters";
  case UNUR_ERR_PARAM:
    return "invalid parameter";
    
  case UNUR_ERR_SET:
    return "SET failed";
  case UNUR_ERR_SET_INVALID:
    return "SET failed (invalid parameter)";
  case UNUR_ERR_SET_NOTREQU:
    return "SET failed (parameter not required)";
  case UNUR_ERR_CHG:
    return "CHG failed";
  case UNUR_ERR_CHG_INVALID:
    return "CHG failed (invalid parameter)";
  case UNUR_ERR_CHG_NOTREQU:
    return "CHG failed (parameter not required)";
  case UNUR_ERR_GET:
    return "GET failed";
  case UNUR_ERR_GET_INVALID:
    return "GET failed (invalid parameter)";
  case UNUR_ERR_GET_NOTREQU:
    return "GET failed (parameter not required)";

  case UNUR_ERR_INIT:
    return "INIT.";
  case UNUR_ERR_INIT_FAILED:
    return "INIT failed";
  case UNUR_ERR_INIT_INVALID:
    return "INIT failed (invalid parameter)";
  case UNUR_ERR_INIT_VIOLATE:
    return "INIT failed (condition for method violated)";
    
  case UNUR_ERR_SAMPLE:
    return "SAMPLing error (condition for method violated)";
    
  case UNUR_ERR_ADAPT:
    return "ADAPTive step failed";
  case UNUR_ERR_ADAPT_VIOLATE:
    return "ADAPTive step failed (condition for method violated)";

  case UNUR_ERR_GENERIC:
    return "";
  case UNUR_ERR_UNIMPLEMENTED:
    return "unimplemented feature";
    
  case UNUR_ERR_UNKNOWN:
  default:
    return "unknown error (report this!)";
  }

} /* end if unur_get_strerror() */

/*---------------------------------------------------------------------------*/

