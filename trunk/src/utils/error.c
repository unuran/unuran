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
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Thu Jul 22 16:07:30 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <time.h>
#include <stdarg.h>

#include <unur_errno.h>
#include <unur_umalloc.h>

/*---------------------------------------------------------------------------*/

static FILE *_unur_logfile_open( const char *filename );  

/*---------------------------------------------------------------------------*/

static FILE *unur_stream = NULL;

/*---------------------------------------------------------------------------*/

FILE* 
unur_set_log( FILE *new_stream )
     /*----------------------------------------------------------------------*/
     /* (re)set output stream for (error) messages                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   new_stream ... pointer to new output stream                        */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old stream                                              */
     /*----------------------------------------------------------------------*/
{
  FILE * previous_stream;
  
  if (unur_stream == NULL) {
    unur_stream = _unur_logfile_open(UNUR_LOG_FILE);
  }
  
  previous_stream = unur_stream;
  unur_stream = new_stream;
  
  return previous_stream;
} /* end of unur_set_log() */

/*---------------------------------------------------------------------------*/

FILE* 
unur_get_log( void )
     /*----------------------------------------------------------------------*/
     /* get output stream for (error) messages                               */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to output stream                                           */
     /*----------------------------------------------------------------------*/
{
  if (unur_stream == NULL) {
    unur_stream = _unur_logfile_open(UNUR_LOG_FILE);
  }

  return unur_stream;
} /* end of unur_get_log() */

/*---------------------------------------------------------------------------*/

static FILE*
_unur_logfile_open( const char *logfilename )
     /*----------------------------------------------------------------------*/
     /* open log file                                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   logfilename ... name of log file                                   */
     /*----------------------------------------------------------------------*/
{
  static FILE* LOG = NULL;
  time_t started;   

  if (LOG) return LOG;  /* log file already open */

  /* open log file */
  if (strcmp("stdout",logfilename))
    LOG = fopen(logfilename,"w");
  else /* use stdout instead of a log file */
    LOG = stdout;

#if UNUR_DEBUG & UNUR_DB_STDERR
  if (!LOG) fprintf(stderr,"warning: cannot open logfile %s\n",logfilename);
  fflush(stderr);   /* in case of a segmentation fault */
#endif

  /* write header into log file */
  fprintf(LOG,"\nUNURAN Universal Non-Uniform RANdom number generator\n\n");

  /* time when created */
  if (time( &started ) != -1)
    fprintf(LOG,"%s",ctime(&started));

  fprintf(LOG,"\n====================================================\n\n");

  /* return file handler */
  return LOG;

} /* end of _unur_open_logfile() */

/*---------------------------------------------------------------------------*/

char* 
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

/*---------------------------------------------------------------------------*/

static char GENID_UNKNOWN[] = "UNURAN";

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

void 
_unur_db_error( const char *genid, int errortype, char *filename, int line, const char *msg, ...)
     /*----------------------------------------------------------------------*/
     /* write error message                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genid     ... identifier of generator object (NULL if not known)   */
     /*   errortype ... type of error                                        */
     /*   filename  ... name of source file (provided by __FILE__)           */
     /*   line      ... line number in source file (provided by __LINE__)    */
     /*   msg       ... additional error message                             */
     /*   ...       ... optional arguments                                   */
     /*----------------------------------------------------------------------*/
{
  va_list ap;
  const char *errormsg;

  /* generator identifier known ? */
  if (!genid) genid = GENID_UNKNOWN;

  /* optional argmuents */
  va_start(ap, msg);

  /* get main warning message */
  errormsg = unur_get_strerror( errortype );

#if UNUR_DEBUG & UNUR_DB_STDERR   /* write warnings and errors on stderr */
  fprintf(stderr,"%s: error in %s, line %d: %s ",genid,filename,line,errormsg);
  vfprintf(stderr,msg,ap);
  fprintf(stderr,"\n");
  fflush(stderr);   /* in case of a segmentation fault */
#endif

#if UNUR_DEBUG & UNUR_DB_LOG      /* write warnings and errors into log file */
  if (!unur_stream) _unur_logfile_open(UNUR_LOG_FILE);
  fprintf(unur_stream,"%s: error in %s, line %d: %s ",genid,filename,line,errormsg);
  vfprintf(unur_stream,msg,ap);
  fprintf(unur_stream,"\n");
  fflush(unur_stream);   /* in case of a segmentation fault */
#endif

  /* terminate list of optional parameters */
  va_end(ap);

} /* end of _unur_db_error() */

/*---------------------------------------------------------------------------*/

void 
_unur_db_warning( const char *genid, int errortype, char *filename, int line, const char *msg, ...)
     /*----------------------------------------------------------------------*/
     /* write warning                                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genid     ... identifier of generator object (NULL if not known)   */
     /*   errortype ... type of error                                        */
     /*   filename  ... name of source file (provided by __FILE__)           */
     /*   line      ... line number in source file (provided by __LINE__)    */
     /*   msg       ... format and additional message                        */
     /*   ...       ... optional arguments                                   */
     /*----------------------------------------------------------------------*/
{
  va_list ap;
  const char *errormsg;

  /* generator identifier known ? */
  if (!genid) genid = GENID_UNKNOWN;

  /* optional argmuents */
  va_start(ap, msg);

  /* get main warning message */
  errormsg = unur_get_strerror( errortype );

#if UNUR_DEBUG & UNUR_DB_STDERR   /* write warnings and errors on stderr */
  fprintf(stderr,"%s: warning in %s, line %d: %s ",genid,filename,line,errormsg);
  vfprintf(stderr,msg,ap);
  fprintf(stderr,"\n");
  fflush(stderr);   /* in case of a segmentation fault */
#endif

#if UNUR_DEBUG & UNUR_DB_LOG      /* write warnings and errors into log file */
  if (!unur_stream) _unur_logfile_open(UNUR_LOG_FILE);
  fprintf(unur_stream,"%s: warning in %s, line %d: %s ",genid,filename,line,errormsg);
  vfprintf(unur_stream,msg,ap);
  fprintf(unur_stream,"\n");
  fflush(unur_stream);   /* in case of a segmentation fault */
#endif

  /* terminate list of optional parameters */
  va_end(ap);

} /* end of _unur_db_warning() */

/*---------------------------------------------------------------------------*/

