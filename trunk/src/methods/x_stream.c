/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_stream.x                                                        *
 *                                                                           *
 *   routines for output streams                                             *
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

#include <source_unuran.h>

#include <time.h>
#include <stdarg.h>

/*---------------------------------------------------------------------------*/

static FILE *_unur_logfile_open( const char *filename );  
inline FILE *unur_get_stream( void );

/*---------------------------------------------------------------------------*/

static FILE *unur_stream = NULL;
const static char GENID_UNKNOWN[] = "UNURAN";

/*---------------------------------------------------------------------------*/

void
_unur_stream_printf( const char *genid, 
		     const char *filename, int line,
		     const char *format, ... )
     /*----------------------------------------------------------------------*/
     /* write messages on output stream(s)                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genid     ... identifier of generator object (NULL if not known)   */
     /*   filename  ... name of source file (provided by __FILE__)           */
     /*   line      ... line number in source file (provided by __LINE__)    */
     /*   format    ... format for fprintf()                                 */
     /*   ...       ... (optional) arguments to be be printed                */
     /*----------------------------------------------------------------------*/
{
  va_list ap;

  /* generator identifier known ? */
  if (!genid) genid = GENID_UNKNOWN;

  va_start(ap, format);

#ifdef UNUR_ENABLE_STDERR
  /* write on stderr */
  fprintf(stderr,"%s: %s:%d: ",genid,filename,line);
  vfprintf(stderr,format,ap);
  fprintf(stderr,"\n");
  fflush(stderr);   /* in case of a segmentation fault */
#endif

#ifdef UNUR_ENABLE_LOGFILE
  /* write onto output stream */
  if (!unur_stream) unur_get_stream();
  fprintf(unur_stream,"%s: %s:%d: ",genid,filename,line);
  vfprintf(unur_stream,format,ap);
  fprintf(unur_stream,"\n");
  fflush(unur_stream);   /* in case of a segmentation fault */
#endif

  va_end(ap);

} /* end of unur_stream_printf() */

/*---------------------------------------------------------------------------*/

void
_unur_stream_printf_simple( const char *format, ... )
     /*----------------------------------------------------------------------*/
     /* write messages on output stream(s)                                   */
     /* (same as _unur_stream_printf() but without file and line number)     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   genid     ... identifier of generator object (NULL if not known)   */
     /*   format    ... format for fprintf()                                 */
     /*   ...       ... (optional) arguments to be be printed                */
     /*----------------------------------------------------------------------*/
{
  va_list ap;

  va_start(ap, format);

#ifdef UNUR_ENABLE_STDERR
  /* write on stderr */
  vfprintf(stderr,format,ap);
  fflush(stderr);   /* in case of a segmentation fault */
#endif

#ifdef UNUR_ENABLE_LOGFILE
  /* write onto output stream */
  if (!unur_stream) unur_get_stream();
  vfprintf(unur_stream,format,ap);
  fflush(unur_stream);   /* in case of a segmentation fault */
#endif

  va_end(ap);

} /* end of unur_stream_printf_simple() */

/*---------------------------------------------------------------------------*/

FILE * 
unur_set_stream( FILE *new_stream )
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

  _unur_check_NULL( GENID_UNKNOWN,new_stream,NULL );

  previous_stream = unur_stream;
  unur_stream = new_stream;
  
  return previous_stream;
} /* end of unur_set_stream() */

/*---------------------------------------------------------------------------*/

inline FILE * 
unur_get_stream( void )
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
} /* end of unur_get_stream() */

/*---------------------------------------------------------------------------*/

static FILE *
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

#ifdef UNUR_ENABLE_STDERR
  if (!LOG) fprintf(stderr,"warning: cannot open logfile %s\n",logfilename);
  fflush(stderr);   /* in case of a segmentation fault */
#endif

  /* write header into log file */
  fprintf(LOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");

  /* time when created */
  if (time( &started ) != -1)
    fprintf(LOG,"%s",ctime(&started));

  fprintf(LOG,"\n====================================================\n\n");

  /* return file handler */
  return LOG;

} /* end of _unur_logfile_open() */

/*---------------------------------------------------------------------------*/

