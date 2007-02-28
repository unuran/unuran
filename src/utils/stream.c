/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: stream.c                                                          *
 *                                                                           *
 *   routines for output streams and reading data                            *
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
#include <ctype.h>

/*---------------------------------------------------------------------------*/

static FILE *_unur_logfile_open( void );  

/*---------------------------------------------------------------------------*/

static FILE *unur_stream = NULL;
static const char GENID_UNKNOWN[] = "UNURAN";

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

FILE * 
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
    unur_stream = _unur_logfile_open();
  }

  return unur_stream;
} /* end of unur_get_stream() */

/*---------------------------------------------------------------------------*/

static FILE *
_unur_logfile_open( void )
     /*----------------------------------------------------------------------*/
     /* open log file                                                        */
     /*----------------------------------------------------------------------*/
{
  static FILE* LOG = NULL;
  time_t started;   

  if (LOG) return LOG;  /* log file already open */

#ifdef UNUR_LOG_FILE 
  /* open log file */
  LOG = fopen(UNUR_LOG_FILE,"w");
#else
  /* use stdout instead of a log file */
  LOG = stdout;
#endif

  /* write header into log file */
  fprintf(LOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");
  fprintf(LOG,"Version: %s\n",PACKAGE_VERSION);

  /* time when created */
  if (time( &started ) != -1)
    fprintf(LOG,"%s",ctime(&started));

  fprintf(LOG,"\n====================================================\n\n");

  /* return file handler */
  return LOG;

} /* end of _unur_logfile_open() */

/*---------------------------------------------------------------------------*/

int
_unur_read_data( const char *filename, int no_of_entries, double **ar )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   filename      ... name of data file                                */
     /*   no_of_entries ... number of entries per line                       */                
     /*   ar            ... to store pointer to array                        */
     /*                                                                      */
     /* return:                                                              */
     /*   number of valid lines read.                                        */
     /*   the pointer to the double array is stored in ar.                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0, ar is set to NULL.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This function takes the name of a data file and read the first     */
     /*   `no_of_entries' double numbers of each line successive into an     */
     /*   array. All lines not starting with a number are skipped            */
     /*   (including lines starting with white space!).                      */
     /*   If a line contains less than `no_of_entries' numbers, the          */
     /*   function reports an error, `ar' is set to NULL, and 0 is returned. */
     /*   The function allocates the required double array as a side effect. */
     /*   In case of an error this array is freed.                           */
     /*----------------------------------------------------------------------*/
{

#define LINELENGTH  1024      /* max length of lines allowed    */

  /* ------------------------------------------------------- */
  /* variable declarations                                   */

  const int datasize = 1000; /* initial size of data array   */
  int i, j;
  int memfactor = 1;

  char line[LINELENGTH];
  char *toline;
  char *chktoline;

  double *data;              /* pointer to data array        */
  int n_data;                /* number of data in array      */

  FILE *fp;

  /* ------------------------------------------------------- */

  /* initialize array ar */
  *ar = NULL;
  n_data = 0;

  /* array must be able to hold at least no_of_entries numbers */
  if (datasize < no_of_entries) {
    _unur_error("read_data",UNUR_ERR_GEN_DATA,"No of entries > max datasize");
    return 0;   
  }

  /* allocate memory for data */
  data = _unur_xmalloc(memfactor * datasize * sizeof(double));

  /* open the file with the data */
  fp = fopen(filename, "r");
  if (fp == NULL) {
    _unur_error("read_data",UNUR_ERR_GENERIC,"cannot open file");
    free(data);
    return 0; 
  }

  /* read lines until eof */
  for ( fgets(line, LINELENGTH, fp), i=0;
        !feof(fp);
        fgets(line, LINELENGTH, fp) ) {

    /* if necessary allocate more memory for array */
    if (i > memfactor*datasize - no_of_entries-2){
      memfactor++;
      data = _unur_xrealloc(data, memfactor*datasize*sizeof(double));
    }

    /* ignore all lines not starting with a double number */
    if ( ! (isdigit(line[0]) || line[0] == '.' || line[0] == '+' 
           || line[0] == '-' ) )
      continue;

    /* increase counter */
    ++n_data;

    /* read data from line */    
    toline = line;  /* pointer to first element of array */    
    for (j=0 ; j<no_of_entries; i++, j++){
      chktoline = toline;
      data[i] = strtod(toline, &toline); /* success -> toline changes */

      /* no success reading a double */
      if (chktoline == toline) {
	_unur_error("read_data",UNUR_ERR_GEN_DATA,"data file not valid");
	free(data);
	return 0;    /* terminate routine */
      }  

    } /* end of for -- read data from line */
  }   /* end of for -- read lines of file  */ 

  /* allocate exactly the memory needed */
  data = _unur_xrealloc( data, (i+1) * sizeof(double) );

  /* o.k. */
  *ar = data;
  return n_data;

  /* ------------------------------------------------------- */

#undef LINELENGTH

} /* end of _unur_read_data() */

/*---------------------------------------------------------------------------*/
