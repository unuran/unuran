/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_misc.c                                                          *
 *                                                                           *
 *   miscellaneous routines                                                  *
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
#include <ctype.h>
/*---------------------------------------------------------------------------*/

int _unur_read_data( const char *file, int no_of_entries, double **array );

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
  data = _unur_malloc(memfactor * datasize * sizeof(double));

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
      data = _unur_realloc(data, memfactor*datasize*sizeof(double));
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
  data = _unur_realloc( data, (i+1) * sizeof(double) );

  /* o.k. */
  *ar = data;
  return n_data;

  /* ------------------------------------------------------- */

#undef LINELENGTH

} /* end of _unur_read_data() */

/*---------------------------------------------------------------------------*/



