/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: histogram.c                                                       *
 *                                                                           *
 *   Routines fordrawing histograms in the log file                          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
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
				                                                                                    
/*--------------------------------------------------------------------------*/

#include <unur_source.h>
#include "histogram_source.h"

/*---------------------------------------------------------------------------*/


void 
_unur_hist ( double *v , int length, int number_of_bins, const char *info, const char *genid )
     /*----------------------------------------------------------------------*/
     /* A histogram of the data-vector v is drawn into the logfile         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   v     ... vector to be histogrammized                              */
     /*   len   ... length of vector v                                       */
     /*   number_of_bins ... between min(v) and max(v)                       */
     /*   info  ... additional info string                                   */
     /*   genid ... id string                                                */
     /*----------------------------------------------------------------------*/
{

  FILE *log;
  int i,j;
  
  double vmin, vmax; /* minimum and maximum of the data-vactor v */
  long *hist; /* histogram data */
  long histmax; /* maximum entry in histogram */
  int scale_factor = 52; /* character-length of histogram */
  
  /* validate input values */
  /* TODO: if (length<=0) ... */
  if (number_of_bins<=0) number_of_bins = length / 100;
  
  /* allocate memory for the histogram */
  hist = _unur_malloc( (number_of_bins+1) * sizeof(long));
  
  log = unur_get_stream();

  fprintf(log,"%s: %s\n", genid, info); 
  fprintf(log,"%s\n", genid); 

  if (v==NULL) {
    fprintf(log,"%s: NULL pointer\n", genid);
  }

  else {
    /* calculate minimum and maximum of the data array v */
    vmin=v[0]; vmax=v[0];
    for (i=1; i<length; i++) {
      if (v[i]<vmin) vmin=v[i];
      if (v[i]>vmax) vmax=v[i];
    }

    /* check data range */ 
    if ((vmax-vmin) < UNUR_EPSILON) {
      vmax=vmin+UNUR_EPSILON;
    }
    
    /* reset histogram */
    for (j=0; j<=number_of_bins; j++) {
      hist[j]=0; 
    } 
    
    /* fill histogramm */
    for (i=0; i<length; i++) {
      hist[ (int) (number_of_bins * (v[i]-vmin)/(vmax-vmin)) ] += 1; 
    }

    /* calculate hist_max */
    int binmax; /* in which bin do we have the (first) maximum ? */
    histmax=0;
    for (j=0; j<=number_of_bins; j++) {
      if (histmax<hist[j]) {
        histmax=hist[j];
        binmax=j;
      }	
    }

    fprintf(log,"%s Data range %f <= x <= %f in bins #0 -> #%d\n", genid, vmin, vmax, number_of_bins); 
    fprintf(log,"%s Histogram maximum in bin #%d i.e. in interval [%f, %f) \n", 
                 genid, binmax, vmin + ((vmax-vmin) * binmax)/number_of_bins,
		 vmin + ((vmax-vmin) * (binmax+1))/number_of_bins ); 
    fprintf(log,"%s\n", genid); 
    for (j=0; j<=number_of_bins; j++) {

      fprintf(log, "%s: %7ld ", genid, hist[j]); 

      for (i=0; i < (scale_factor * hist[j])/histmax; i++) {
        fprintf(log, "*"); 
      }
      fprintf(log, "\n");
    }
  }

  fprintf(log,"%s\n", genid); 

  if (hist) free(hist);
  
} /* end of _unur_hist() */  

/*---------------------------------------------------------------------------*/

