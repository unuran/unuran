/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cemp.c                                                       *
 *                                                                           *
 *   manipulate empirical univariate continuous distribution objects         *
 *   (i.e. samples)                                                          *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "cemp.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cemp

/*---------------------------------------------------------------------------*/

static void _unur_distr_cemp_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** empirical univariate continuous distributions                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cemp_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: empirical univariate continuous                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CEMP);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CEMP;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* name for distribution */
  distr->name = "(empirical)";
  distr->name_str = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cemp_free;

  /* clone */
  distr->clone = _unur_distr_cemp_clone;

  /* set defaults */

  /* observed sample */
  DISTR.sample    = NULL;   /* sample      */
  DISTR.n_sample  = 0;      /* sample size */

  /* histogram */
  DISTR.n_hist    = 0;          /* number of bins          */
  DISTR.hist_prob = NULL;       /* probabilities for bins  */
  DISTR.hmin      = -INFINITY;  /* lower ...               */
  DISTR.hmax      = INFINITY;   /* ... and upper bound     */       
  DISTR.hist_bins = NULL;       /* boundary between bins   */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_cemp_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cemp_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.cemp

  struct unur_distr *clone;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CEMP, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy data about sample into generator object (when there is one) */
  if (DISTR.sample) {
    CLONE.sample = _unur_xmalloc( DISTR.n_sample * sizeof(double) );
    memcpy( CLONE.sample, DISTR.sample, DISTR.n_sample * sizeof(double) );
  }

  /* copy histogram into generator object (when there is one) */
  if (DISTR.hist_prob) {
    CLONE.hist_prob = _unur_xmalloc( DISTR.n_hist * sizeof(double) );
    memcpy( CLONE.hist_prob, DISTR.hist_prob, DISTR.n_hist * sizeof(double) );
  }

  /* copy user name for distribution */
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_cemp_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cemp_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;

  COOKIE_CHECK(distr,CK_DISTR_CEMP,RETURN_VOID);

  if (DISTR.sample)    free( DISTR.sample );
  if (DISTR.hist_prob) free( DISTR.hist_prob );
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_cemp_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cemp_set_data( struct unur_distr *distr, const double *sample, int n_sample )
     /*----------------------------------------------------------------------*/
     /* set observed sample for distribution                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   sample   ... pointer to array of observations                      */
     /*   n_sample ... number of observations (sample size)                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, sample, UNUR_ERR_NULL );

  /* check new parameter for generator */
  if (n_sample <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"sample size");
    return UNUR_ERR_DISTR_SET;
  }

  /* allocate memory for sample */
  DISTR.sample = _unur_xmalloc( n_sample * sizeof(double) );
  if (!DISTR.sample) return UNUR_ERR_MALLOC;

  /* copy observed sample */
  memcpy( DISTR.sample, sample, n_sample * sizeof(double) );
  DISTR.n_sample = n_sample;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cemp_set_data() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cemp_read_data( struct unur_distr *distr, const char *filename )
     /*----------------------------------------------------------------------*/
     /* Read data from file `filename'.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   filename ... name of data file                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );

  /* read data from file */
  DISTR.n_sample = _unur_read_data( filename, 1, &(DISTR.sample) );

  /* o.k. ? */
  return (DISTR.n_sample > 0) ? UNUR_SUCCESS : UNUR_ERR_DISTR_DATA;

} /* end of unur_distr_cemp_read_data() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cemp_get_data( const struct unur_distr *distr, const double **sample )
     /*----------------------------------------------------------------------*/
     /* get number of observations and set pointer to array observations     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   sample   ... pointer to observed sampletor                         */
     /*                                                                      */
     /* return:                                                              */
     /*   sample size                                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CEMP, 0 );

  *sample = (DISTR.sample) ? DISTR.sample : NULL;
  return DISTR.n_sample;

} /* end of unur_distr_cemp_get_data() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cemp_set_hist( struct unur_distr *distr, const double *prob, 
			  int n_hist, double xmin, double xmax )
     /*----------------------------------------------------------------------*/
     /* set histogram with bins of equal length for distribution             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   prob     ... pointer to array of distributions                     */
     /*   n_hist   ... number of bins                                        */
     /*   xmin     ... lower bound of histogram                              */
     /*   xmax     ... upper bound of histogram                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CEMP, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, prob, UNUR_ERR_NULL );

  /* check new parameter for generator */
  if (n_hist <= 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram size");
    return UNUR_ERR_DISTR_SET;
  }
  if (xmin >= xmax) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"histogram, min >= max");
    return UNUR_ERR_DISTR_SET;
  }

  /* allocate memory for hist */
  DISTR.hist_prob = _unur_xmalloc( n_hist * sizeof(double) );
  if (!DISTR.hist_prob) return UNUR_ERR_MALLOC;

  /* copy probabilities */
  memcpy( DISTR.hist_prob, prob, n_hist * sizeof(double) );
  DISTR.n_hist = n_hist;

  /* store boundaries of histogram */
  DISTR.hmin = xmin;
  DISTR.hmax = xmax;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cemp_set_hist() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cemp_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*   printvector ... print observed sample if not 0                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CEMP,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous univariate distribution (ie. a sample)\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  if (DISTR.n_sample>0) {
    /* observed samples */
    fprintf(log,"%s:\tsample size = %d",genid,DISTR.n_sample);
    if (printvector) {
      for (i=0; i<DISTR.n_sample; i++) {
	if (i%10 == 0)
	  fprintf(log,"\n%s:\t",genid);
	fprintf(log,"  %.5f",DISTR.sample[i]);
      }
    }
    fprintf(log,"\n%s:\n",genid);
  }

  if (DISTR.n_hist>0) {
    /* histogram */
    fprintf(log,"%s:\thistogram: #bins = %d, ",genid,DISTR.n_hist);
    fprintf(log,"min = %g, max = %g, length = %g\n",
	    DISTR.hmin, DISTR.hmax, (DISTR.hmax-DISTR.hmin)/DISTR.n_hist);
    if (printvector) {
      fprintf(log,"%s:\t> bin probabilities = ",genid);
      for (i=0; i<DISTR.n_hist; i++) {
	if (i%10 == 0)
	  fprintf(log,"\n%s:\t",genid);
	fprintf(log,"  %.5f",DISTR.hist_prob[i]);
      }
    }
    fprintf(log,"\n%s:\n",genid);
  }
} /* end of _unur_distr_cemp_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
