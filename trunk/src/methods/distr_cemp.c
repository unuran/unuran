/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr_cemp.c                                                 *
 *                                                                           *
 *   manipulate empirical univariate continuous distribution objects         *
 *   (i.e. samples)                                                          *
 *                                                                           *
 *   return:                                                                 *
 *     1 ... on success                                                      *
 *     0 ... on error                                                        *
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

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
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

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cemp_free;

  /* set defaults                                                            */

  /* observed sample */
  DISTR.sample    = NULL;          /* sample                                 */
  DISTR.n_sample  = 0;             /* sample size                            */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cemp_new() */

/*---------------------------------------------------------------------------*/

int
_unur_distr_cemp_copy( struct unur_distr *to, const struct unur_distr *from )
     /*----------------------------------------------------------------------*/
     /* copy distribution object 'from' into distribution object 'to'.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   to   ... pointer to target distribution object                     */
     /*   from ... pointer to source distribution object                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
#define FROM from->data.cemp
#define TO   to->data.cemp

  int len;

  /* check arguments */
  _unur_check_NULL( NULL,from,0 );
  _unur_check_distr_object( from, CEMP, 0 );

  /* copy distribution object into generator object */
  memcpy( to, from, sizeof( struct unur_distr ) );

  /* copy data about sample into generator object (when there is one) */
  if (FROM.sample) {
    TO.sample = _unur_malloc( FROM.n_sample * sizeof(double) );
    memcpy( TO.sample, FROM.sample, FROM.n_sample * sizeof(double) );
  }

  /* copy user name for distribution */
  if (from->name_str) {
    len = strlen(from->name_str) + 1;
    to->name_str = _unur_malloc(len);
    memcpy( to->name_str, from->name_str, len );
    to->name = to->name_str;
  }

  return 1;

#undef FROM
#undef TO
} /* end of _unur_distr_cemp_copy() */

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

  COOKIE_CHECK(distr,CK_DISTR_CEMP,/*void*/);

  if (DISTR.sample) free( DISTR.sample );
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_cemp_free() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cemp_clear( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* frees all memory blocks in distribution object                       */
     /* object.                                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *distr = &(gen->distr);

  /* check arguments */
  COOKIE_CHECK(distr,CK_DISTR_CEMP,/*void*/);

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

} /* end of unur_distr_cemp_clear() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CEMP, 0 );

  /* check new parameter for generator */
  if (n_sample < 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"sample size");
    return 0;
  }

  /* allocate memory for sample */
  DISTR.sample = _unur_malloc( n_sample * sizeof(double) );
  if (!DISTR.sample) return 0;

  /* copy observed sample */
  memcpy( DISTR.sample, sample, n_sample * sizeof(double) );
  DISTR.n_sample = n_sample;

  /* o.k. */
  return 1;

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CEMP, 0 );

  /* read data from file */
  DISTR.n_sample = _unur_read_data( filename, 1, &(DISTR.sample) );

  /* o.k. ? */
  return (DISTR.n_sample > 0) ? 1 : 0;

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

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cemp_debug( const struct unur_distr *distr, const char *genid, int printvector )
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
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_CEMP,/*void*/);

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
} /* end of _unur_distr_cemp_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
