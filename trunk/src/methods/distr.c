/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr.c                                                      *
 *                                                                           *
 *   manipulate distribution objects                                         *
 *                                                                           *
 *   PARAMETER: struct unur_distr *                                          *
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
#include <unur_stddistr.h>

/*---------------------------------------------------------------------------*/

/* create empty distribution object for ...                                  */
inline struct unur_distr *_unur_distr_cont_new( void );  /* univ. continuous */
inline struct unur_distr *_unur_distr_discr_new( void ); /* univ. discrete   */
inline struct unur_distr *_unur_distr_demp_new( void );  /* emp. univ. discr.*/
inline struct unur_distr *_unur_distr_cemp_new( void );  /* emp. univ. cont. */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** routines for all distribution objects                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_new( unsigned int type )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   type ... type of distribution                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  switch (type) {
  case UNUR_DISTR_CONT:
    return unur_distr_cont_new();
  case UNUR_DISTR_DISCR:
    return unur_distr_discr_new();
  case UNUR_DISTR_DEMP:
    return unur_distr_demp_new();
  case UNUR_DISTR_CEMP:
    return unur_distr_cemp_new();
  default:
    _unur_error(NULL,UNUR_ERR_DISTR_UNKNOWN,"");
    return NULL;
  }

} /* end of unur_distr_new() */

/*---------------------------------------------------------------------------*/

void
unur_distr_free( struct unur_distr *distr )
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

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,/*void*/);
    break;
  case UNUR_DISTR_DISCR:
    COOKIE_CHECK(distr,CK_DISTR_DISCR,/*void*/);
    break;
  case UNUR_DISTR_DEMP:
    COOKIE_CHECK(distr,CK_DISTR_DEMP,/*void*/);
    if (distr->data.demp.prob) free( distr->data.demp.prob );
    break;
  case UNUR_DISTR_CEMP:
    COOKIE_CHECK(distr,CK_DISTR_CEMP,/*void*/);
    if (distr->data.cemp.sample) free( distr->data.cemp.sample );
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_DISTR_UNKNOWN,"");
  }

  /* derived distribution ? */
  free( distr->base );

  free( distr );

} /* end of unur_distr_free() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_set_name( struct unur_distr *distr, const char *name )
     /*----------------------------------------------------------------------*/
     /* set name of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   name  ... name of distribution                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  distr->name = name;
  return 1;
} /* end of unur_distr_set_name() */

/*---------------------------------------------------------------------------*/

const char *
unur_distr_get_name( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get name of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   name of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return distr->name;
} /* end of unur_distr_get_name() */

/*---------------------------------------------------------------------------*/

unsigned int unur_distr_get_type( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get type of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   type of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return (distr->type);
} /* end of unur_distr_get_type() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cont( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate continuous.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate continuous                                     */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_CONT) ? 1 : 0);
} /* end of unur_distr_is_cont() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_discr( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate discrete.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate discrete                                       */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_DISCR) ? 1 : 0);
} /* end of unur_distr_is_discr() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_demp( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is empirical univariate discrete.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate discrete                                       */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_DEMP) ? 1 : 0);
} /* end of unur_distr_is_demp() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cemp( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is empirical univariate continuous,             */
     /* i.e. a sample.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate discrete                                       */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_CEMP) ? 1 : 0);
} /* end of unur_distr_is_cemp() */

/*---------------------------------------------------------------------------*/



