/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr.c                                                           *
 *                                                                           *
 *   manipulate distribution objects                                         *
 *                                                                           *
 *   PARAMETER: struct unur_distr *                                          *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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

#include <unur_source.h>
#include "distr.h"
#include "distr_source.h"

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** routines for all distribution objects                                   **/
/**                                                                         **/
/*****************************************************************************/

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
  if (distr) _unur_distr_free( distr );
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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int len;
  char *name_str;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );

  /* allocate memory for storing string */
  len = strlen(name) + 1;
  name_str = _unur_xrealloc(distr->name_str,len);

  /* copy string */
  memcpy( name_str, name, len );

  /* store string in distribution object */
  distr->name_str = name_str;
  distr->name = name_str;

  return UNUR_SUCCESS;
} /* end of unur_distr_set_name() */

/*---------------------------------------------------------------------------*/

const char *
unur_distr_get_name( const struct unur_distr *distr )
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
  _unur_check_NULL( NULL, distr, NULL );

  return distr->name;
} /* end of unur_distr_get_name() */

/*---------------------------------------------------------------------------*/

int
unur_distr_get_dim( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get number of components of random vector (i.e. its dimension)       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   dimension                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );

  return distr->dim;
} /* end of unur_distr_get_dim() */

/*---------------------------------------------------------------------------*/

unsigned int 
unur_distr_get_type( const struct unur_distr *distr )
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
  _unur_check_NULL( NULL, distr, 0u );

  return (distr->type);
} /* end of unur_distr_get_type() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cont( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate continuous.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if continuous distribution                               */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_CONT) ? TRUE : FALSE);
} /* end of unur_distr_is_cont() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cvec( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is multivariate continuous.                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if multivariate continuous                               */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_CVEC) ? TRUE : FALSE);
} /* end of unur_distr_is_cvec() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cvemp( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is empirical multivariate continuous.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if empirical multivariate continuous                     */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_CVEMP) ? TRUE : FALSE);
} /* end of unur_distr_is_cvemp() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_matr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is matrix distribution.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if matrix distribution                                   */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_MATR) ? TRUE : FALSE);
} /* end of unur_distr_is_matr() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_discr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate discrete.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TREU  ... if univariate discrete                                   */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_DISCR) ? TRUE : FALSE);
} /* end of unur_distr_is_discr() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cemp( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is empirical univariate continuous,             */
     /* i.e. a sample.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE  ... if univariate discrete                                   */
     /*   FALSE ... otherwise (and in case of an error)                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, FALSE );

  return ((distr->type == UNUR_DISTR_CEMP) ? TRUE : FALSE);
} /* end of unur_distr_is_cemp() */

/*---------------------------------------------------------------------------*/



