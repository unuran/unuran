/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen.c                                                           *
 *                                                                           *
 *   miscelleanous routines for manipulation generator objects               *
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

/*****************************************************************************/
/**                                                                         **/
/**  Call Init, Sampling, and Free functions                                **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

UNUR_GEN *unur_init( UNUR_PAR *par )
{                
  _unur_check_NULL(NULL,par,NULL);
  return (par->init(par));
} /* end of unur_init() */

/*---------------------------------------------------------------------------*/

int unur_reinit( UNUR_GEN *gen )
{               
  /* if gen->reinit == NULL   --> nothing to do */
  return (gen ? gen->reinit(gen) : 1 );
} /* end of unur_reinit() */

int _unur_reinit_error( UNUR_GEN *gen )
{
  /* no re-initialization routine available */
  _unur_warning(gen->genid,UNUR_ERR_GENERIC,"no re-initialization routine available");
  return 0;
} /* end of _unur_reinit_error() */

/*---------------------------------------------------------------------------*/

int unur_sample_discr(UNUR_GEN *gen)
{
  CHECK_NULL(gen,0);
  return (gen->sample.discr(gen));
} /* end of unur_sample_discr() */

double unur_sample_cont(UNUR_GEN *gen)
{
  CHECK_NULL(gen,0.);
  return (gen->sample.cont(gen));
} /* end of unur_sample_cont() */

void unur_sample_vec(UNUR_GEN *gen, double *vector)
{
  CHECK_NULL(gen,/*void*/);
  gen->sample.vec(gen,vector);
} /* end of unur_sample_vec() */

/*---------------------------------------------------------------------------*/

void unur_free( UNUR_GEN *gen )
{                
  CHECK_NULL(gen,/*void*/);
  gen->destroy(gen);
} /* end of unur_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Get data about generator object                                        **/
/**                                                                         **/
/*****************************************************************************/

int
unur_get_dimension( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get dimension of generator for multivariate distribution             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   dimension of distribution                                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  switch (gen->method) {
  case UNUR_METH_RECT:
    COOKIE_CHECK(gen,CK_RECT_GEN,0);
    return gen->data.rect.dim;
    break;
  default:
    /* method unknown assume dim = 1 */
    /** TODO: distinguish between univariate --> dim = 1 and
	unknown multivariate **/
    return 1.;
  }

} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/

const char *
unur_get_genid( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get generator id                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator id                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  return gen->genid;
} /* end of unur_get_genid() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_get_distr( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object from generator object             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);

  return &(gen->distr);
} /* end of unur_get_distr() */

/*---------------------------------------------------------------------------*/
