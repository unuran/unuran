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

#include <unur_source.h>
#include <distr/distr_source.h>
#include "x_gen.h"
#include "x_gen_source.h"

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
  CHECK_NULL(gen,RETURN_VOID);
  gen->sample.cvec(gen,vector);
} /* end of unur_sample_vec() */

void unur_sample_matr(UNUR_GEN *gen, double *matrix)
{
  CHECK_NULL(gen,RETURN_VOID);
  gen->sample.matr(gen,matrix);
} /* end of unur_sample_matr() */

/*---------------------------------------------------------------------------*/
/* aux routines when no sampling routine is available                         */

double _unur_sample_cont_error( UNUR_GEN *gen )
{
  /* no sampling routine available */
  return INFINITY;
} /* end of _unur_sample_cont_error() */

/*---------------------------------------------------------------------------*/

void unur_free( UNUR_GEN *gen )
{                
  if (gen) gen->destroy(gen);
} /* end of unur_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Get data about generator object                                        **/
/**                                                                         **/
/*****************************************************************************/

int
unur_get_dimension( const struct unur_gen *gen )
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

  return (gen->distr->dim);
} /* end of unur_get_dimension() */

/*---------------------------------------------------------------------------*/

const char *
unur_get_genid( const struct unur_gen *gen )
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
  CHECK_NULL(gen,NULL);

  return gen->genid;
} /* end of unur_get_genid() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_get_distr( const struct unur_gen *gen )
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
  CHECK_NULL(gen,NULL);

  return gen->distr;
} /* end of unur_get_distr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Copy (clone) generator object                                          **/
/**                                                                         **/
/*****************************************************************************/

struct unur_gen *
unur_gen_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "Clone", gen, NULL );
  _unur_check_NULL( "Clone", gen->clone, NULL );

  return (gen->clone(gen));
} /* end of unur_get_clone() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  create and Copy (clone) generator objects                              **/
/**                                                                         **/
/*****************************************************************************/

struct unur_gen *
_unur_generic_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* create new generic generator object                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* copy distribution object into generator object */
  gen->distr = (par->distr) ? _unur_distr_clone(par->distr) : NULL;

  /* copy some parameters into generator object */
  gen->method = par->method;        /* indicates method and variant          */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */
  gen->urng_aux = par->urng_aux;    /* pointer to auxilliary URNG            */

  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_list = NULL;         /* no auxilliary generator objects       */

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_generic_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_generic_clone( const struct unur_gen *gen, const char *type )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generic generator object                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   type ... type of generator (string)                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *clone;

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(type);

  /* copy distribution object into generator object */
  if (gen->distr) clone->distr = _unur_distr_clone(gen->distr);

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = _unur_gen_clone( gen->gen_aux );

  /* finished clone */
  return clone;
} /* _unur_generic_clone() */
/*---------------------------------------------------------------------------*/
