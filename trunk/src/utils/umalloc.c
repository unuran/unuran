/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: umalloc.c                                                         *
 *                                                                           *
 *   allocate memory                                                         *
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

/*---------------------------------------------------------------------------*/
#ifdef WITH_DMALLOC
/*---------------------------------------------------------------------------*/
/* we use a macro, see umalloc_source.h                                      */
/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/

void*
_unur_malloc(size_t size)
     /*----------------------------------------------------------------------*/
     /* allocate memory                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   size ... size of allocated block                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   abort program                                                      */
     /*----------------------------------------------------------------------*/
{
  register void *ptr;

  /* allocate memory */
  ptr = malloc( size );

  /* successful ? */
  if (ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }

  return ptr;

} /* end of _unur_malloc() */

/*---------------------------------------------------------------------------*/

void*
_unur_realloc(void *ptr, size_t size)
     /*----------------------------------------------------------------------*/
     /* reallocate memory                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ptr  ... address of memory block previously allocated by malloc.   */
     /*   size ... size of reallocated block                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   abort program                                                      */
     /*----------------------------------------------------------------------*/
{
  register void *new_ptr;

  /* reallocate memory */
  
  new_ptr = realloc( ptr, size );

  /* successful ? */
  if (new_ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }

  return new_ptr;

} /* end of _unur_realloc() */

/*---------------------------------------------------------------------------*/
#endif   /* WITH_DMALLOC */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Function prototypes for allocating new generator objects                  */
struct unur_gen *
_unur_malloc_gen( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for new generator object                             */
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

} /* end of _unur_malloc_gen() */

/*---------------------------------------------------------------------------*/
