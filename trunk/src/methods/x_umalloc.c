/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_umalloc.c                                                       *
 *                                                                           *
 *   allocate memory                                                         *
 *   store allocated blocks in linked list                                   *
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

#include <stdlib.h>

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

void
_unur_add_mblocks( struct unur_mblock **mblocks, void *ptr )
     /*----------------------------------------------------------------------*/
     /* add pointer to allocated block to list                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   mblocks ... pointer to generator object                                */
     /*   ptr ... allocated block                                            */
     /*----------------------------------------------------------------------*/
{
  register struct unur_mblock *new;

  /* get new entry */
  /* (only once. maybe it we could allocate a larger block.
     but this would require some additional book keeping.) */
  new = _unur_malloc( sizeof(struct unur_mblock) );
  COOKIE_SET(new,CK_MBLOCK);

  /* store allocated block */
  new->memptr = ptr;

  /* link new entry into list */
  new->next = *mblocks;
  *mblocks = new;

} /* end of _unur_add_mblocks() */

/*---------------------------------------------------------------------------*/

void
_unur_free_mblocks( struct unur_mblock *block )
     /*----------------------------------------------------------------------*/
     /* free all blocks in list of allocate memory and free (destroy) list   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   blocks ... pointer to list of blocks                               */
     /*----------------------------------------------------------------------*/
{
  register struct unur_mblock *next;

  while (block != NULL) {
    COOKIE_CHECK(block,CK_MBLOCK,/*void*/);
    next = block->next;
    free(block->memptr);  /* free memory block */
    free(block);          /* free postition in list */
    block = next;
  }
} /* end of _unur_free_mblocks() */

/*---------------------------------------------------------------------------*/
