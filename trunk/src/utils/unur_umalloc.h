/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_umalloc.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for allocating        *
 *         memory and storing allocated blocks in a linked list.             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file.                                             *
 *         included in umalloc.c, error.c and ../methods/unur_methods.h.     *
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
#ifndef __UNUR_UMALLOC_H_SEEN
#define __UNUR_UMALLOC_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>

/*---------------------------------------------------------------------------*/
/* linked list to store pointers to allocated blocks                         */

struct unur_mblock {
  void *memptr;                 /* store pointer to allocated block          */
  struct unur_mblock *next;     /* pointer to next entry in list             */
#if UNUR_DEBUG & UNUR_DB_COOKIES
  unsigned cookie;              /* magic cookie                              */
#endif
};

/*---------------------------------------------------------------------------*/
/* Function prototypes                                                       */

void *_unur_malloc(size_t size);
void  _unur_add_mblocks( struct unur_mblock **mblocks, void *ptr );
void  _unur_free_mblocks( struct unur_mblock *mblocks );

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_UMALLOC_H_SEEN */
/*---------------------------------------------------------------------------*/






