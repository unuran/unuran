/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: umalloc_source.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         prototypes for allocating memory blocks                           *
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
/* Function prototypes for allocating memory blocks                          */

/*---------------------------------------------------------------------------*/
#ifdef WITH_DMALLOC
/*---------------------------------------------------------------------------*/
#define _unur_malloc(size)        xmalloc((size))
#define _unur_realloc(ptr,size)   xrealloc((ptr),(size))
/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
void *_unur_malloc(size_t size);
void *_unur_realloc(void *ptr, size_t size);
/*---------------------------------------------------------------------------*/
#endif   /* WITH_DMALLOC */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Function prototypes for allocating new generator objects                  */
struct unur_gen *_unur_malloc_gen( struct unur_par *par );

/*---------------------------------------------------------------------------*/





