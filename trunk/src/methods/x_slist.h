/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_slist.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines function prototypes for simple list                       *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
#ifndef X_SLIST_H_SEEN
#define X_SLIST_H_SEEN
/*---------------------------------------------------------------------------*/
/* Not part of manual!                                                       */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* A simple list can be used to store an arbitrary numbers of pointers       */
/* to allocated memory in a list.                                            */
/*                                                                           */
/* IMPORTANT: These elements must be allocated via (c|m|re)alloc()!!         */
/*                                                                           */
/*---------------------------------------------------------------------------*/

struct unur_slist *_unur_slist_new( void );
/*---------------------------------------------------------------------------*/
/* Make new simple list.                                                     */
/*---------------------------------------------------------------------------*/

void _unur_slist_append( struct unur_slist *slist, void *element );
/*---------------------------------------------------------------------------*/
/* Append pointer to element to simple list.                                 */
/*---------------------------------------------------------------------------*/

int _unur_slist_length( struct unur_slist *slist );
/*---------------------------------------------------------------------------*/
/* Get length if list (number of list entries).                              */
/*---------------------------------------------------------------------------*/

void *_unur_slist_get( struct unur_slist *slist, int n );
/*---------------------------------------------------------------------------*/
/* Get pointer to n-th element.                                              */
/*---------------------------------------------------------------------------*/

void _unur_slist_free( struct unur_slist *slist );
/*---------------------------------------------------------------------------*/
/* Free all elements and list in simple list.                                */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* X_SLIST_H_SEEN */
/*---------------------------------------------------------------------------*/
