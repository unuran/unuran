/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: functparser_source.h                                              *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares prototypes for function parser                           *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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
#ifndef __FSTR_SOURCE_H_SEEN
#define __FSTR_SOURCE_H_SEEN

/*---------------------------------------------------------------------------*/
/* Function prototypes for function string parser                            */
/*---------------------------------------------------------------------------*/

struct ftreenode *_unur_fstr2tree( char *functstring );
/*---------------------------------------------------------------------------*/
/* Compute funtion tree from string.                                         */
/*---------------------------------------------------------------------------*/

double _unur_fstr_eval_tree( struct ftreenode *functtree_root, double x );
/*---------------------------------------------------------------------------*/
/* Evalutes function given by a function tree at x.                          */
/*---------------------------------------------------------------------------*/

void _unur_fstr_free( struct ftreenode *functtree_root);
/*---------------------------------------------------------------------------*/
/* Destroys function tree and frees memory.                                  */
/*---------------------------------------------------------------------------*/

char *Ntree2string( struct ftreenode *functtree_root, char *ret_str );
/*---------------------------------------------------------------------------*/
/* Creates string for program code.                                          */
/*---------------------------------------------------------------------------*/

struct ftreenode *_unur_fstr_make_derivative( struct ftreenode *functtree_root );
/*---------------------------------------------------------------------------*/
/* Make function tree for derivate of given function (tree).                 */ 
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* __FSTR_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/

