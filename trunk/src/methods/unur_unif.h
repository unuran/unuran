/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_unif.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method UNIF       *
 *         (passes UNIForm random numbers through UNURAN framework;          *
 *         for testing only)                                                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only used in unur_methods.h                                       *
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
/* Information for constructing the generator                                */

struct unur_unif_par { 
  int start;            /* starting point in sequence (0 = first number)     */
  int skip;             /* skip for subsequence (1 = full sequence)          */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_unif_gen { 
  int start;            /* starting point in sequence                        */
  int skip;             /* skip for subsequence                              */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_unif_new( int start, int skip );
/* get default parameters for generator                                      */

struct unur_gen *unur_unif_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_unif_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_unif_free( struct unur_gen *generator );
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/

