/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_rect.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method RECT       *
 *         (uniformly distributed in (multidimensional) RECTangle            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only used in unur_methods.h                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_rect_par { 
  int      dim;         /* dimension                                         */
  double **domain;      /* vertices of rectangle                             */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_rect_gen { 
  int      dim;         /* dimension                                         */
  double **domain;      /* vertices of rectangle                             */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_rect_new( int dim );
/* get default parameters for generator                                      */

struct unur_gen *unur_rect_init( struct unur_par *parameters );
/* initialize new generator                                                  */

void unur_rect_sample_vec( struct unur_gen *gen, double *vec );
/* sample from generator                                                     */

void unur_rect_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/

