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

/*...........................................................................*/

int unur_rect_set_domain_vec( struct unur_par *par, double **domain );
/* set coordinates for domain boundary                                       */

#define unur_rect_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

