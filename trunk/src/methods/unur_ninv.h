/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_ninv.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method NINV       *
 *         (Numerical inversion of cumulative distribution function)         *
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

struct unur_ninv_par { 
  int max_inter;             /* maximal number of iterations                 */
  int rel_x_resolution;      /* maximal relative error in x                  */
  double sl, sr;             /* interval boundaries at start (left/right)    */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ninv_gen { 
  double  um;                /* height of rectangle: square root of f(mode)  */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_ninv_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_ninv_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_ninv_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_ninv_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_ninv_set_usemirror( struct unur_par *par, int usemirror );
/* set flag for using mirror principle (default: off)                        */

#define unur_ninv_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

