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
  int max_iter;              /* maximal number of iterations                 */
  double rel_x_resolution;   /* maximal relative error in x                  */
  double s[3];               /* interval boundaries at start (left/right)    */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ninv_gen { 
  int max_iter;              /* maximal number of iterations                 */
  double rel_x_resolution;   /* maximal relative error in x                  */
  double s[3];               /* interval boundaries at start (left/right)    */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_ninv_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_ninv_init( struct unur_par *parameters );
/* initialize new generator                                                  */

double unur_ninv_sample_regula( struct unur_gen *generator );
double unur_ninv_sample_newton( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_ninv_free( struct unur_gen *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_ninv_use_newton( struct unur_par *par );
/* use Newton's method                                                       */

int unur_ninv_use_regula( struct unur_par *par );
/* use regula falsi                                                          */

int unur_ninv_set_max_iter( struct unur_par *par, int max_iter );
/* set number of maximal iterations                                          */

int unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution);
/* set maximal relative error in x                                           */

int unur_ninv_set_start( struct unur_par *par, double s1, double s2, double s3 );
/*---------------------------------------------------------------------------*/
/* set starting points.                                                      */
/*   Newton:        s1           starting point                              */
/*   regular falsi: s1, s2       boundary of starting interval               */
/*   Muller/Brent:  s1. s2, s3   starting points                             */
/* arguments that are used by method are ignored.                            */

#define unur_ninv_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

