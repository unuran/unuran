/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: walk.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method WALK                               *
 *         Metropolis random walk sampler                                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************

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

/*
   =METHOD  WALK   Metropolis random-walk sampler

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: fast, Sampling: fast

   =DESCRIPTION
      WALK is an implementation of Metropolis random walk sampler.
      
   =HOWTOUSE
      For using the WALK method UNURAN needs the PDF of the
      distribution. 

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_walk_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/



int unur_walk_set_ball_radius( UNUR_PAR *par, double ball_radius );
/* 
   Sets radius of ball used for the random walk sampler. 

   Default: @code{1.0}
*/


int unur_walk_set_thinning( UNUR_PAR *par, long thinning );
/*
   Sets the parameter @var{thinning}.
   A thinning parameter of n means, that the interval size between two 
   sampled points (returned random vectors) is n, whereas the internally 
   sampled points are being a unit interval apart.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=1.

   Default: @code{1}.
*/


/* =END */

/*---------------------------------------------------------------------------*/

void _unur_walk_set_point_current( UNUR_GEN *gen, double *x);
/* set the current point (dimension=dim) */


void _unur_walk_get_point_current( UNUR_GEN *gen, double *x);
/* get the current point (dimension=dim) */
