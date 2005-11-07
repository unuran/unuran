/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: gibbs.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method GIBBS                              *
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
   =TEST-METHOD  GIBBS   Implementation of the Gibbs sampler using the full-conditionals.

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode

   =SPEED Set-up: fast, Sampling: moderate

   =DESCRIPTION
      The Gibbs sampler is an implementation of a coordinate sampler,
      wherein, the sampling coordinate directions are changed cyclically.
      For each direction, the full conditional through the actual
      sampling point is being used in connection with the Gilks&Wild
      variant of the Transformed Density Rejection method.
      
   =REF  [HLD04: Cha.14] [GWa92]      
      
   =HOWTOUSE
      For using the GIBBS method UNURAN needs the PDF of the
      distribution. 
      The Gibbs sampler should only be used for weak correlation structure 
      of the PDF. Strong correlations between the independent variates
      often lead to very weak ergodicity properties when using this 
      method.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_gibbs_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_gibbs_set_skip( UNUR_PAR *parameters, long skip );
/*
   Sets the parameter @var{skip} i.e. the number of sampling
   steps between two points that will be used as returned 
   random numbers

   @emph{Notice}: This parameter must satisfy @var{skip}>=0.

   Default: @code{0}.
*/

/*...........................................................................*/

int unur_gibbs_set_variant_coordinate( UNUR_PAR *par );
/* 
   Coordinate Sampler :
   Sampling along the coordinate directions (cyclic).
   This is the default.
*/

/*...........................................................................*/

int unur_gibbs_set_variant_random_direction( UNUR_PAR *par );
/* 
   Random Direction Sampler :
   Sampling along the random directions.
   Not implemented at the moment.
*/

/* =END */

/*...........................................................................*/

void _unur_gibbs_set_point_current( UNUR_GEN *gen, double *x);
/* set the current point (dimension=dim) */

/*---------------------------------------------------------------------------*/

void _unur_gibbs_get_point_current( UNUR_GEN *gen, double *x);
/* get the current point (dimension=dim) */
