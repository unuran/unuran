/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mcgibbs.h                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MCGIBBS                            *
 *         (Markov Chain - GIBBS sampler)                                    *
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

/* 
   =TEST-METHOD  MCGIBBS   Markov Chain - GIBBS sampler

   =UP  Methods_for_CVEC

   =REQUIRED log-concave PDF, dPDF

   =SPEED Set-up: fast, Sampling: moderate

   =REF  [HLD04: Sect.14.1.2]

   =DESCRIPTION
      The Gibbs sampler is an implementation of a coordinate sampler,
      wherein, the sampling coordinate directions are changed cyclically.
      For each direction, the full conditional through the actual
      sampling point is being used in connection with the Gilks&Wild
      variant of the Transformed Density Rejection method.

      TODO.

   =HOWTOUSE
      For using the GIBBS method UNURAN needs the PDF of the
      distribution. 
      The Gibbs sampler should only be used for weak correlation structure 
      of the PDF. Strong correlations between the independent variates
      often lead to very weak ergodicity properties when using this 
      method.

      TODO.

      @strong{Warning:} If the ... debugging flag is set to verbose (other than @code{1u}) 
      every Gibbs step might produces a long entry in the log file ...

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mcgibbs_new( const UNUR_DISTR *distribution );

/*...........................................................................*/

int unur_mcgibbs_set_variant_coordinate( UNUR_PAR *parameters );
/* 
   Coordinate Sampler :
   Sampling along the coordinate directions (cyclic).
   This is the default.
*/

int unur_mcgibbs_set_variant_random_direction( UNUR_PAR *parameters );
/* 
   Random Direction Sampler :
   Sampling along the random directions.
   Not implemented at the moment.
*/

int unur_mcgibbs_set_c( UNUR_PAR *parameters, double c );
/* 
   TODO.

   Set parameter @var{c} for transformation T. 
   Currently only values between 0 and -0.5 are allowed.
   If @code{c} is between 0 and -0.5 it is set to -0.5.

   For @code{c=0} much more robust against badly normalized PDFs
   than $code{c=-0.5}.

   @code{c=0}: method TDRGW
   @code{c=-0.5}: method TDR

   other values not implemented yet.


   Default is @code{0}.
*/


int unur_mcgibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0);
/* 
   Starting point.
   TODO.
*/

int unur_mcgibbs_set_thinning( UNUR_PAR *parameters, int thinning );
/*
   Sets the parameter @var{thinning}.
   A thinning parameter of n means, that the interval size between two 
   sampled points (returned random vectors) is n, whereas the internally 
   sampled points are being a unit interval apart.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=1.

   Default: @code{1}.
*/

/*...........................................................................*/

/* =END */
/*---------------------------------------------------------------------------*/


