/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: sinv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method SINV                               *
 *         (Spline approximation for INVerse of CDF)                         *
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
   =METHOD SINV  Spline approximation for INVerse CDF

   =UP  Methods_for_CONT

   =REQUIRED CDF, PDF

   =SPEED Set-up: (very) slow, Sampling: fast

   =DESCRIPTION
      SINV is a variant of numerical inversion, where the inverse CDF
      is approximated using splines. These splines have to be computed 
      in a setup step. However, it only works for distributions with 
      bounded domain (i.e., for distributions with unbounded domain the 
      tails have to be chopped off using the unur_sinv_set_boundary()
      call. If these boundaries are not set, SINV tries to find out
      proper boundaries.
      
      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_sinv_chg_truncated() call.

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      can be made as small as desired by means of the
      unur_sinv_set_u_resolution(). Notice that this increased the cost 
      of the setup step.
      As a rule of thumb reducing the error of a factor 100 increased
      the computed table and thus the setup time by a factor of 10.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_sinv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_sinv_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set the left and right boundary of the computation interval.
   The probability outside of this region must not be of
   computational relevance.
   Of course @code{+/- UNUR_INFINITY} is not allowed.

   Default is @code{1.e20}.
*/

int unur_sinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal error in .... direction.
   Default is @code{10^-8}.
*/

int unur_sinv_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Changes the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call. The tables splines are not recomputed. Thus
   it might happen that the relative error for the generated variates
   from the truncated distribution is greater than the bound given by a
   unur_sinv_chg_truncated() call.
*/

/* =END */
/*---------------------------------------------------------------------------*/
