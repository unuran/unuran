/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mvtdr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MVTDR                              *
 *         (Multi-Variate Transformed Density Rejection)                     *
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
   =EXPERIMENTAL-METHOD  MVTDR  Multi-Variate Transformed Density Rejection

   =UP  Methods_for_CVEC

   =REQUIRED log-concave (log)PDF, gradient of (log)PDF

   =OPTIONAL mode

   =SPEED Set-up: slow,
          Sampling: depends on dimension

   =REF  [HLD04: Sect.11.3.4, Alg.11.15.] [LJa98]

   =DESCRIPTION
      MVTDR

   =HOWTOUSE
      Create a multivariate generator object ...

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mvtdr_new( const UNUR_DISTR *distribution );
/* 
   Get parameters for generator.
*/

/*...........................................................................*/

int unur_mvtdr_set_stepsmin( UNUR_PAR *par, int stepsmin );
/* 
   Set minimum number of triangulation step for each starting cone.
   @var{stepsmin} must be nonnegative.

   Default: @code{5}.
*/

int unur_mvtdr_set_maxcones( UNUR_PAR *par, int maxcones );
/* 
   Set maximum number of cones. 

   Notice that this number is always increased to 
   @unurmath{2^{dim+stepsmin}} where @i{dim} is the dimension of the
   distribution object and @i{stepsmin} the given mimimum number of
   triangulation steps.

   Notice: For higher dimensions and/or higher correlations between the
   coordinates of the random vector the required number of cones can
   be very high. A too small maximum number of cones can lead to 
   a very high rejection constant.

   Default: @code{10000}.
*/

/* =END */
/*---------------------------------------------------------------------------*/


