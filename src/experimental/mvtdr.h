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
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

   =REINIT not implemented

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

int unur_mvtdr_set_boundsplitting( UNUR_PAR *par, double boundsplitting );
/* 
   Set bound for splitting cones. All cones are split which have a
   volume below the hat is greater than @var{bound_splitting} times
   the average over all volumes. However, the number given by the 
   unur_mvtdr_set_maxcones() is not exceeded.
   Notice that the later number is always reached 
   if @var{bound_splitting} is less than 1.

   Default: @code{1.5}
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

int unur_mvtdr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_mvtdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/


