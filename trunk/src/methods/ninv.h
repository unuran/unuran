/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ninv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NINV                               *
 *         (Numerical INVersion of cumulative distribution function)         *
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
   =METHOD NINV Numerical INVersion

   =UP  Methods_for_CONT

   =DESCRIPTION
      NINV is the implementation of numerical inversion.
      For finding the root it is possible to choose between
      Newton's method and the regula falsi.
      
      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_ninv_chg_truncated() call.
      
      To speed up the marginal generation time a table with suitable
      starting points can be computed in the setup. Using such a table can be 
      switched on by means of a unur_ninv_set_table() call where the table
      size is given as a parameter. The table is still useful when the
      (truncated) domain is changed often, since it is computed for the
      domain of the given distribution. (It is not possible to enlarge
      this domain.)
      
      As a rule of thumb using such a table is approriate when the number of
      generated points exceeds the table size by a factor of 100.
      
      It is also possible to change the parameters of the given distribution
      by a unur_ninv_chg_pdfparams() call. If a table exists, it will be
      recomputed immediately.

   =END
*/

/*
 =REQUIRED NINV

 cdf, (in addition pdf only in case of Newton's method) 
*/


/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_ninv_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_ninv_set_usenewton( UNUR_PAR *parameters );
/* 
   Switch Newton's method.
*/

int unur_ninv_set_useregula( UNUR_PAR *parameters );
/* 
   Switch regula falsi. (This the default.)
*/

int unur_ninv_set_max_iter( UNUR_PAR *parameters, int max_iter );
/* 
   Set number of maximal iterations.
*/

int unur_ninv_set_x_resolution( UNUR_PAR *parameters, double x_resolution);
/* 
   Set maximal relative error in x.
*/

int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
/* 
   Set starting points.
   If not set, suitable values are chosen automatically.

   @multitable @columnfractions 0.2 0.2 0.6
   @item Newton:       @tab @var{left}:              @tab starting point
   @item Regula falsi: @tab @var{left}, @var{right}: @tab boundary of starting interval
   @end multitable

   If the starting points are not set then the follwing points are used by
   default:
   @multitable @columnfractions 0.2 0.2 0.6
   @item Newton:       @tab @var{left}:  @tab CDF(@var{left}) = 0.5
   @item Regula falsi: @tab @var{left}:  @tab CDF(@var{left}) = 0.1
   @item               @tab @var{right}: @tab CDF(@var{right}) = 0.9
   @end multitable

   If @var{left} == @var{right}, then  UNURAN always uses the default
   starting points!
*/
    
int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
/* 
   Generates a table with @var{no_of_points} points containing
   suitable starting values for the iteration. The value of
   @var{no_of_points} must be at least 10 (otherwise it will be set
   to 10 automatically).

   The table points are chosen such that the CDF at these points
   form an equidistance sequence in the interval (0,1).

   If a table is used, then the starting points given by
   unur_ninv_set_start() are ignored.
 */

int unur_ninv_chg_max_iter(UNUR_GEN *generator, int max_iter);
/* 
   Change the maximum number of iterations of an inversion step.
*/

int unur_ninv_chg_x_resolution(UNUR_GEN *generator, double x_resolution);
/*
  Change the maximal relative error in x.
*/

int unur_ninv_chg_start(UNUR_GEN *gen, double left, double right);
/* 
   Change the starting points for numerical inversion. 

   If left==right, then UNURAN uses the default starting points 
   (see unur_ninv_set_start()).
*/

int unur_ninv_chg_table(UNUR_GEN *gen, int no_of_points);
/*
   Recomputes a table as described in unur_ninv_set_table().
*/

int unur_ninv_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Change the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call.
   Moreover the starting point(s) will not be changed.
*/

int unur_ninv_chg_pdfparams(UNUR_GEN *generator, double *params, int n_params);
/*
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of
   distribution parameters. If a table is used, it will be computed
   immediately.

   @emph{IMPORTANT:} The given parameters are not checked against
   domain errors; in opposition to the 
   @command{unur_<distr>_new} calls.
*/ 


/* =END */
/*---------------------------------------------------------------------------*/
