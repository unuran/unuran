/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: condi.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         id   CONDI  (continuous full conditional distribution)            *
 *         type CONT   (continuous univariate distribution)                  *
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

/* 
   =NODEX   CONDI   Continuous univariate full conditional distribution

   =UP Distribution_objects [35]

   =DESCRIPTION
      Full conditional distribution for a given continuous
      multivariate distributiion.

      This is a special case of a continuous univariate distribution
      and thus they have most of these parameters (with the exception
      that functions cannot be changed). Additionally,

      @itemize @minus
      @item there is a call to extract the underlying multivariate
            distribution,

      @item and a call to handle the variables that are fixed.

      @end itemize

      This distibution type is primarily used for evaluation the
      conditional distribution and its derivative (as required for,
      e.g., the Gibbs sampler). The density is not normalized (does
      not integrate to one). Mode and area are not available and it
      does not make sense to use any call to set or change parameters
      except the ones given below.

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous full conditional
   distributions (CONDI).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_condi_new( const UNUR_DISTR *distribution, const double *pos, int k );
/* 
   Create an object for full conditional distribution for the given
   @var{distribution} for the @var{k}-thvariable and the other
   variables fixed to @var{pos}.

   @var{distribution} must be a pointer to a univariate continuous
   distribution. 
   @var{pos} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension of the underlying distribution object.
   @var{k} must be in the range @code{0, @dots{}, dim-1}.

   The resulting generator object is of the same type as of a
   unur_distr_cont_new() call.
*/

int unur_distr_condi_set_condition( struct unur_distr *distribution, const double *pos, int k );
/* 
   Set/change condition for conditional @var{distribution}. 
   Change values of fixed variables to @var{pos} and use @var{k}-th
   variable of conditional @var{distribution}.

   @var{pos} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension of the underlying distribution object.
   @var{k} must be in the range @code{0, @dots{}, dim-1}.
*/

int unur_distr_condi_get_condition( struct unur_distr *distribution, const double **pos );
/* 
   Get condition for conditional @var{distribution}. 
   The values for the fixed variables are stored in @var{pos}, which
   must be a pointer to an array of size @code{dim}.
   The function returns the number of the variable for which the 
   conditional distribution is computed.

   @emph{Important:} Do @strong{not} change the entries in @var{pos}!
*/

const UNUR_DISTR *unur_distr_condi_get_distribution( const UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/

/* =END */

/*---------------------------------------------------------------------------*/
