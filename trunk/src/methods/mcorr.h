/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mcorr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MCORR                              *
 *         (Matrix - CORRelation matrix)                                     *
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
   =METHOD  MCORR   Random CORRelation matrix

   =UP  Methods_for_MATR

   =REQUIRED  Distribution object for random correlation matrix

   =OPTIONAL 

   =SPEED Set-up: fast,
          Sampling: depends on dimension

   =REF  [DLa86: Sect.6.1; p.605]

   =DESCRIPTION
      MCORR generates a random correlation matrix.
      Thus a matrix @unurmath{H} is generated where all rows are
      independent random vectors of unit length uniformly on a sphere.
      Then @unurmath{HH'} is a correlation matrix (and vice versa if
      @unurmath{HH'} is a correlation matrix then the rows of
      @unurmath{H} are random vectors on a sphere).
      There are many other possibilites (distributions) of sampling
      the random rows from a sphere. The chosen one is simple but
      does in not result in a uniform distriubution of the random
      correlation matrices.

      It only works with distribution objects of random correlation
      matrices (@pxref{correlation,,Random Correlation Matrix}).

   =HOWTOUSE
      Create a distibution object for random correlation matrices by a
      @code{unur_distr_correlation} call
      (@pxref{correlation,,Random Correlation Matrix}).
      Notice that due to round-off errors,
      there is an (extremely small) chance that the resulting matrix is
      not positive definite for a Cholesky decomposition algorithm.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mcorr_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/* =END */
/*---------------------------------------------------------------------------*/


