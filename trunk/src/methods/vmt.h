/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vmt.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VMT                                *
 *         (Vector Matrix Transformation)                                    *
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
   =METHOD  VMT   Vector Matrix Transformation

   =UP  Methods_for_CVEC

   =REQUIRED mean vector, covariance matrix, standardized marginal distributions

   =SPEED Set-up: slow,
          Sampling: depends on dimension

   =DESCRIPTION
      VMT generates random vectors for distributions with given mean
      vector mu and covariance matrix Sigma. It produces random vectors
      of the form X = L Y + mu, where L is the Cholesky factor of Sigma,
      i.e. L L^t = Sigma, and Y has independent components of the same
      distribution with mean 0 and standard deviation 1.
      
      The method VMT has been implemented especially to sample from a
      multinormal distribution. Nevertheless, it can also be used (or
      abused) for other distributions. However, notice that the given
      standardized marginal distributions are not checked; i.e.
      if the given distributions do not have mean 0 and variance 1
      then mu and Sigma are not the mean vector and covariance matrix,
      respectively, of the resulting distribution. 

      @strong{Important:} Notice that except for the multinormal
      distribution the given marginal distribution are distorted by
      the transformation using the Cholesky matrix. Thus for other
      (non-multinormal) distributions this method should only be used
      when everything else fails and some approximate results which
      might even be not entirely correct are better than no results.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vmt_new( const UNUR_DISTR *distribution );
/* 
   Get parameters for generator.
*/

/* =END */
/*---------------------------------------------------------------------------*/


