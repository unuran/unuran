/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cmat.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CMAT  (continuous matrix distribution)                      *
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
   =NODE   CMAT   Continuous MATrix distributions

   =UP Distribution_objects [45]

   =DESCRIPTION
      Distributions for random matrices. Notice that UNURAN uses
      arrays of @code{double}s to handle matrices. There the rows of
      the matrix are stored consecutively.

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling continuous matrix distributions (CMAT).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cmat_new( int n_rows, int n_cols );
/* 
   Create a new (empty) object for a continuous matrix
   distribution. @var{n_rows} and @var{n_cols} are the respective
   numbers of rows and columns of the random matrix (i.e. its
   dimensions). Each must be at least 2; otherwise
   unur_distr_cont_new() or unur_distr_cvec_new() should be used to
   create an object for a univariate distribution and a multivariate
   (vector) distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cmat_get_dim( const UNUR_DISTR *distribution, int *n_rows, int *n_cols );
/* 
   Get number of rows and columns of random matrix (its dimension).
   It returns the total number of components. In case of an error
   @code{0} is returned.
*/

/* =END */


/*---------------------------------------------------------------------------*/

/* not implemented: */
/* DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the
   essential parameters has been set or changed (and the parameter is
   required for the chosen method).
*/

/*---------------------------------------------------------------------------*/
