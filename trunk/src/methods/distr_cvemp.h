/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_cvemp.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CVEMP  (continuous empirical multivariate distribution)     *
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
   =DISTRIBUTION  CVEMP  [40]  continuous empirical multivariate distribution


   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical multivariate continuous distributions (VCEMP).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cvemp_new( int dim ); 
/* 
   Create a new (empty) object for empirical multivariate continuous
   distribution. @var{dim} is the number of components of the random
   vector (i.e. its dimension). It must be at least 2; otherwise
   unur_distr_cemp_new() should be used to create an object for an
   empirical univariate distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cvemp_set_data( UNUR_DISTR *distribution, double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
   @var{sample} is an array of double arrays of size 
   @code{dim}x@var{n_sample}, where
   @code{dim} is the dimension of the distribution returned by
   unur_distr_get_dim(). 
   The data points must be stored consecutively in @var{sample}.
*/


int unur_distr_cvemp_get_data( UNUR_DISTR *distribution, double **sample );
/* 
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given,
   @code{0} is returned and @var{sample} is set to NULL.
   If successful @var{sample} points to an array of length
   @code{dim}x@code{n_sample}, where
   @code{dim} is the dimension of the distribution returned by
   unur_distr_get_dim() and @code{n_sample} the return value of the
   function.
*/

/* =END */

/*---------------------------------------------------------------------------*/
