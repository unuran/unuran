/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_cemp.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CEMP  (continuous empirical univariate distribution)        *
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
   Routines for handling empirical univariate continuous distributions (CEMP).
*/

UNUR_DISTR *unur_distr_cemp_new( void );
/* 
   Create a new (empty) object for empirical univariate continuous distribution.
*/

/* Essential parameters */

int unur_distr_cemp_set_data( UNUR_DISTR *distribution, double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
*/


int unur_distr_cemp_get_data( UNUR_DISTR *distribution, double **sample );
/* 
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given,
   @code{0} is returned and @code{sample} is set to NULL.
*/

/*---------------------------------------------------------------------------*/

