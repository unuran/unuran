/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cemp.h                                                            *
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
   =NODEX   CEMP   Continuous empirical univariate distributions

   =UP Distribution_objects [20]

   =DESCRIPTION
      Empirical univariate distributions are just lists of numbers.
      Thus there are only calls to insert these data.
      How these data are used to sample from the empirical distribution
      depends from the chosen generation method.  

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical univariate continuous distributions (CEMP).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cemp_new( void );
/* 
   Create a new (empty) object for empirical univariate continuous distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cemp_set_data( UNUR_DISTR *distribution, const double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
*/

int unur_distr_cemp_read_data( UNUR_DISTR *distribution, const char *filename );
/* 
   Read data from file @file{filename}.
   It reads the first number from each line. 
   Numbers are parsed by means of the C standard routine @command{strtod}.
   Lines that do not start with @code{+}, @code{-}, @code{.}, or a
   digit are ignored. (Beware of lines starting with a blank!)

   In case of an error (file cannot be opened, invalid string for
   double in line) no data are copied into the distribution object
   and an error code is returned.
*/

int unur_distr_cemp_get_data( const UNUR_DISTR *distribution, const double **sample );
/* 
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given, an error code 
   is returned and @code{sample} is set to NULL.

   @emph{Important:} Do @strong{not} change the entries in @var{sample}!
*/

/* =END */

/*---------------------------------------------------------------------------*/
