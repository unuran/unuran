/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: gibbs.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method GIBBS                              *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************

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
   =TEST-METHOD  GIBBS   Implementation of the Gibbs sampler using the full-conditionals.

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode

   =SPEED Set-up: fast, Sampling: slow

   =DESCRIPTION
      TODO ...   
   
   =HOWTOUSE
      For using the GIBBS method UNURAN needs the PDF of the
      distribution. Additionally, the parameter @i{r} can be set via
      a unur_gibbs_set_r() call.


   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_gibbs_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_gibbs_set_skip( UNUR_PAR *parameters, long skip );
/*
   Sets the parameter @var{skip} i.e. the number of sampling
   steps between two points that will be used as returned 
   random numbers

   @emph{Notice}: This parameter must satisfy @var{skip}>=0.

   Default: @code{0}.
*/

/*...........................................................................*/

long _unur_gibbs_get_pdfcount( UNUR_GEN *gen);
/* Return the number of PDF calls */

void _unur_gibbs_reset_pdfcount( UNUR_GEN *gen);
/* Reset the number of PDF calls to 0 */

