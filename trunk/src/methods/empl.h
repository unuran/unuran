/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: empl.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method EMPD                               *
 *         (EMPirical distribution with Linear interpolation)                *
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
   =METHOD  EMPL   EMPirical distribution with Linear interpolation

   =UP  Methods_for_CEMP

   =REQUIRED observed sample

   =SPEED Set-up: slow (as sample is sorted),
          Sampling: very fast (inversion) 

   =REF  [HLa00]

   =DESCRIPTION
      EMPL generates random variates from an empirical distribution
      that is given by an observed sample. This is done by linear
      interpolation of the empirical CDF. Although this
      method is suggested in the books of Law and Keltn (2000) and
      Bratly, Fox, and Schrage (1987) we do not recommend this method at
      all since it has many theoretical drawbacks:
      The variance of empirical distribution function does not
      coincide with the variance of the given sample. Moreover,
      when the sample increases the empirical density function
      does not converge to the density of the underlying random
      variate. Notice that the range of the generated point set is
      always given by the range of the given sample. 

      This method is provided in UNURAN for the sake of
      completeness. We always recommend to use method EMPK
      (@pxref{EMPK,,EMPirical distribution with Kernel smoothing}).

      If the data seem to be far away from having a bell shaped
      histogram, then we think that naive resampling is still better
      than linear interpolation.

      @noindent
      @emph{Important}: Using this method is not recommended!

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_empl_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/* =END */

/*---------------------------------------------------------------------------*/
