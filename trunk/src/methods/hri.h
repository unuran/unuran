/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hri.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HRI                                *
 *         (Hazard Rate Increasing)                                          *
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
   =METHOD  HRI   Hazard Rate Increasing

   =UP  Methods_for_CONT

   =REQUIRED increasing (non-decreasing) hazard rate 

   =SPEED Set-up: fast, Sampling: slow

   =DESCRIPTION
      Generates random variate with given non-increasing hazard rate. 
      It is necessary that the distribution object contains this hazard rate.
      Increasing hazard rate implies that the corresponding PDF of the
      distribution has heavier tails than the exponential distribution
      (which has constant hazard rate).

      It is important to note that the domain of the distribution can
      be set via a unur_distr_cont_set_domain() call. However, only
      the left hand boundary is used. For computational reasons the
      right hand boundary is always reset to @code{UNUR_INFINITY}.
      If no domain is given by the user then the left hand boundary is
      set to @code{0}.
      
      For distributions with decreasing hazard rate method HRD is
      required.
      For distributions which do not have increasing or decreasing
      hazard rates but are bounded from above use method HRB. 

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hri_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hri_set_p0( UNUR_PAR *parameters, double p0 );
/* 
   Set design point for algorithm. It is used to split the domain of the 
   distribution. Values for @var{p0} close to the expectation of the
   distribution results in a relatively good performance of the algorithm.
   It is important that the hazard rate at this point must be greater
   than @code{0} and less than @code{UNUR_INFINITY}.
   
   Default: left boundary of domain + @code{1.}
*/

int unur_hri_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_hri_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the hazard rate is not bounded by the given bound, then
   @code{unur_errno} is set to @code{UNUR_ERR_GEN_CONDITION}. 

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/
