/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hrd.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HRD                                *
 *         (Hazard Rate Decreasing)                                          *
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
   =METHOD  HRD   Hazard Rate Decreasing

   =UP  Methods_for_CONT

   =REQUIRED decreasing (non-increasing) hazard rate 

   =SPEED Set-up: fast, Sampling: slow

   =DESCRIPTION
      Generates random variate with given non-increasing hazard rate. 
      It is necessary that the distribution object contains this hazard rate.
      Decreasing hazard rate implies that the corresponding PDF of the
      distribution has heavier tails than the exponential distribution
      (which has constant hazard rate).

      It is important to note that the domain of the distribution can
      be set via a unur_distr_cont_set_domain() call. However, only
      the left hand boundary is used. For computational reasons the
      right hand boundary is always reset to UNUR_INFINITY.
      If no domain is given by the user then the left hand boundary is
      set to @code{0}.
      
      For distributions which do not have decreasing hazard rates but
      are bounded from above use method HRB. 
      For distributions with increasing hazard rate method HRI is
      required.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hrd_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hrd_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_hrd_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the hazard rate is not bounded by the given bound, then
   @code{unur_errno} is set to @code{UNUR_ERR_GEN_CONDITION}. 

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/
