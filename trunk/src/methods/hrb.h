/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hrb.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HRB                                *
 *         (Hazard Rate Bounded)                                             *
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
   =METHOD  HRB   Hazard Rate Bounded

   =UP  Methods_for_CONT

   =REQUIRED bounded hazard rate 

   =OPTIONAL upper bound for hazard rate

   =SPEED Set-up: fast, Sampling: slow

   =DESCRIPTION
      Generates random variate with given hazard rate. It requires that
      the distribution object contains a hazard rate and it requires an
      upper bound for the hazard rate which must be set using
      unur_hrb_set_upperbound() call. If no such upper bound is given
      it is assumed that the upper bound can be achieved by evaluating 
      the hazard rate at the left hand boundary of the domain of the 
      distribution.

      It is important to note that the domain of the distribution can
      be set via a unur_distr_cont_set_domain() call.
      However, the left border must not be negative. Otherwise it is
      set to @code{0}. This is also the default if no domain is
      given at all. For computational reasons the right border is
      always set to @code{UNUR_INFINITY} independently of the given 
      domain. Thus for domains bounded from right the function for
      computing the hazard rate should return @code{UNUR_INFINITY}
      right of this domain.

      For distributions with decreasing hazard rates use method 
      HRD, which is faster.
      For distributions with increasing hazard rate method HRI is
      required.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hrb_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hrb_set_upperbound( UNUR_PAR *parameters, double upperbound );
/* 
   Set upper bound for hazard rate. If this call is not used it is
   assumed that the the maximum of the hazard rate is achieved at the
   left hand boundary of the domain of the distribution.
*/

int unur_hrb_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_hrb_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the hazard rate is not bounded by the given bound, then
   @code{unur_errno} is set to @code{UNUR_ERR_GEN_CONDITION}. 

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/
