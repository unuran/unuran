/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: itdr.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method ITDR                               *
 *         (Inverse Transformed Density Rejection)                           *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
     $Id: nrou.h 2124 2004-12-02 14:56:49Z leydold $
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
   =EXPMETHOD  ITDR   Inverse Transformed Density Rejection

   =UP  Methods_for_CONT

   =REQUIRED monotone PDF, pole

   =OPTIONAL splitting point between pole and tail region, c-values

   =SPEED Set-up: moderate, Sampling: moderate

   =REF  [HLD04: Sect.??]

   =DESCRIPTION
      ITDR is ...

   =HOWTOUSE
      For using the ITDR method UNURAN needs the PDF of the
      distribution an the pole.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_itdr_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_itdr_set_bx( UNUR_PAR *parameters, double bx );
/* 
   Sets splitting point @var{bx} between pole and tail region.
   If no such point is provided it will be set automatically.

   Default: not set.
*/

int unur_itdr_set_cp( UNUR_PAR *parameters, double cp );
/* 
   Sets parameter @var{cp} for transformation T for inverse 
   density in pole region.
   It must be at most 0 and greater than -1. 
   If no @var{cp}-value is given it is estimated automatically.

   Default: not set.
*/

int unur_itdr_set_ct( UNUR_PAR *parameters, double ct );
/* 
   Sets parameter @var{ct} for transformation T for 
   density in tail region.
   It must be at most 0. For densities with unbounded domain
   it must be greater than -1. 
   If no @var{ct}-value is given it is estimated automatically.

   Default: not set.
*/

int unur_itdr_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition @unurmath{PDF(x) \leq hat(x)} is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However, notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_itdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
