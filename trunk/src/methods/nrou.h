/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: nrou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NROU                               *
 *         (Naive Ratio-Of-Uniforms method)                                  *
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
   =METHOD  NROU   Naive Ratio-Of-Uniforms method

   =UP  Methods_for_CONT

   =REQUIRED PDF

   =OPTIONAL mode, bounding rectangle for acceptance region

   =SPEED Set-up: fast, Sampling: slow

   =REF  [HDL04: Sect.2.4]

   =DESCRIPTION
      NROU is an implementation of the ratio-of-uniforms method
      which uses (minimal) bounding rectangles.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_nrou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_nrou_set_rect_u( UNUR_PAR *parameters, double umin, double umax );
/* 
   Sets left and right boundary of bounding rectangle.
   If no values are given, the boundary of the minimal bounding
   rectangle is computed numerically.
   
   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. In particular for multimodal
   distributions this might fail.
   For @unurmath{T_c}-concave distributions with @unurmath{c=-1/2} it
   should work.

   Default: not set.
*/

int unur_nrou_set_rect_v( UNUR_PAR *parameters, double vmax );
/* 
   Set upper boundary for bounding rectangle. If this value is not
   given then @unurmath{\sqrt{PDF(mode)}} is used instead.

   @emph{Notice}: When the mode is not given for the distribution
   object, then it will be computed numerically.

   Default: not set.
*/

int unur_nrou_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (approximate mode) of the PDF.
   For distributions like the gamma distribution with large shape
   parameters the acceptance region becomes a long inclined skinny
   oval with a large bounding rectangle and thus an extremely large
   rejection constant. Using the @var{center} shifts the mode of the
   distribution near the origin and thus makes the bounding box of the
   acception region smaller.

   Default: Mode of the distribution if neither unur_nrou_set_rect_u()
   nor unur_nrou_set_rect_v() are called (and the mode is available
   for the distribution).
   Otherwise (if at least one of these two calls has been
   used @code{0} is used as center.
*/

int unur_nrou_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_nrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/
