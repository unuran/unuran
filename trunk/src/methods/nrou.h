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

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: slow or fast, Sampling: moderate

   =REF  [HLD04: Sect.2.4]

   =DESCRIPTION
      NROU is an implementation of the ratio-of-uniforms method
      which uses (minimal) bounding rectangles, see
      @ref{Ratio-of-Uniforms}. The coordinates of this rectangles are
      given by 

      @display
      @unurmath{v^+ = \sup\limits_{x}       \sqrt{f(x)},} @*
      @unurmath{u^- = \inf\limits_{x} (x-c) \sqrt{f(x)},} @*
      @unurmath{u^+ = \sup\limits_{x} (x-c) \sqrt{f(x)}.}
      @end display

      where @i{c} is the center of the distribution.
      These bounds can either be given directly, or these are computed
      automatically by means of an numerical routine.
      Of course this can fail, especially when this rectangle is not
      bounded. 

      It is important to note that the algorithm works with 
      @unurmath{PDF(x-center)} instead of 
      @unurmath{PDF(x)}, i.e. the bounding rectangle that have to be
      provided are for the @unurmath{PDF(x-center)}.
      This is important as otherwise the acceptance region can become
      a very long and skinny ellipsoid along a diagonal of the (huge)
      bounding rectangle.

   =HOWTOUSE
      For using the NROU method UNURAN needs the PDF of the
      distribution. 
      The bounding rectangle can be given by the unur_vnrou_set_u()
      and unur_vnrou_set_v() calls. If these are not called then the
      minimal bounding rectangle is computed automatically. Using 
      unur_vnrou_set_verify() and unur_vnrou_chg_verify() one can run
      the sampling algorithm in a checking mode, i.e., in every cycle
      of the rejection loop it is checked whether the used 
      rectangle indeed enclosed the acceptance region of the
      distribution. When in doubt (e.g., when it is not clear whether
      the numerical routine has worked correctly) this can be used to
      run a small Monte Carlo study.

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

int unur_nrou_set_u( UNUR_PAR *parameters, double umin, double umax );
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

int unur_nrou_set_v( UNUR_PAR *parameters, double vmax );
/* 
   Set upper boundary for bounding rectangle. If this value is not
   given then @unurmath{\sqrt{PDF(mode)}} is used instead.

   @emph{Notice}: When the mode is not given for the distribution
   object, then it will be computed numerically.

   Default: not set.
*/

int unur_nrou_set_r( UNUR_PAR *parameters, double r );
/*
   Sets the parameter @var{r} of the generalized ratio-of-uniforms 
   method.
   
   @emph{Notice}: This parameter must satisfy @var{r}>0.
   Setting to a nonpositive value is ignored and in this case the
   default value value is used instead.

   Default: @code{1}.
*/

int unur_nrou_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (@unurmath{\mu}) of the PDF.
   For distributions like the gamma distribution with large shape
   parameters the acceptance region becomes a long inclined skinny
   oval with a large bounding rectangle and thus an extremely large
   rejection constant. Using the @var{center} shifts the mode of the
   distribution near the origin and thus makes the bounding box of the
   acception region smaller.

   Default: Mode if known, else @code{0}.
*/

int unur_nrou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_nrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
