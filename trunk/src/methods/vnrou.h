/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vnrou.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VNROU                              *
 *         (Vector Naive Ratio-Of-Uniforms method)                           *
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
   =METHOD  VNROU   Naive Ratio-Of-Uniforms method

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, bounding rectangle for acceptance region

   =SPEED Set-up: fast, Sampling: slow

   =REF  [HLD04: Sect.2.4]

   =DESCRIPTION
      VNROU is an implementation of the ratio-of-uniforms method
      which uses a (minimal) bounding hyper-rectangle.

      See the Appendix for more details on the ratio-of-uniforms method.
   =END
*/


/* 
   TODO:

   Set the center (@unurmath{\mu}) of the PDF.
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
   used) a @code{0}-vector is used as center.
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vnrou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_vnrou_set_u( UNUR_PAR *parameters, double *umin, double *umax );
/* 
   Sets left and right boundaries of bounding hyper-rectangle.
   If no values are given, the boundary of the minimal bounding
   hyper-rectangle is computed numerically.
  
   @strong{Important}: The boundaries are those of the density shifted
   by the center of the distribution.

   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. In particular for multimodal
   distributions this might fail.

   Default: not set (i.e. computed automatically)
*/

int unur_vnrou_set_v( UNUR_PAR *parameters, double vmax );
/* 
   Set upper boundary for bounding hyper-rectangle. 
   If no values are given, the density at the mode is evaluated.
   If no mode is given for the distribution it is computed
   numercally (and might fail).
  
   Default: not set (i.e. computed automatically)
*/

int unur_vnrou_set_r( UNUR_PAR *parameters, double r );
/* 
   Sets the parameter @var{r} of the generalized multivariate 
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{r}>0. 
   Setting to a nonpositive value is ignored and in this case the
   default value value is used instead.

   Default: @code{1}.
*/

int unur_vnrou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_vnrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
