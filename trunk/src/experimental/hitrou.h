/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hitrou.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HITROU                             *
 *         (HIT and run Ratio-Of-Uniforms method)                            *
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
   =METHOD  HITROU   Multivariate HIT and run Ratio-Of-Uniforms method

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: fast, Sampling: fast

   =DESCRIPTION
      HITROU is an implementation of the multivariate
      ratio-of-uniforms method which uses a (minimal) bounding
      hyper-rectangle, see also @ref{Ratio-of-Uniforms}.  It uses an
      additional parameter @i{r} that can be used for adjusting the
      algorithm to the given distribution to improve performance
      and/or to make this method applicable.  Larger values of
      @i{r} increase the class of distributions for which the
      method works at the expense of higher rejection
      constants. Moreover, this implementation uses the center
      @unurmath{\mu} of the distribution (which is set to the mode or
      mean by default, see unur_distr_cvec_get_center() for details of
      its default values).

      The minimal bounding has then the coordinates

      @unurmathdisplay{
      v^+   = \sup\limits_{x}               (f(x))^{1/r\,d+1}, \\
      u^-_i = \inf\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, \\
      u^+_i = \sup\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, }

      where @unurmath{x_i} is the @i{i}-th coordinate of point @i{x};
      @unurmath{\mu_i} is the @i{i}-th coordinate of the center
      @unurmath{\mu.}
      @unurmath{d} denotes the dimension of the distribution.
      These bounds can either be given directly, or are computed
      automatically by means of an numerical routine
      by Hooke and Jeeves @unurbibref{HJa61} called direct search
      (see @file{src/utils/hooke.c} for further references and
      details). Of course this algorithm can fail, especially when
      this rectangle is not bounded.

      It is important to note that the algorithm works with
      @unurmath{PDF(x-center)} instead of
      @unurmath{PDF(x),} i.e. the bounding rectangle has to be
      provided for @unurmath{PDF(x-center).}
      This is important as otherwise the acceptance region can become
      a very long and skinny ellipsoid along a diagonal of the (huge)
      bounding rectangle.

   =HOWTOUSE
      For using the HITROU method UNURAN needs the PDF of the
      distribution. Additionally, the parameter @i{r} can be set via
      a unur_vnrou_set_r() call.

      A bounding rectangle can be given by the
      unur_vnrou_set_u() and unur_vnrou_set_v() calls.

      @emph{Important:} The bounding rectangle has to be
      provided for the function @unurmath{PDF(x-center)!}
      Notice that @code{center} is the center of the given
      distribution, see unur_distr_cvec_set_center().
      If in doubt or if this value is not optimal, it can be changed
      (overridden) by a unur_distr_cvec_set_center() call.

      If the coordinates of the bounding rectangle are not provided by
      the user then the minimal bounding rectangle is computed
      automatically.

      By means of unur_vnrou_set_verify() and unur_vnrou_chg_verify()
      one can run the sampling algorithm in a checking mode, i.e., in
      every cycle of the rejection loop it is checked whether the used
      rectangle indeed enclosed the acceptance region of the
      distribution. When in doubt (e.g., when it is not clear whether
      the numerical routine has worked correctly) this can be used to
      run a small Monte Carlo study.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hitrou_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hitrou_set_u( UNUR_PAR *parameters, double *umin, double *umax );
/*
   Sets left and right boundaries of bounding hyper-rectangle.
   If no values are given, the boundary of the minimal bounding
   hyper-rectangle is computed numerically.

   @strong{Important}: The boundaries are those of the density shifted
   by the center of the distribution, i.e., for the
   function @unurmath{PDF(x-center)!}

   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. Moreover, for multimodal distributions
   the bounds might be too small as only local extrema are computed.
   Nevertheless, for log-concave distributions it should work.

   Default: not set (i.e. computed automatically)
*/

int unur_hitrou_set_v( UNUR_PAR *parameters, double vmax );
/*
   Set upper boundary for bounding hyper-rectangle.
   If no values are given, the density at the mode is evaluated.
   If no mode is given for the distribution it is computed
   numerically (and might fail).

   Default: not set (i.e. computed automatically)
*/

int unur_hitrou_set_r( UNUR_PAR *parameters, double r );
/*
   Sets the parameter @var{r} of the generalized multivariate
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{r}>0.

   Default: @code{1}.
*/

int unur_hitrou_set_skip( UNUR_PAR *parameters, long skip );
/*
   Sets the parameter @var{skip} i.e. the number of hit and
   run steps between two points in the (U,V)-space that will
   be used as random numbers @unurmath{X=X_0+U/V^r}

   @emph{Notice}: This parameter must satisfy @var{skip}>=0.

   Default: @code{0}.
*/

int unur_hitrou_set_u_planes( UNUR_PAR *parameters, int u_planes );
/*
   Sets the @var{u_planes} flag i.e. if we should calculate
   and use all the u-planes of the bounding rectangle.
   In case, that this flag is set to 0, we will only calculate
   and use the planes perpendicular to the v-coordinate - hence
   in that case, our bounding shape is not a finite rectangle
   but an infinite strip.

   Default: @code{0}.
*/

int unur_hitrou_set_adaptive( UNUR_PAR *parameters, int adaptive_flag );
/*
   When calculating a random point along the given direction,
   this adaptive flag controls how to proceed, when the random
   point is located outside the shape.
   In case, that the adaptive_flag is set to 0, we sample again
   without changing the ends of the direction-line-segment.
   In case, that the adaptive_flag is set to 1, we sample again
   reusing the previous outside point as new end of the direction-
   line-segment.

   Default: @code{1}.
*/

int unur_hitrou_set_variant_coordinate( UNUR_PAR *par );
/* 
   Coordinate Sampler :
   Sampling along the coordinate directions (cyclic).
*/

/*...........................................................................*/

int unur_hitrou_set_variant_random_direction( UNUR_PAR *par );
/* 
   Random Direction Sampler :
   Sampling along the random directions.
   
   This is the default.
*/

/*...........................................................................*/



/* =END */
/*---------------------------------------------------------------------------*/

long _unur_hitrou_get_pdfcount( UNUR_GEN *gen);
/* Return the number of PDF calls */

void _unur_hitrou_reset_pdfcount( UNUR_GEN *gen);
/* Reset the number of PDF calls to 0 */

void _unur_hitrou_set_shape( UNUR_GEN *gen, int shape_flag);
/* shape_flag : 0=normal, 1=rectangle, 2=single simplex, 3=stacked simplex */

void _unur_hitrou_set_testrectangle( UNUR_GEN *gen, double *relative_size);
/* set the relative size of the test rectangle (relative to bounding rect) */

void _unur_hitrou_set_point( UNUR_GEN *gen, double *uv);
/* set the current point (dimension=dim+1) inside the testrectangle */

void _unur_hitrou_get_point( UNUR_GEN *gen, double *uv);
/* get the current point (dimension=dim+1) inside the testrectangle */

long _unur_hitrou_get_simplex_jumps( UNUR_GEN *gen);
/* Return the number of simplex jumps for double-simplex shape  */

void _unur_hitrou_reset_simplex_jumps( UNUR_GEN *gen);
/* Reset the number of simplex jumps to 0 */

