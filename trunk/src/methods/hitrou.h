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

   =OPTIONAL mode, center, bounding rectangle

   =SPEED Set-up: fast, Sampling: fast 

   =REF  

   =DESCRIPTION 
      HITROU is a variation of the multivariate
      ratio-of-uniforms method which uses a hit-and-run random walk
      algorithm to sample points from the enclosed volume. see also 
      @ref{Ratio-of-Uniforms}.  
      

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

int unur_hitrou_set_skip( UNUR_PAR *parameters, long r );
/* 
   Sets the parameter @var{skip} of the generalized multivariate 
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{skip}>=0. 

   Default: @code{10}.
*/

int unur_hitrou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_hitrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
