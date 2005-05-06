/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ball.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method BALL                               *
 *         Ball-sampler from area below PDF or inside the RoU-shape          *
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
   =METHOD  BALL   ball-sampler

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: fast, Sampling: fast

   =DESCRIPTION
      BALL is an implementation of a ball-sampler.
      
      It is possible to sample directly from the area below the pdf
      or from the area inside the Ratio of Uniforms shape as chosen
      by the appropriate set_variant calls.
      
      The variant to sample directly below the pdf is currently
      not implemented.
      
      In case of the ratio-of-uniform variant an additional parameter 
      @i{r} that can be used for adjusting the algorithm to the 
      given distribution to improve performance and/or to make 
      this method applicable.  
      Moreover, this implementation uses the center @unurmath{\mu} 
      of the distribution (which is set to the mode or mean by default, 
      see unur_distr_cvec_get_center() for details of
      its default values).

   =HOWTOUSE
      For using the BALL method UNURAN needs the PDF of the
      distribution. Additionally, the parameter @i{r} can be set via
      a unur_vnrou_set_r() call.

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

UNUR_PAR *unur_ball_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/


int unur_ball_set_variant_pdf( UNUR_PAR *par );
/* 
   Sampler variant :
   Sampling from the area below pdf. (not yet implemented)
*/


int unur_ball_set_variant_rou( UNUR_PAR *par );
/* 
   Sampler variant :
   Sampling from the area inside the RoU shape.
   
   This is the default.
*/

int unur_ball_set_r( UNUR_PAR *parameters, double r );
/*
   Sets the parameter @var{r} of the generalized multivariate
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{r}>0.

   Default: @code{1}.
*/

int unur_ball_set_ball_radius( struct unur_par *par, double ball_radius );
/* 
   Sets (initial) radius of ball used for the ball-sampler. 

   Default: @code{1.0}
*/


int unur_ball_set_adaptive_ball( UNUR_PAR *parameters, int adaptive_flag );
/*
   Increasing/Decreasing ball radius whenever candidate point is inside/outside 
   the RoU shape.
   
   Default: @code{0}.
*/

int unur_ball_set_skip( UNUR_PAR *parameters, long skip );
/*
   Sets the parameter @var{skip} i.e. the number steps 
   between two points that will be used as random numbers.
   
   @emph{Notice}: This parameter must satisfy @var{skip}>=0.

   Default: @code{0}.
*/


/* =END */


/*---------------------------------------------------------------------------*/


void _unur_ball_set_point( UNUR_GEN *gen, double *uv);
/* set the current point (dimension=dim+1) inside the RoU-shape */

void _unur_ball_get_point( UNUR_GEN *gen, double *uv);
/* get the current point (dimension=dim+1) inside the RoU-shape */

