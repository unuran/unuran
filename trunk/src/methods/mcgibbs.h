/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mcgibbs.h                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MCGIBBS                            *
 *         (Markov Chain - GIBBS sampler)                                    *
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
   =METHOD  MCGIBBS   Markov Chain - GIBBS sampler

   =UP  Methods_for_CVEC

   =REQUIRED T-concave logPDF, derivatives of logPDF

   =SPEED Set-up: fast, Sampling: moderate

   =REF  [HLD04: Sect.14.1.2]

   =DESCRIPTION
      Methods MCGIBBS implements a Gibbs sampler for multivariate
      distributions with given joint density and its gradient.
      When running such a Markov chain all coordinates are updated
      cyclically using full conditional distributions. Each step
      is returned (i.e., a random point is returned whenever a single
      coordinate has been updated). 
      It is also possible to return only points after all coordinates
      have been updated by "thinning" the chain.
      Moreover, to reduce autocorrelation this thinning factor can be
      any integer. Notice, however, that the sampling time for a chain
      of given length is also increased by the same factor.

      MCGIBBS also provides a variant of the Gibbs sampler where in
      each step a point from the full conditional distribution along
      some random direction is sampled. This direction is chosen
      uniformly from the sphere in each step.
      This method is also known as Hit-and-Run algorithm for
      non-uniform distributions.

      Our experiences shows that the original Gibbs sampler with
      sampling along coordinate axes is superior to random direction
      sampling as long as the correlations between the components of
      the random vector are not too high.

      For both variants TDR ((@pxref{TDR}, @pxref{TDRGW}) is used to
      sample from the full conditional distributions. In opposition to
      the univariate case, it is important that the factor @code{c} is
      as large as possible. I.e., for a log-concave density @code{c}
      must be set to @code{0.}, since otherwise numerical underflow
      might stop the algorithm.

      @emph{Important:} MCGIBBS does not generate independent random
      points.

   =HOWTOUSE
      For using the MCGIBBS method UNURAN needs the logarithm of the
      PDF of the multivariate joint distribution and its gradient or
      partial derivatives. 

      It provides two variants:
      @table @emph 
      @item coordinate direction sampling (Gibbs sampling) [default]
      The coordinates are updated cyclically.
      It requires the partial derivatives of the (logarithm of the)
      PDF of the target distribution, 
      see unur_distr_cvec_set_pdlogpdf(). 
      Otherwise, the gradient of the logPDF 
      (see unur_distr_cvec_set_dlogpdf()) 
      is used, which is more expensive. 

      This variant can be selected using
      unur_mcgibbs_set_variant_coordinate().

      @item random direction sampling (nonuniform Hit-and-Run algorithm)
      In each step is a direction is sampled uniformly from the sphere
      and the next point in the chain is sampled from the full
      conditional distribution along this direction.

      It requires the gradient of the logPDF and thus each step is
      more expensive than each step for coordinate direction sampling.

      This variant can be selected using
      unur_mcgibbs_set_variant_random_direction().
      @end table

      It is important that the @code{c} parameter for the TDR method
      is as large as possible. For logconcave distribution it must be
      set to @code{0.}, since otherwise numerical underflow can cause
      the algorithm to stop.
      
      In case of a fatal error in the generator for conditional
      distributions the methods generates points that contain
      UNUR_INFINITY.

      @strong{Warning:} The algorithm requires that all full
      conditionals for the given distribution object are
      @i{T}-concave. However, this property is not checked. 
      If this property is not satisfied, then generation from the
      conditional distributions becomes (very) slow and might fail or
      (even worse) produces random vectors from an incorrect
      distribution. 

      @strong{Warning:} Be carefull with debugging flags. If it
      contains flag @code{0x01000000u} it produces a lot of output for
      each step in the algorithm.
      (This flag is switched of in the default debugging flags).

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mcgibbs_new( const UNUR_DISTR *distribution );

/*...........................................................................*/

int unur_mcgibbs_set_variant_coordinate( UNUR_PAR *parameters );
/* 
   Coordinate Direction Sampling:
   Sampling along the coordinate directions (cyclic).

   This is the default.
*/

int unur_mcgibbs_set_variant_random_direction( UNUR_PAR *parameters );
/* 
   Random Direction Sampling:
   Sampling along the random directions.
*/

int unur_mcgibbs_set_c( UNUR_PAR *parameters, double c );
/* 
   Set parameter @var{c} for transformation @unurmath{T} of the
   transformed density rejection method.
   Currently only values between @code{0} and @code{-0.5} are
   allowed. If @code{c} is between @code{0} and @code{-0.5} it is set
   to @code{-0.5}. 

   For @var{c} @code{=0} (for logconcave densities) method TDRGW
   (@pxref{TDRGW}) is used which is very robust against badly
   normalized PDFs. For other values method TDR (@pxref{TDR}) is used.

   The value for @var{c} should be as large as possible to avoid 
   fatal numerical underflows. Thus for log-concave distributions
   @var{c} must be set to @code{0.}
 
   Default is @code{0}.
*/


int unur_mcgibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0);
/* 
   Sets the starting point of the Gibbs sampler. @var{x0} must be 
   a "typical" point of the given distribution.

   Default is 0.
*/

int unur_mcgibbs_set_thinning( UNUR_PAR *parameters, int thinning );
/*
   Sets the @var{thinning} parameter. When @var{thinning} is set to
   @i{k} then every @i{k}-th point from the iteration is returned by
   the sampling algorithm.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=1.

   Default: @code{1}.
*/

/*...........................................................................*/

/* =END */
/*---------------------------------------------------------------------------*/


