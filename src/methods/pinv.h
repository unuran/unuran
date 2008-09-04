/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: pinv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method PINV                               *
 *         (Polynomial interpolation based INVersion of CDF)                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
   =METHOD PINV  Polynomial interpolation based INVersion of CDF

   =UP  Methods_for_CONT

   =REQUIRED PDF, center

   =OPTIONAL domain

   =REF [HLD04: FIXME]

   =SPEED Set-up: (very) slow, Sampling: (very) fast

   =REINIT not implemented

   =DESCRIPTION
      PINV is a variant of numerical inversion, where the inverse CDF
      is approximated using Newton's recursion for interpolating
      polynomials. The interval [0,1] is split into several
      subintervals and in each subinterval the inverse CDF is
      approximated by polynomials constructed by means of values of
      the CDF. The CDF is numerically computed from the given PDF
      using adaptive Gauss-Lobatto integration with 5 points.
      Subintervals are splitted until the requested accuracy goal is
      reached.

      The interpolating polynomials have to be computed in a setup
      step. However, it only works for distributions with bounded
      domain; for distributions with unbounded domain the tails are
      chopped off such that the probability for the tail regions is
      small compared to the given u-resolution. 

      The method is not exact, as it only produces random variates of
      the approximated distribution. Nevertheless, the maximal
      numerical error in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be set to the
      required resolution (within machine precision).  
      Notice that very small values of the u-resolution are possible
      but may increase the cost for the setup step.

      The construction of the interpolation polynomial only works when
      the PDF is unimodal or when the PDF does not vanish between two
      modes. 


   =HOWTOUSE
      PINV works for continuous univariate distribution objects with
      given PDF. It uses Gauss-Lobatto integration with 5 points to
      estimate the CDF and Newton's recursion for the interpolating
      polynomial. The order of the polynomial can be set by means of
      the unur_pinv_set_order() call. 

      For distributions with unbounded domains the tails are chopped
      off such that the probability for the tail regions is small
      compared to the given u-resulution. For finding these cut points
      the algorithm starts with the region @code{[-1.e100,1.e100]}. For
      the exceptional case where this does not work it these starting
      points can be changed via a unur_pinv_set_boundary() call.
      It is important that the given distribution object contains the
      mode or a typical point of the distribution. The latter can be
      set for distribution by means of a unur_distr_cont_set_center()
      call. (It neither the mode nor the center is set for a
      distribution object, then @code{0} is assumed.)

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be controlled
      by means of unur_pinv_set_u_resolution().

      The maximal error of this approximation is only estimated. If
      this error is crucial for an application we recommend to compute
      this error using unur_pinv_estimate_error() which runs a small
      Monte Carlo simulation.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_pinv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_pinv_set_order( UNUR_PAR *parameters, int order);
/* 
   Set order of interpolation. Valid orders are between @code{2} and
   @code{19}. Higher orders result in fewer intervals for the
   approximations. 

   Default: @code{5}.
*/

int unur_pinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal error in u-direction. However, the given u-error must not
   be smaller than machine epsilon (@code{DBL_EPSILON}) and should not be
   too close to this value. As the resolution of most uniform random
   number sources is 2^(-32) = @code{2.3e-10}, a value of @code{1.e-10}
   leads to an inversion algorithm that could be called exact. For most
   simulations slightly bigger values for the maximal error are enough
   as well. 

   Default is @code{1.e-10}.
*/

int unur_pinv_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set @var{left} and @var{right} point for finding the cut-off points
   for the "computational domain", i.e., the domain that covers the
   essential part of the distribution.
   The cut-off points are computed such that the tail probabilities
   are smaller than given by unur_pinv_set_u_resolution().
   It is usually safe to use a large interval.
   However, @code{+/- UNUR_INFINITY} is not allowed.

   @emph{Important}: This call does not change the domain of the
   given distribution itself. But it restricts the domain for the
   resulting random variates.

   Default: intersection of @code{[-1.e100,+1.e100]} and the given
   domain of the distribution.
*/

int unur_pinv_set_searchboundary( UNUR_PAR *parameters, int left, int right );
/* 
   If @var{left} or @var{right} is set to FALSE then the respective
   boundary as given by a unur_pinv_set_boundary() call is used
   without any further computations.
   However, these boundary points might cause numerical problems
   during the setup when PDF returns @code{0.} ``almost everywhere''.
   If set to TRUE (the default) then the computational interval is
   shortened to a more sensible region by means of a search algorithm.
   Switching off this search is useful, e.g. for the Gamma(2)
   distribution where the left border 0 is fixed and finite.

   @emph{Remark:}
   The searching algorithm assumes that the support of the distribution
   is connected.

   @emph{Remark:}
   Do not set this parameter to FALSE except when searching for
   cut-off points fails and one wants to try with precomputed values.

   Default: TRUE.
*/

int unur_pinv_get_n_intervals( const UNUR_GEN *generator ); 
/* 
   Get number of intervals used for interpolation in 
   the generator object.
   It returns @code{0} in case of an error.
*/

double unur_pinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
/*
   Evaluate interpolation of inverse CDF at @var{u}.
   If @var{u} is out of the domain (0,1) then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).
*/

int unur_pinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
/*
   Estimate maximal u-error and mean absolute error (MAE) for @var{generator}
   by means of Monte-Carlo simulation with sample size @var{samplesize}.
   The results are stored in @var{max_error} and @var{MAE}, respectively.

   It returns @code{UNUR_SUCCESS} if successful. 
*/

/* =END */
/*---------------------------------------------------------------------------*/
