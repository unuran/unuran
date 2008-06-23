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
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
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

   =REQUIRED CDF

   =OPTIONAL PDF?

   =REF [HLD04: TODO]

   =SPEED Set-up: (very) slow, Sampling: (very) fast

   =REINIT not implemented

   =DESCRIPTION
      TODO!!

      PINV is a variant of numerical inversion, where the inverse CDF
      is approximated using Hermite interpolation, i.e., the interval 
      [0,1] is split into several intervals and in each interval the
      inverse CDF is approximated by polynomials constructed by means
      of values of the CDF and PDF at interval boundaries. This makes
      it possible to improve the accuracy by splitting a particular
      interval without recomputations in unaffected intervals. Three
      types of splines are implemented: linear, cubic, and quintic
      interpolation. For linear interpolation only the CDF is
      required. Cubic interpolation also requires PDF and quintic
      interpolation PDF and its derivative. 

      These splines have to be computed in a setup step. However, it
      only works for distributions with bounded domain; for
      distributions with unbounded domain the tails are chopped off
      such that the probability for the tail regions is small compared
      to the given u-resolution. 

      The method is not exact, as it only produces random variates of
      the approximated distribution. Nevertheless, the maximal
      numerical error in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be set to the
      required resolution (within machine precision).  
      Notice that very small values of the u-resolution are possible
      but may increase the cost for the setup step.

      As the possible maximal error is only estimated in the setup it
      may be necessary to set some special design points for computing
      the Hermite interpolation to guarantee that the maximal u-error
      can not be bigger than desired. Such points are points where the
      density is not differentiable or has a local extremum. Notice
      that there is no necessity to do so. However, if you do not
      provide these points to the algorithm there might be a small
      chance that the approximation error is larger than the given
      u-resolution, or that the required number of intervals is larger
      than necessary.

   =HOWTOUSE
      PINV works for continuous univariate distribution objects with
      given CDF and (optional) PDF. It uses Hermite interpolation of
      order 1, 3 [default] or 5. The order can be set by means of
      unur_pinv_set_order().

      For distributions with unbounded domains the tails are chopped
      off such that the probability for the tail regions is small
      compared to the given u-resulution. For finding these cut points
      the algorithm starts with the region @code{[-1.e20,1.e20]}. For
      the exceptional case where this might be too small (or one knows
      this region and wants to avoid this search heuristics) it can be
      directly set via a unur_pinv_set_boundary() call.

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be controlled
      by means of unur_pinv_set_u_resolution().

      As already mentioned the maximal error of this approximation is 
      only estimated. If this error is crucial for an application we
      recommend to compute this error using unur_pinv_estimate_error()
      which runs a small Monte Carlo simulation.

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
   Set order of Hermite interpolation. Valid orders are
   @code{1}, @code{3}, and @code{5}.
   Notice that @var{order} greater than @code{1} requires the density 
   of the distribution, and @var{order} greater than @code{3} even
   requires the derivative of the density. Using @var{order} @code{1}
   results for most distributions in a huge number of intervals
   and is therefore not recommended. If the maximal error in
   u-direction is very small (say smaller than @code{1.e-10}),
   @var{order} @code{5} is recommended as it leads to considerably 
   fewer design points, as long there are no poles or heavy tails.

   @emph{Remark:} When the target distribution has poles or (very) heavy
   tails @var{order} @code{5} (i.e., quintic interpolation) is 
   numerically less stable and more sensitive to round-off errors than
   @var{order} @code{3} (i.e., cubic interpolation).

   Default is @code{3} if the density is given and @code{1} otherwise.
*/

int unur_pinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal error in u-direction. However, the given u-error must not
   be smaller than machine epsilon (@code{DBL_EPSILON}) and should not be
   too close to this value. As the resoultion of most uniform random
   number sources is 2^(-32) = @code{2.3e-10}, a value of @code{1.e-10}
   leads to an inversion algorithm that could be called exact. For most
   simulations slighly bigger values for the maximal error are enough
   as well. 

   Default is @code{1.e-10}.
*/

int unur_pinv_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set the left and right boundary of the computational interval.
   The interval must cover the essential part of the distribution.
   When unur_pinv_set_searchboundary() is called with TRUE then
   this given domain is shortened to a domain of "computational
   relevance" such that the tail probabilities are smaller than given
   by unur_pinv_set_u_resolution().
   Thus it is usually safe to use a large interval.
   However, @code{+/- UNUR_INFINITY} is not allowed.

   @emph{Important}: This call does not change the domain of the
   given distribution itself. But it restricts the domain for the
   resulting random variates.

   Default is @code{1.e100}.
*/

int unur_pinv_set_searchboundary( UNUR_PAR *parameters, int left, int right );
/* 
   If @var{left} or @var{right} is set to FALSE then the respective
   boundary is used as given by a unur_pinv_set_boundary() call.
   However, these boundary points might cause numerical problems
   during the setup when PDF returns @code{0.} "almost everywhere".
   If set to TRUE (the default) then the computational interval is
   shortened to a more sensible region by means of a search algorithm.
   Switching off this search is useful, e.g. for the Gamma(2)
   distribution where the left border 0 is fixed and finite.

   @emph{Remark:}
   The searching algorithm assumes that the support of the distribution
   is connected.

   Default is TRUE.
*/

int unur_pinv_get_n_intervals( const UNUR_GEN *generator ); 
/* 
   Get number of nodes (design points) used for Hermite interpolation in 
   the generator object. The number of intervals is the number of
   nodes minus 1.
   It returns an error code in case of an error.
*/

double unur_pinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
/*
   Evaluate Hermite interpolation of inverse CDF at @var{u}.
   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
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
