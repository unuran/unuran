/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tabl.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TABL                               *
 *         (Ahren's TABLe method: piecewise constant hat)                    *
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
   =METHOD  TABL   a TABLe method with piecewise constant hats

   =UP  Methods_for_CONT

   =REQUIRED PDF, all local extrema

   =OPTIONAL approximate area

   =SPEED Set-up: slow, Sampling: fast

   =REF  [AJa93]  [AJa95]

   =ABSTRACT Not suited for heavy tailed distributions

   =DESCRIPTION
      TABL is an acceptance/rejection method that uses piecewise
      constant hat and squeezes. Immediate acceptance of points below the
      squeeze reduces the expected number of uniform random numbers to
      less than two and makes this method extremely fast.
      
      The method only works for distributions with bounded domain. Thus
      for unbounded domains the left and right tails are cut off
      (the cutting points can be set by a unur_tabl_set_boundary() call).
      This is no problem when the probability of falling into these tail
      regions is beyond computational relevance.
      
      The method works for all probability density functions where the
      regions of monotonicity (called slopes) are given. This can be done
      explicitly by a unur_tabl_set_slopes() call. If (and only if) no
      slopes are given, the domain and the mode of the PDF are used to
      compute the slopes. If neither slopes nor the mode and the domain
      are given initializing of the generator fails.
      
      In the setup first the equal area rule is used to construct a hat
      function, i.e., the interval boundaries are chosen such that the
      area below each interval is equal to a given fraction of the total
      area below the given PDF. This fraction can be set by a
      unur_tabl_set_areafraction() call. Additionally these intervals are
      split until the maximum number of intervals is reached or the
      ratio between the area below squeeze and the area below the hat is
      exceeded. 
      
      It is possible to switch off this second setup step (called 
      derandomized adaptive rejection sampling -- DARS). Then adaptive
      rejection sampling is used to split these intervals. There are
      three variants for adaptive rejection sampling. These differ in the
      way how an interval is split: 
      @enumerate
      @item use the generated point to split the interval;
      @item use the mean point of the interval; or
      @item use the arcmean point.
      @end enumerate

      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_tabl_set_verify() and unur_tabl_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tabl_new( const UNUR_DISTR* distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_tabl_set_usedars( UNUR_PAR *parameters, int usedars );
/*
   If @var{usedars} is set to TRUE, ``derandomized adaptive rejection
   sampling'' (DARS) is used in setup.
   Intervals, where the area between hat and squeeze is too
   large compared to the average area between hat and squeeze
   over all intervals, are splitted.
   This procedure is repeated until the ratio between squeeze and hat
   exceeds the bound given by unur_tabl_set_max_sqhratio() call or the
   maximum number of intervals is reached. Moreover, it also aborts
   when no more intervals can be found for splitting.

   For finding splitting points the arc-mean rule (a mixture of
   arithmetic mean and harmonic mean) is used.

   Default is TRUE.
*/

int unur_tabl_set_darsfactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for ``derandomized adaptive rejection sampling''.
   This factor is used to determine the segments that are ``too
   large'', that is, all segments where the area between squeeze and
   hat is larger than @var{factor} times the average area over all
   intervals between squeeze and hat.
   Notice that all segments are split when @var{factor} is set to
   @code{0.}, and that there is no splitting at all when @var{factor}
   is set to UNUR_INFINITY.

   Default is @code{0.99}. There is no need to change this parameter.
*/

int unur_tabl_set_variant_splitmode( UNUR_PAR *parameters, unsigned splitmode );
/* 
   There are three variants for adaptive rejection sampling. These
   differ in the way how an interval is split:
   @table @r
   @item splitmode @code{1}
   use the generated point to split the interval.
   @item splitmode @code{2}
   use the mean point of the interval.
   @item splitmode @code{3}
   use the arcmean point;
   suggested for distributions with heavy tails.
   
   @end table
   Default is splitmode @code{2}.
*/

int unur_tabl_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area below squeeze) / (area below hat).
   It must be a number between 0 and 1.
   When the ratio exceeds the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.
   If @var{max_ratio} is close to one, many construction points are used.

   Default is @code{0.9}.
*/

double unur_tabl_get_sqhratio( const UNUR_GEN *generator );
/* 
   Get the current ratio (area below squeeze) / (area below hat)
   for the generator. (In case of an error @code{0} is returned.)
*/

double unur_tabl_get_hatarea( const UNUR_GEN *generator );
/* 
   Get the area below the hat for the generator.
   (In case of an error @code{0} is returned.)
*/

double unur_tabl_get_squeezearea( const UNUR_GEN *generator );
/* 
   Get the area below the squeeze for the generator.
   (In case of an error @code{0} is returned.)
*/

int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added after the setup when the number of
   intervals suceeds @var{max_ivs}.

   Default is @code{1000}.
*/

int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
/* 
   Set parameter for equal area rule. During the setup a piecewise
   constant hat is constructed, such that the area below each of these
   pieces (strips) is the same and equal to the (given) area below the
   distribution times @var{fraction} (which must be greater than
   zero).

   @emph{Important:} It the area below the PDF is not set, then 1 is
   assumed. 

   Default is @code{0.1}.
*/

int unur_tabl_set_nstp( UNUR_PAR *parameters, int n_stp );
/* 
   Set number of construction points for the hat function. @var{n_stp}
   must be greater than zero. After the setup there are about
   @var{n_stp} construction points. However it might be larger when a
   small fraction is given by the unur_tabl_set_areafraction() call.
   It also might be smaller for some variants.

   Default is @code{30}.
*/

int unur_tabl_set_slopes( UNUR_PAR *parameters, const double *slopes, int n_slopes );
/* 
   Set slopes for the PDF.
   A slope <a,b> is an interval [a,b] or [b,a] where the PDF is
   monotone and PDF(a) >= PDF(b). 
   The list of slopes are given by an array @var{slopes} where each
   consecutive tuples (i.e. @code{(slopes[0], slopes[1])}, 
   @code{(slopes[2], slopes[3])}, etc.) is one slopes.
   Slopes must be sorted (i.e. both @code{slopes[0]} and
   @code{slopes[1]} must not be greater than any entry of the slope
   @code{(slopes[2], slopes[3])}, etc.)
   and must not overlapping. Otherwise no slopes are set and
   @var{unur_errno} is set to @code{UNUR_ERR_PAR_SET}.

   @emph{Notice:} @var{n_slopes} is the number of slopes (and not the
   length of the array @var{slopes}).

   @emph{Notice} that setting slopes resets the given domain for the
   distribution. However in case of a standard distribution the area
   below the PDF is not updated.
*/

int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is @code{1}.
*/

int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set the left and right boundary of the computation interval.
   The piecewise hat is only constructed inside this interval. The
   probability outside of this region must not be of
   computational relevance.
   Of course @code{+/- UNUR_INFINITY} is not allowed.

   Default is @code{1.e20}.
*/

int unur_tabl_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_tabl_chg_verify( UNUR_GEN *generator, int verify );
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


