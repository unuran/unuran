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

   =REQUIRED PDF, all local extrema, cut-off values for the tails

   =OPTIONAL approximate area

   =SPEED Set-up: (very) slow, Sampling: fast

   =REF  [AJa93] [AJa95] [Book,Cha5.1]  

   =ABSTRACT Large tables necessary for heavy tailed distributions 

   =DESCRIPTION
      TABL (called Ahrens method in [Book]) is an acceptance/rejection
      method (see [ARVRej]) that uses a decomposition of the domain of
      the distribution into many short subintervals. Inside of these
      subintervals constant hat and squeeze functions are
      utilized. Thus it is easy to use the idea of immediate
      acceptance (see [ARVRej]) for points below the squeeze. This
      reduces the expected number of uniform random numbers per
      generated random variate to less than two. Using a large number
      of subintervals only little more than one random number is
      necessary on average. Thus this method becomes very fast.

      Due to the constant hat function this method only works for
      distributions with bounded domains. Thus for unbounded domains
      the left and right tails have to be cut off.  This is no problem
      when the probability of falling into these tail regions is
      beyond computational relevance (e.g. smaller than @code{1.e-12}).

      For easy construction of hat and squeeze functions it is necessary
      to know the regions of monotonicity (called @emph{slopes}) or
      equivalently all local maxima and minima of the density.  
      The main problem for this method in the setup is the choice of the
      subintervals. A simple and close to optimal approach is the 
      "equal area rule" [Book,Cha5.1]. There the subintervals are
      selected such that the area below the hat is the same for
      each subinterval which can be realized with a simple recursion.
      If more subintervals are necessary it is possible to split
      either randomly chosen intervals (adaptive rejection sampling, ARS)
      or those intervals, where the ratio between squeeze and hat is
      smallest. This version of the setup is called derandomized ARS
      (DARS). With the default settings TABL is first calculating
      approximately 30 subintervals with the equal area rule. Then
      DARS is used till the desired fit of the hat is reached.

      A convenient measure to control the quality of the fit of hat
      and squeeze is the ratio (area below squeeze)/(area below hat)
      called "sqhratio" which must be smaller or equal to one.
      The expected number of iterations in the rejection algorithm
      is known to be smaller than 1/sqhratio and the expected number
      of evaluations of the density is bounded by 1/sqhratio - 1.
      So values of the sqhratio close to one (e.g. @code{0.95} or
      @code{0.99}) lead to many subintervals. Thus a better fitting
      hat is constructed and the sampling algorithm becomes fast; on
      the other hand large tables are needed and the setup is very
      slow. For moderate values of sqhratio (e.g. @code{0.9} or
      @code{0.8}) the sampling is slower but the required tables are 
      smaller and the setup is not so slow.
      
      It follows from the above explanations that TABL is always
      requiring a slow setup and that it is not very well suited for
      heavy-tailed distributions. 

   
   HOWTOUSE

      For using the TABL method UNURAN needs a bounded interval to
      which the generated variates can be restricted and information
      about all local extrema of the distribution. For unimodal
      densities is is sufficient to provide the mode of the
      distribution. For the case of a built-in unimodal distribution
      with bounded domain all these information is present in the
      distribution object and thus no extra input is necessary (see
      example_TABL1 below).

      For a built-in unimodal distribution with unbounded domain we
      should specify the cut-off values for the tails. This can be
      done with the unur_tabl_set_boundary() call (see example_TABL2
      below). For the case that we do not set these boundaries the
      default values of @code{+/- 1.e20} are used. We can see in
      example_TABL1 that this still works fine for many standard
      distributions. 

      For the case of a multimodal distribution we have to set the
      regions of monotonicity (called slopes) explicitly using the
      unur_tabl_set_slopes() command (see example_TABL3 below).

      To controll the fit of the hat and the size of the tables and thus the
      speed of the setup and the sampling it is most convenient to use the
      unur_tabl_set_max_sqhratio() call. The default is 0.9 which is a sensible
      value for most distributions and applications. If very large samples of
      a distribution are required or the evaluation of a density is very slow
      it may be useful to increase the sqhratio to eg. 0.95 or even 0.99. With
      the unur_tabl_get_sqhratio() call we can check which sqhratio was really
      reached. If that value is below the desired value it is necessary to
      increase the maximal number of subintervals, which defaults to 1000,
      using the unur_tabl_set_max_intervals() call. 
      The unur_tabl_get_intervals() call can be used to find out the 
      number of subintervals the setup calculated.

      The usage of the commands mentioned here are demonstrated in
      example_TABL1, example_TABL2 and example_TABL3 below.

   =END

  include the C-code of the files example_TABL1, example_TABL2 and example_TABL3 here.
  TODO: Soll man example_TABL3 von stringinterface auf das andere GUI umstellen? 
  TODO: Soll man die Reihenfolge der function references aendern??

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
   sampling'' (DARS) is used in the setup.
   Intervals, where the area between hat and squeeze is too
   large compared to the average area between hat and squeeze
   over all intervals, are split.
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
   is set to @code{UNUR_INFINITY}.

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
   points are inserted via DARS in the setup.

   For the case of ARS (set_usedars() must be set to FALSE):
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.
   If @var{max_ratio} is close to one, many construction points are used.

   Default is @code{0.9}.
*/

double unur_tabl_get_sqhratio( const UNUR_GEN *generator );
/* 
   Get the current ratio (area below squeeze) / (area below hat)
   for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_tabl_get_hatarea( const UNUR_GEN *generator );
/* 
   Get the area below the hat for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_tabl_get_squeezearea( const UNUR_GEN *generator );
/* 
   Get the area below the squeeze for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added in or after the setup when the
   number of intervals suceeds @var{max_ivs}.

   Default is @code{1000}.
*/

int unur_tabl_get_intervals( const UNUR_GEN *generator );
/* 
   Get the current number of intervals.
   (In case of an error 0 is returned.)
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
   The list of slopes is given by an array @var{slopes} where each
   consecutive tuple (i.e. @code{(slopes[0], slopes[1])},
   @code{(slopes[2], slopes[3])}, etc.) defines one slope.
   Slopes must be sorted (i.e. both @code{slopes[0]} and
   @code{slopes[1]} must not be greater than any entry of the slope
   @code{(slopes[2], slopes[3])}, etc.)
   and must not be overlapping. Otherwise no slopes are set and
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

   Default is @code{-1.e20,1.e20}.
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


