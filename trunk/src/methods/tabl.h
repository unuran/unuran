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

   TABL is an acceptance/rejection method that uses piecewise
   constant hat and squeezes. Immediate acceptance of points below the
   squeeze reduces the expected number of uniform random numbers to
   less than two and makes this method extremely fast.

   The method only works for distributions with bounded domain. Thus
   for unbounded domains the left and right tails are cut off
   (the cutting points can be set by a unur_tabl_set_boundary() call).
   This no problem when the probability of falling into these tail
   regions is beyond computational relevance.

   The method works for all probability density functions where the
   regions of monotonicity (called slopes) are given. This can be done
   explicitly by a unur_tabl_set_slopes() call. If (and only if) no
   slopes are given, the domain and the mode of the p.d.f. are used to
   compute the slopes. If neither slopes nor the mode and the domain
   are given initializing of the generator fails.

   In the setup first the equal area rule is used to construct a hat
   function, i.e., the interval boundaries are chosen such that the
   area below each interval is equal to a given fraction of the total
   area below the given p.d.f. This fraction can be set by a
   unur_tabl_set_areafraction() call. Additionally these intervals are
   split by until the maximum number of intervals is reached or the
   ration between the area below squeeze and the area below the hat is
   exceeded. 

   It is possible to switch off this second setup step. Then adaptive
   rejection sampling is used to split these intervals. There are
   three variants for adaptive rejection sampling. These differ in the
   way how an interval is split: 
   (1) use the generated point to split the interval;
   (2) use the mean point of the interval; or
   (3) use the arcmean point.
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tabl_new( UNUR_DISTR* distribution );
/* Get default parameters for generator.                                     */

/*...........................................................................*/

int unur_tabl_set_variant_setup( UNUR_PAR *parameters, unsigned variant );
/* 
   Set variant for setup. Two modes are possible for @var{variant}:
   @code{1}: only use the equal area rule to construct hat.
   @code{2}: additionally split the intervals created by the equal
   area rule until the maximum number of intervals is reached or the
   ration between the area below squeeze and the area below the hat is
   exceeded.
   Default is variant @code{2}.
*/

int unur_tabl_set_variant_splitmode( UNUR_PAR *parameters, unsigned splitmode );
/* 
   There are three variants for adaptive rejection sampling. These
   differ in the way how an interval is split:
   splitmode @code{1}: use the generated point to split the interval.
   splitmode @code{2}: use the mean point of the interval.
   splitmode @code{3}: use the arcmean point.
   Default is splitmode @code{3}.
*/

int unur_tabl_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area below squeeze) / (area below hat).
   It must be a number between 0 and 1.
   When the ratio exceed the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.
   Default is ??.
*/

double unur_tabl_get_sqhratio( UNUR_GEN *generator );
/* 
   Get the current ratio (area below squeeze) / (area below hat)
   for the generator. (In case of an error 0 is returned.)
*/

int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals (default is ??).
   No construction points are added after the setup when the number of
   intervals suceeds @code{max_ivs}.
*/

int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
/* 
   Set parameter for equal area rule. During the setup a piecewise
   constant hat is constructed, such that the area below each of these
   pieces (strips) is the same and equal to the (given) area below the
   distribution times @var{fraction} (which must be greater than
   zero).
   Default is @code{0.25}.
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

int unur_tabl_set_slopes( UNUR_PAR *parameters, double *slopes, int n_slopes );
/* 
   Set slopes for the PDF.
   A slope <a,b> is an interval [a,b] or [b,a] where the PDF is
   monotone and PDF(a) >= PDF(b). 
   The list of slopes are given by an array @var{slopes} where each
   consecutive duples (i.e. @code{(slopes[0], slopes[1])}, 
   @code{(slopes[2], slopes[3])}, etc.) is one slopes.
   Slopes must be sorted (i.e. both @code{slopes[0]} and
   @code{slopes[1]} must not be greater than any entry of the slope
   @code{(slopes[2], slopes[3])}, etc.)
   and must not overlapping. Otherwise no slopes are set and
   @var{unur_errno} is set to @code{UNUR_ERR_PAR_SET}.

   Notice: @var{n_slopes} is the number of slopes (and not the length
   of the array @var{slopes}).

   Notice that setting slopes resets the given domain for the
   distribution.
*/

int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT). It must be greater than or equal to 0.
   If it is set to 0, then sequential search is used.
   Default is ??.
*/

int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set the left and right boundary of the computation interval.
   The piecewise hat is only constructed inside this interval. The
   region outside of this region must/should not be should be of
   computational importance.
   Of course +/- @code{UNUR_INFINITY} is not allowed.
   Default is ??.
*/

int unur_tabl_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/


