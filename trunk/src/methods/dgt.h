/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dis.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DIS                                *
 *         ((Discrete) Indexed Search (guide table))                         *
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
   =METHOD  DGT  (Discrete) Guide Table method (indexed search)

   =UP  Methods_for_DISCR

   =REQUIRED probability vector (PV)

   =SPEED Set-up: slow (linear with the vector-length), Sample: very fast

   =DESCRIPTION
      DGT samples from arbitrary but finite probability vectors. Random
      numbers are generated by the inversion method, i.e.,

      @enumerate
      @item
      Generate a random number U ~ U(0,1).
      @item
      Find largest integer I such that F(I) = P(X<=I) <= U.
      @end enumerate

      Step (2) is the crucial step. Using sequential search requires
      @i{O(E(X))} comparisons, where @i{E(X)} is the expectation of
      the distribution. Indexed search however uses a guide table to
      jump to some @i{I'} <= @i{I} near @i{I} to find @i{X} in constant
      time. Indeed the expected 
      number of comparisons is reduced to 2, when the guide table has the
      same size as the probability vector (this is the default). For
      larger guide tables this number becomes smaller (but is always
      larger than 1), for smaller tables it becomes larger. For the limit
      case of table size 1 the algorithm simply does sequential
      search. On the other hand the setup time for guide table is
      @i{O(N)} (for size 1 no preprocessing is required).
      Moreover for very large guide tables memory effects might
      even reduce the speed of the algorithm. So we do not recommend to
      use guide tables that are more than three times larger than the
      given probability vector. If only a few random numbers have to be
      generated, (much) smaller table sizes are better.
      The size of the guide table relative to the length of the given
      probability vector can be set by a unur_dgt_set_guidefactor() call.

      There exist two variants for the setup step which can be set by a 
      unur_dgt_set_variant() call: Variants 1 and 2.
      Variant 2 is faster but more sensitive to roundoff errors when the
      guide table is large. By default variant 2 is used for short
      probability vectors (@i{N}<1000) and variant 1 otherwise.
      
      By default the probability vector is indexed starting at
      @code{0}. However this can be changed in the distribution object by
      a unur_distr_discr_set_domain() call.

      The method also works when no probability vector but a PMF is
      given. However then additionally a bounded (not too large) domain
      must be given or the sum over the PMF (see
      unur_distr_discr_make_pv() for details).

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dgt_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_dgt_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set size of guide table relative to length of PV.
   Larger guide tables result in faster generation time but require a
   more expensive setup. Sizes larger than 3 are not recommended.
   If the relative size is set to 0, sequential search is used.

   Default is @code{1}. 
*/

int unur_dgt_set_variant( UNUR_PAR *parameters, unsigned variant );
/* 
   Set variant for setup step. Possible values are @code{1} or
   @code{2}.
   Variant @code{2} is faster but more sensitive to roundoff errors
   when the guide table is large. 
   By default variant @code{2} is used for short probability
   vectors (@i{N}<1000) and variant @code{1} otherwise.
*/

/* =END */

/*---------------------------------------------------------------------------*/





