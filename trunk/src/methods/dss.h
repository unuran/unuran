/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dss.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DSS                                *
 *         ((Discrete) Sequential Search (guide table))                      *
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
   =METHOD  DSS  (Discrete) Sequential Search method

   =UP  Methods_for_DISCR

   =REQUIRED probability vector (PV) and sum over PV; or
             probability mass function(PMF), sum over PV and domain; or
             or cumulative distribution function (CDF)

   =SPEED Set-up: fast, Sampling: very slow (linear in expectation)

   =REF

   =DESCRIPTION
      DSS samples from arbitrary discrete distributions.
      Random numbers are generated by the inversion method, i.e.,

      @enumerate
      @item
      Generate a random number U ~ U(0,1).
      @item
      Find largest integer I such that F(I) = P(X<=I) <= U.
      @end enumerate

      Step (2) is the crucial step. Using sequential search requires
      @i{O(E(X))} comparisons, where @i{E(X)} is the expectation of
      the distribution. Thus this method is only recommended when only
      a few random variates from the given distribution are required.
      Otherwise, table methods like DGT (@pxref{DGT}) or DAU (@pxref{DAU})
      are much faster. These methods also need not the sum over the
      PMF (or PV) as input. On the other hand, however, these methods
      always compute a table.

      DSS runs with the PV, the PMF, or the CDF of the distribution.
      It uses actually uses the first one in this list (in this
      ordering) that could be found.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dss_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/* =END */

/*---------------------------------------------------------------------------*/
