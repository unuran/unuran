/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hinv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HINV                               *
 *         (Hermite interpolation based INVersion of CDF)                    *
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
   =METHOD HINV  Hermite interpolation based INVersion of CDF

   =UP  Methods_for_CONT

   =REQUIRED CDF

   =OPTIONAL PDF, dPDF

   =SPEED Set-up: (very) slow, Sampling: (very) fast

   =DESCRIPTION
      HINV is a variant of numerical inversion, where the inverse CDF
      is approximated using Hermite interpolation. These splines have
      to be computed in a setup step. However, it only works for
      distributions with bounded domain; for distributions with
      unbounded domain the tails are chopped off such that the
      probability for the tail regions is small compared to the given
      u-resulution. For finding these cut points the algorithm starts 
      with the region @code{[-1.e20,1.e20]}. For the exceptional case
      where this might be too small (or one knows this region and
      wants to avoid this search heuristics) it can be chanced using
      the unur_hinv_set_cpoints() call.
      
      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_hinv_chg_truncated() call.

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      can be made as small as desired by means of the
      unur_hinv_set_u_resolution(). Notice that small values of the
      u-resultion increases the cost for the setup step.
      
      Sometimes it might be necessary to set some design points for
      the computing the Hermite interpolation. Such points may be
      points where the density is not differentialble or is an
      inflection points. Notice that there is no necessity to do
      so. However, if you do not provide these points to the algorithm
      there might be a (very) small chance that the approximation
      error is larger than the given u-resolution, or that the
      required number of intervals is larger than necessary.
      Setting such design points can be done using the
      unur_hinv_set_cpoints() call.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_hinv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hinv_set_order( UNUR_PAR *parameters, int order);
/* 
   Set order of Hermite interpolation. Valid orders are
   @code{1}, @code{3}, and @code{5}.
   Notice that @var{order} greater than 1 requires the density 
   of the distribution, and @var{order} greater than 3 even
   requires the derivative of the density.

   Default is @code{3} if the density is given
   and @code{1} otherwise.
*/

int unur_hinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal error in u-direction.
   Default is @code{10^-8}.
*/

int unur_hinv_set_cpoints( UNUR_PAR *parameters, const double *stp, int n_stp );
/* 
   Set starting construction points (nodes) for Hermite interpolation. 

   @emph{Important}: Notice that the given points must be in
   increasing order and they must be disjoint. There also must be at
   least two such points. Otherwise the points cannot be used and
   @code{unur_errno} is set to @code{UNUR_ERR_PAR_SET}.

   Notice that the leftmost and the rightmost point are used as
   boundary point of the region where the Hermite approximation is
   computated. Thus this call can also be used to define the region
   where computational important, i.e., to define the cut points for
   the tails.

   @emph{Important}: The boundary point of the computational region
   must be given in this list!
*/

int unur_hinv_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is @code{1}.
*/

int unur_hinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals. No generator object is created if
   the necessary number of intervals for the Hermite interpolation 
   exceed @var{max_ivs}. It is used to prevent the algorithm to eat up
   all memory in badly shaped CDFs.

   Default is @code{1.e6}.
*/

int unur_hinv_get_n_intervals( const UNUR_GEN *generator );
/* 
   Get number of intervals used for Hermite interpolation in 
   generator object.
   It return @code{0} in case of an error.
*/

int unur_hinv_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Changes the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call. The tables splines are not recomputed. Thus
   it might happen that the relative error for the generated variates
   from the truncated distribution is greater than the bound given by a
   unur_hinv_chg_truncated() call.
*/

/* =END */
/*---------------------------------------------------------------------------*/
