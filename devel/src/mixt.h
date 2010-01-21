/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mixt.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for meta method MIXT                          *
 *         (MIXTure of distributions)                                        *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2010 Wolfgang Hoermann and Josef Leydold                  *
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
   =xxxMETHOD  MIXT  MIXTure of distributions

   =UP Meta_Methods [xy]

   =REINIT not implemented

   =DESCRIPTION
      MIXT allows to sample from a mixture of univariate
      distributions.

      Let @unurmath{f_1,\ldots,f_n} be PDFs of various distributions
      and @unurmath{(p_1,\ldots,p_n)} be a probability vector.
      Then by @unurmath{f(x) = \sum_{i=1}^n p_i\cdot f_i(x)}
      get the PDF of the so called mixture of these distributions.

      The methods takes generator objects for these distributions as
      well as probability vector and creates a generator object for
      this mixture.

      The sampling part works as follows:
      @enumerate
      @item
      Generate an index @i{J} as the realisation of a discrete
      random variate with given probability vector.
      @item
      Generate a random variate @i{X} with PDF @unurmath{f_J.}
      @end enumerate
   
      It is possible to recycle the uniform random number used for
      sampling index @i{J}. It is then used as the first uniform
      random number required for the sampling from the corresponding 
      distribution.
      This can be used for sampling from the mixture by inversion.
      However, for this feature the following conditions must be
      satisfied:
      @itemized
      @item
      An algorithm that implements the inversion method must be used
      for all distributions.
      @item
      The domain of the PDFs @unurmath{f_i} must be non-overlapping.
      @item
      The distributions must be order with respect to their domains.
      @end itemize
      (However, these conditions are never verified.)
      
   =HOWTOUSE
      Create generator objects for the components of the mixture and
      store the corresponding pointers in an array.
      Store all probabilities an a double array of the same size.
      Create the parameter object for the generator of the mixture
      distribution by means of unur_mixt_new().
      Recycling of uniform random random numbers can be enabled by
      means of unur_mixt_set_userecycle().
      However, we do not recommend recycling except for enabling
      inversion.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mixt_new( int n, const double *prob, const UNUR_GEN **genarray );
/* 
   Get default parameters for the generator for a mixture of the
   distributions given in the array @var{genarray} of length @var{n}.
   The probabilities are given by @var{prob}.

   The generators in @var{genarray} must be objects for univariate
   distributions (continuous or discrete).
*/

/*...........................................................................*/

int unur_mixt_set_userecycle( UNUR_PAR *parameters, int recycle );
/* 
   If @var{recycle} is TRUE, then the uniform random number used for
   sampling index @i{J} is recycled and used as the first uniform
   random number required for the sampling from the corresponding 
   distribution.
   This can be used for sampling from the mixture by inversion.

   We do not recommend recycling except for enabling inversion.

   Default is FALSE.
*/

/* =END */

/*---------------------------------------------------------------------------*/
