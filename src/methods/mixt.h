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
   =METHOD  MIXT  MIXTure of distributions

   =UP Meta_Methods

   =SPEED Set-up: fast, Sampling: depends on components

   =REINIT not implemented

   =DESCRIPTION
      MIXT allows to sample from a mixture of univariate
      distributions.

      Let @unurmath{f_1,\ldots,f_n} be PDFs of various distributions
      called the components and @unurmath{(p_1,\ldots,p_n)} be a
      probability vector. Then
      @unurmath{f(x) = p_1\cdot f_1(x) + \ldots + p_n\cdot f_n(x)}
      is the PDF of the so called mixture of these distributions.

      Method MIXT takes generator objects for the components and 
      a probability vector and creates a generator object for
      this mixture.

      The sampling part works as follows:

      @enumerate
      @item
      Generate an index @i{J} as the realisation of a discrete
      random variate with the given probability vector.
      @item
      Generate a random variate @i{X} with PDF @unurmath{f_J.}
      @end enumerate
   
      When the (interior of the) domains of the the components are
      disjoint then it is possible to sample from the mixture by
      inversion, provided that the following conditions are met:

      @itemize @minus
      @item
      The generator objects must use an inversion method for each
      component.
      @item
      The domain of the PDFs @unurmath{f_i} must be non-overlapping.
      @item
      The components must be order with respect to their domains.
      @end itemize
      
   =HOWTOUSE
      Create generator objects for the components of the mixture and
      store the corresponding pointers in an array.
      Store all probabilities an a double array of the same size.
      Create the parameter object for the generator of the mixture
      distribution by means of unur_mixt_new().

      The inversion method can be switched on by means of
      unur_mixt_set_useinversion() call.
      However, the conditions for this method must then be met. 
      Otherwise, initialization of the mixture object fails.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mixt_new( int n, const double *prob, UNUR_GEN **comp );
/* 
   Get default parameters for the generator for a mixture of the
   distributions given in the array @var{comp} (components) of length
   @var{n}. The probabilities are given by @var{prob}.

   The generators in @var{comp} must be objects for (continuous or
   discrete) univariate distributions 
*/

/*...........................................................................*/

int unur_mixt_set_useinversion( UNUR_PAR *parameters, int useinv );
/* 
   If @var{useinv} is TRUE, then the inversion method is used for
   sampling from the mixture distribution.

   However, the following conditions must be satisfied:

   @itemize @minus
   @item
   The generator objects must use an inversion method for each
   component.
   @item
   The domain of components must be non-overlapping.
   @item
   The components must be order with respect to their domains.
   @end itemize

   If one of these conditions is violated, then initialization of the
   mixture object fails.

   Default is FALSE.
*/

/* =END */

/*---------------------------------------------------------------------------*/
