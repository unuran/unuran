/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: norta.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NORTA                              *
 *         (NORmal To Anything)                                              *
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
   =METHOD  NORTA  NORmal To Anything

   =UP  Methods_for_CVEC

   =REQUIRED  rank correlation matrix, marginal distributions

   =SPEED Set-up: slow,
          Sampling: depends on dimension

   =REF  [HLD04: Sect.12.5.2, Alg.12.11.]

   =DESCRIPTION
      NORTA (NORmal to anything) is a model to get random vectors with
      given marginal distributions and rank correlation. 

      @strong{Important:} Notice that marginal distribution and (rank)
      correlation structure do not uniquely define a multivariate
      distribution. Thus the are many other (more or less sensible)
      models.

      In the NORTA model multinormal random variates with the given
      rank (Spearman's) correlations are generated. 
      In a second step the (standard normal distributed) marginal variates
      are transformed by means of the CDF of the normal distribution to get
      uniform marginals. The resulting random vectors have uniform
      marginals and the desired rank correlation between its components.
      Such a random vector is called 'copula'.

      By means of the inverse CDF the uniform marginals are then
      transformed into the target marginal distributions. This
      transformation does not change the rank correlation.
      
      For the generation of the multinormal distribution the
      (Spearman's) rank correlation matrix is transformed into the
      corresponding (Pearson) correlation matrix and the VMT method
      (see @ref{VMT}) is used to sample from the resulting multinormal
      distribution (by means of the Cholesky decomposition of the
      covariance matrix).
      It can happen that the desired rank correlation matrix is not
      feasible, i.e., it cannot occur as rank correlation matrix of a
      multinormal distribution. The resulting "covariance" matrix is
      not positive definite. In this an eigenvector correction
      method is used. Then all non-positive eigenvalues are set to a
      small positive value and hence the rank correlation matrix of the
      generated random vectors is "close" to the desired matrix.
      
   =HOWTOUSE
      Create a multivariate generator object and set marginal 
      distributions using unur_distr_cvec_set_marginals() , 
      unur_distr_cvec_set_marginal_array() , or 
      unur_distr_cvec_set_marginal_list().
      (Do not use the corresponding calls for the standard
      marginal distributions).
      
      If copulae are required (i.e. multivariate distributions with
      uniform marginals) such a generator object can be created by
      means of unur_distr_copula() .

      There are no optional parameters for this method.
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_norta_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/


/* =END */
/*---------------------------------------------------------------------------*/


