/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_discr.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  DISCR  (discrete univariate distribution)                   *
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

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate discrete distributions (DISCR).
*/

UNUR_DISTR *unur_distr_discr_new( void );
/* 
   Create a new (empty) object for univariate discrete distribution.
*/


/* Essential parameters. */

/* 
   There are two interfaces for discrete univariate distributions:
   Either provide a (finite) probability vector.
   Or provide a probability mass function (p.m.f.). For the latter
   case there exist also a couple of derived parameters that are not
   required when a probability vector is given.

   If both a probability vector and a p.m.f. is given it depends on
   the generation method which of these is used.
   Notice that there might exist some confusion if both are given but
   describe a different distribution! (There is no checking against
   this inconsistency!)
*/

int unur_distr_discr_set_prob( UNUR_DISTR *distribution, double *prob, int n_prob );
/* 
   Set finite probability vector for a distribution. It is not
   necessary that the entries in the given probability vector sum to
   1. @code{n_prob} must be positive. However there is no testing
   whether all entries in @code{prob} are non-negative. 
*/


int unur_distr_discr_get_prob( UNUR_DISTR *distribution, double **prob );
/* 
   Get length of probability vector of the distribution and set pointer
   @code{prob} to array of probabilities. If no probability vector is given,
   @code{0} is returned and @code{prob} is set to NULL.
*/

int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
/* See @code{unur_distr_discr_set_cdf}           */

int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
/* 
   Set respective pointer to the probability mass function (pmf) and the
   cumulative distribution function (cdf) of the distribution.
   The type of each of these functions must be of type
   double funct(int k, UNUR_DISTR *distr).

   It is important to note that all these functions must return a
   result for all integers @var{k}. Eg., if the domain of a given
   p.m.f. is the interval @{1,2,3,...,100@}, than the given function
   must return @code{0.0} for all points outside this interval.

   It is not possible to change such a function. Once the p.m.f. or
   CDF is set it cannot be overwritten. A new distribution object
   has to be used instead.

*/

UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( UNUR_DISTR *distribution );
/* See @code{unur_distr_discr_get_cdf}           */

UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the p.m.f. and the CDF of the 
   distribution. The pointer is of type
   double funct(int k, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_discr_eval_pmf( int k, UNUR_DISTR *distribution );
/* See @code{unur_distr_discr_eval_cdf}           */

double unur_distr_discr_eval_cdf( int k, UNUR_DISTR *distribution );
/* 
   Evaluate the p.m.f., and the CDF, respectively, at k.
   Notice that @code{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/


int unur_distr_discr_set_pmfparams( UNUR_DISTR *distribution, double *params, int n_params );
/* 
   Set array of parameters for distribution. There is an upper limit
   for the number of parameters @code{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in unuran_config.h. (It is set to
   5 but can be changed to any appropriate nonnegative number.)
   If @code{n_params} is negative or exceeds this limit no parameters
   are copied into the distribution object and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}. There is no way to use parameters of type int.
*/

int unur_distr_discr_get_pmfparams(UNUR_DISTR *distribution,double **params);
/* 
   Get number of parameters of the p.m.f and set pointer
   @code{params} to array of parameters. If no parameters are stored
   in the object, @code{0} is returned and @code{params} is set to
   NULL.
*/

int unur_distr_discr_set_domain( UNUR_DISTR *distribution, int left, int right );
/* 
   Set the left and right borders of the domain of the
   distribution. This can also be used truncate an existing
   distribution. For setting the boundary to +/- infinity use
   @code{INT_MAX} and @code{INT_MIN}, respectively.
   If @code{right} is not strictly greater than @code{left} no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
   It is allowed to use this call to increase the domain.
   If the probability vector of the discrete distribution is used,
   than the right boudary is ignored (and internally set to 
   left + size of PV - 1).
   Notice that @code{INT_MAX} and @code{INT_MIN} are interpreted as
   (minus) infinity.
   Default is [INT_MIN, INT_MAX] when a PMF is used for generation,
   and [0, size of PV - 1] when a probability vector is used.
*/

int unur_distr_discr_get_domain( UNUR_DISTR *distribution, int *left, int *right );
/* 
   Get the left and right borders of the domain of the
   distribution. If the domain is not set explicitly 
   the interval [INT_MIN, INT_MAX] is assumed and returned.
   There is no error reported in this case.
*/


/* 
   Derived parameters.
*/   
/* 
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_discr_set_pmfsum(UNUR_DISTR *distribution, double sum);
/* 
   Set the sum over the p.m.f. If @code{sum} is non-positive, no
   sum is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 

   For a distribution object created by the 
   (=>) UNURAN library of standard distributions you always should use
   the unur_distr_discr_upd_pmfsum(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution);
/*
   Recompute the sum over the p.m.f. of the distribution. 
   In most cases the normalization constant is recompute and thus the
   sum is 1. This call only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 
*/

double unur_distr_discr_get_pmfsum(UNUR_DISTR *distribution);
/* 
   Get the sum over the p.m.f of the distribution. If this sum is
   not known, unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/*---------------------------------------------------------------------------*/
