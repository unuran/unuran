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
   =DISTR   DISCR   Discrete univariate distributions

   =UP Distribution_objects [50]

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate discrete distributions (DISCR).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_discr_new( void );
/* 
   Create a new (empty) object for univariate discrete distribution.
*/


/* ==DOC
   @subsubheading Essential parameters

   There are two interfaces for discrete univariate distributions:
   Either provide a (finite) probability vector.
   Or provide a probability mass function (PMF). For the latter
   case there exist also a couple of derived parameters that are not
   required when a probability vector is given.

   If both a probability vector and a PMF is given it depends on
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

   If no domain has not been set, then the left boundary is set to
   @code{0}, by default. If @var{n_prob} is too large, e.g. because
   left boundary + @var{n_prob} exceeds the range of integers, 
   then the call fails. 

   Notice it not possible to set both a PV and a PMF.
   (I.e., it is not possible to set a PV for a distribution from
   (=>) UNURAN library of standard distributions.)
*/

int unur_distr_discr_make_prob( UNUR_DISTR *distribution );
/* 
   Compute a probability vector when a PMF is given. However it only
   works when the domain of the distribution is bounded and not too
   large or the sum over the PMF is given. The maximal size of the
   created PV is bounded by the macro @code{UNUR_MAX_AUTO_PV} that is
   defined in @file{unuran_config.h}.
   If successful the length of the generated probablity vector is
   returned.
   If the sum over the PMF on the chopped tail is not neglible small,
   then the last entry of the computed PV contains this sum and the
   negative of the length of the PV is returned.

   Notice that when a discrete distribution object is created from
   scratch, then the left boundary is set to @code{INT_MIN}. Therefore
   it is necessary to change the this domain by a
   unur_distr_discr_set_domain() call, since otherwise it is not
   possible to compute the PV.

   If computing a PV fails for some reasons, @code{0} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.

   Notice that it is not possible to execute this call when the
   distribution object already contains a PV.
*/

int unur_distr_discr_get_prob( UNUR_DISTR *distribution, double **prob );
/* 
   Get length of probability vector of the distribution and set pointer
   @code{prob} to array of probabilities. If no probability vector is given,
   @code{0} is returned and @code{prob} is set to NULL.
   (It does not call unur_distr_discr_make_prob()!)
*/


int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
/* */

int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
/* 
   Set respective pointer to the probability mass function (pmf) and the
   cumulative distribution function (cdf) of the distribution.
   The type of each of these functions must be of type
   double funct(int k, UNUR_DISTR *distr).

   It is important to note that all these functions must return a
   result for all integers @var{k}. Eg., if the domain of a given
   PMF is the interval @{1,2,3,...,100@}, than the given function
   must return @code{0.0} for all points outside this interval.

   It is not possible to change such a function. Once the PMF or
   CDF is set it cannot be overwritten. A new distribution object
   has to be used instead.

   Notice it not possible to set both a PMF and a PV.
*/


UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PMF and the CDF of the 
   distribution. The pointer is of type
   double funct(int k, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/


double unur_distr_discr_eval_prob(int k, UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_pmf( int k, UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_cdf( int k, UNUR_DISTR *distribution );
/* 
   Evaluate the probabilty vector, PMF, and the CDF, respectively, at k.
   Notice that @code{distribution} must not be the NULL pointer.
   If no probability vector is available @code{unur_distr_discr_eval_prob}
   will evaluate the PMF instead.
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
   Get number of parameters of the PMF and set pointer
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


/* ==DOC
   @subsubheading Derived parameters

   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_discr_set_mode( UNUR_DISTR *distribution, int mode );
/* 
   Set mode of distribution.
*/

int unur_distr_discr_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the distribution. This call works properly
   for distribution objects from the 
   (=>) UNURAN library of standard distributions 
   when the corresponding function is available.
   Otherwise a (slow) numerical mode finder is used. If it failes
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_discr_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of distribution. If the mode is not marked as known, 
   unur_distr_discr_upd_mode() is called to compute the mode. If this
   is not successful @code{INT_MAX} is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)
*/


int unur_distr_discr_set_pmfsum(UNUR_DISTR *distribution, double sum);
/* 
   Set the sum over the PMF. If @code{sum} is non-positive, no
   sum is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 

   For a distribution object created by the 
   (=>) UNURAN library of standard distributions you always should use
   the unur_distr_discr_upd_pmfsum(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution);
/*
   Recompute the sum over the PMF of the distribution. 
   In most cases the normalization constant is recompute and thus the
   sum is 1. This call only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   The call does not work for distributions from the 
   (=>) UNURAN library of standard distributions whith truncated
   domain when the CDF is not available.
*/

double unur_distr_discr_get_pmfsum(UNUR_DISTR *distribution);
/* 
   Get the sum over the PMF of the distribution. If this sum is
   not known, unur_distr_discr_upd_pmfsum() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/


/* REMOVE */
int _unur_distr_discr_find_mode( struct unur_distr *distr );
