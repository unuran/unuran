/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: discr.h                                                           *
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
   =NODE   DISCR   Discrete univariate distributions

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
   Create a new (empty) object for a univariate discrete distribution.
*/


/* ==DOC
   @subsubheading Essential parameters

   There are two interfaces for discrete univariate distributions:
   Either provide a (finite) probability vector (PV).
   Or provide a probability mass function (PMF). For the latter
   case there are also a couple of derived parameters that are not
   required when a PV is given.

   It is not possible to set both a PMF and a PV directly. However, a
   PV can be computed from a PMF by means of a
   unur_distr_discr_make_pv() call.
   If both a PV and a PMF are given in the distribution object it
   depends on the generation method which of these is used.
*/

int unur_distr_discr_set_pv( UNUR_DISTR *distribution, const double *pv, int n_pv );
/* 
   Set finite probability vector (PV) for a @var{distribution}. It is not
   necessary that the entries in the given PV sum to 1.
   @var{n_pv} must be positive. However, there is no testing
   whether all entries in @var{pv} are non-negative. 

   If no domain has been set, then the left boundary is set to
   @code{0}, by default. If @var{n_pv} is too large, e.g. because
   left boundary + @var{n_pv} exceeds the range of integers, 
   then the call fails. 

   Notice it is not possible to set both a PV and a PMF.
   (E.g., it is not possible to set a PV for a @var{distribution} from
   UNURAN library of standard distributions.)
*/

int unur_distr_discr_make_pv( UNUR_DISTR *distribution );
/* 
   Compute a PV when a PMF is given. However, when the
   domain is not given or is too large and the sum over the PMF is given
   then the (right) tail of the @var{distribution} is chopped off such that
   the probability for the tail region is less than 1.e-8.
   If the sum over the PMF is not given a PV of maximal length is
   computed.

   The maximal size of the created PV is bounded by the macro
   @code{UNUR_MAX_AUTO_PV} that is defined in @file{unuran_config.h}.

   If successful, the length of the generated PV is returned.
   If the sum over the PMF on the chopped tail is not neglible small
   (i.e. greater than 1.e-8 or unknown) than the 
   negative of the length of the PV is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.

   Notice that when a discrete distribution object is created from
   scratch, then the left boundary of the PV is set to @code{0} by
   default.

   If computing a PV fails for some reasons, @code{0} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
*/

int unur_distr_discr_get_pv( const UNUR_DISTR *distribution, const double **pv );
/* 
   Get length of PV of the @var{distribution} and set pointer
   @var{pv} to array of probabilities. If no PV is given,
   @code{0} is returned and @var{pv} is set to NULL.@*
   (It does not call unur_distr_discr_make_pv()!)
*/

int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
/* */

int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
/* 
   Set respective pointer to the PMF and the CDF of the @var{distribution}.
   These functions must be of type
   @code{double funct(int k, const UNUR_DISTR *distr)}.

   It is important to note that all these functions must return a
   result for all integers @var{k}. E.g., if the domain of a given
   PMF is the interval @{1,2,3,@dots{},100@}, than the given function
   must return @code{0.0} for all points outside this interval.

   It is not possible to change such a function. Once the PMF or
   CDF is set it cannot be overwritten. A new distribution object
   has to be used instead.

   Notice that it not possible to set both a PV and a PMF, i.e. it is not
   possible to use this call after a unur_distr_discr_set_pv() call.
*/


UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( const UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PMF and the CDF of the 
   @var{distribution}. The pointer is of type
   @code{double funct(int k, const UNUR_DISTR *distr)}.
   If the corresponding function is not available for the @var{distribution},
   the NULL pointer is returned.
*/


double unur_distr_discr_eval_pv(int k, const UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_pmf( int k, const UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_cdf( int k, const UNUR_DISTR *distribution );
/* 
   Evaluate the PV, PMF, and the CDF, respectively, at k.
   Notice that @var{distribution} must not be the NULL pointer.
   If no PV is set for the @var{distribution}, then
   unur_distr_discr_eval_pv() behaves like unur_distr_discr_eval_pmf().
   If the corresponding function is not available for the @var{distribution},
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.

   @emph{IMPORTANT:}
   In the case of a truncated standard distribution these calls always
   return the respective values of the @emph{untruncated} distribution!
*/


int unur_distr_discr_set_pmfstr( UNUR_DISTR *distribution, const char *pmfstr );
/* 
   This function provides an alternative way to set a PMF of the 
   @var{distribution}.
   @var{pmfstr} is a character string that contains the formula
   for the PMF, see @ref{StringFunct,,Function String}, for details.
   See also the remarks for the unur_distr_discr_set_pmf() call.

   It is not possible to call this funtion twice or to call this
   function after a unur_distr_discr_set_pmf() call.
*/

int unur_distr_discr_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
/* 
   This function provides an alternative way to set a CDF; analogously
   to the unur_distr_discr_set_pmfstr() call.
*/

char *unur_distr_discr_get_pmfstr( const UNUR_DISTR *distribution );
/* */

char *unur_distr_discr_get_cdfstr( const UNUR_DISTR *distribution );
/* 
   Get pointer to respective string for PMF and CDF
   of @var{distribution} that is given via the string interface.
   This call allocates memory to produce this string. It should be
   freed when it is not used any more.
*/

int unur_distr_discr_set_pmfparams( UNUR_DISTR *distribution, const double *params, int n_params );
/* 
   Set array of parameters for @var{distribution}. There is an upper limit
   for the number of parameters @var{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in @file{unuran_config.h}. (It is set to
   5 but can be changed to any appropriate nonnegative number.)
   If @var{n_params} is negative or exceeds this limit no parameters
   are copied into the @var{distribution} object and @code{unur_errno}
   is set to @code{UNUR_ERR_DISTR_NPARAMS}. 

   For standard distributions from the UNURAN library the parameters
   are checked. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_cont_set_domain() call.
   It these parameters are invalid, then no parameters are set and @code{0} 
   is returned.
   Notice that optional parameters are (re-)set to their default values if 
   not given for UNURAN standard distributions.

   @emph{Important:} Integer parameter must be given as doubles.
*/

int unur_distr_discr_get_pmfparams( const UNUR_DISTR *distribution, const double **params );
/* 
   Get number of parameters of the PMF and set pointer
   @var{params} to array of parameters. If no parameters are stored
   in the object, @code{0} is returned and @code{params} is set to
   NULL.
*/

int unur_distr_discr_set_domain( UNUR_DISTR *distribution, int left, int right );
/* 
   Set the left and right borders of the domain of the
   @var{distribution}. This can also be used to truncate an existing
   distribution. For setting the boundary to +/- infinity use
   @code{INT_MAX} and @code{INT_MIN}, respectively.
   If @var{right} is not strictly greater than @var{left} no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
   It is allowed to use this call to increase the domain.
   If the PV of the discrete distribution is used,
   than the right boudary is ignored (and internally set to 
   @var{left} + size of PV - 1).
   Notice that @code{INT_MAX} and @code{INT_MIN} are interpreted as
   (minus) infinity.
   Default is [@code{INT_MIN}, @code{INT_MAX}] when a PMF is used for generation,
   and [0, size of PV - 1] when a PV is used.
*/

int unur_distr_discr_get_domain( const UNUR_DISTR *distribution, int *left, int *right );
/* 
   Get the left and right borders of the domain of the
   @var{distribution}. If the domain is not set explicitly 
   the interval [@code{INT_MIN}, @code{INT_MAX}] is assumed and returned.
   When a PV is given then the domain is set automatically to 
   [@code{0},size of PV - 1].
*/


/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the essential
   parameters has been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_discr_set_mode( UNUR_DISTR *distribution, int mode );
/* 
   Set mode of @var{distribution}.
*/

int unur_distr_discr_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the @var{distribution}. This call works properly
   for distribution objects from the 
   UNURAN library of standard distributions 
   when the corresponding function is available.
   Otherwise a (slow) numerical mode finder is used. If it failes
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_discr_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of @var{distribution}. If the mode is not marked as known, 
   unur_distr_discr_upd_mode() is called to compute the mode. If this
   is not successful @code{INT_MAX} is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)
*/


int unur_distr_discr_set_pmfsum( UNUR_DISTR *distribution, double sum );
/* 
   Set the sum over the PMF. If @code{sum} is non-positive, no
   sum is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 

   For a distribution object created by the 
   UNURAN library of standard distributions you always should use
   the unur_distr_discr_upd_pmfsum(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution );
/*
   Recompute the sum over the PMF of the @var{distribution}. 
   In most cases the normalization constant is recomputed and thus the
   sum is 1. This call only works for distribution objects from the
   UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   The call does not work for distributions from the 
   UNURAN library of standard distributions with truncated
   domain when the CDF is not available.
*/

double unur_distr_discr_get_pmfsum( UNUR_DISTR *distribution );
/* 
   Get the sum over the PMF of the @var{distribution}. If this sum is
   not known, unur_distr_discr_upd_pmfsum() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
