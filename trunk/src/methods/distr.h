/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects         *
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
  For handling distributions objects of type @code{UNUR_DISTR} are
  used. All data about a distribution are stored on this object. 
  UNURAN provides functions that return such objects for standard
  distributions ((=>) UNURAN library of standard distributions).
  It is then possible to change this distribution object by various
  set calls. Moreover it is possible to build a distribution object
  entirely by yourself. For this purpose unur_distr_new() returna an
  empty object of a given type (eg. univariate contiuous) which can be
  filled with the appropriate set calls.

  Notice that there are essential data about a distribution, eg. the
  p.d.f., a list of (shape, scale, location) parameters for the
  distribution, and the domain of (truncated) distributions.
  And there exist parameters that are/can be derived from these,
  eg. the mode of the distribution or the normalization constant.
  UNURAN keeps track of parameters which are known. Thus if one of the
  essential parameters is changed all derived parameters are marked as
  unknown and must be set again if these are required for a the
  chosen generator method.
  For the objects provided by the ((=>) UNURAN library of standard
  distributions calls for updating these parameters exist (one for
  each parameter to avoid computational overhead since not all
  parameters are required for all generator methods).

  All the following calls only handle distribution object. Since every
  generator object has its own copy of a distribution object, there is
  no need for calls that change parameters of such distributions. So
  it is there exist such calls on a per method basis. Each if the
  domain of the distribution of a generator object (eg. of method
  NINV) should be changed then use (in this example)
  unur_ninv_chg_domain(). NEVER extract the distribution object out of
  the generator object and run unur_distr_cont_set_domain() instead!
  (How should the poor generator object know what you have done?)
*/

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x001u,        /* univariate continuous distribution    */ 
  UNUR_DISTR_DISCR = 0x002u,        /* univariate discrete distribution      */ 
};

/*---------------------------------------------------------------------------*/

/*
  Parameters common to all distributions.
*/

UNUR_DISTR *unur_distr_new( unsigned int type );
/* 
   Create a new (empty) distribution object.
   @code{type} indicates the type of the distribution. Currently the
   following types are available:

   UNUR_DISTR_CONT    ... univariate continuous distribution
   UNUR_DISTR_DISCR   ... univariate discrete distribution

*/


void unur_distr_free( UNUR_DISTR *distribution );
/* 
   Destroy a distribution object.
*/


int unur_distr_set_name( UNUR_DISTR *distribution, const char *name );
const char *unur_distr_get_name( UNUR_DISTR *distribution );
/* 
   Set and get name of distribution.
*/


unsigned int unur_distr_get_type( UNUR_DISTR *distribution );
/*
  Get type of distribution. See description of unur_distr_new() for
  possible types. Alternatively the unur_distr_is_<TYPE>() calls can
  be used.
*/

int unur_distr_is_cont( UNUR_DISTR *distribution );
/* 
   Test if distribution is a univariate continuous distribution.
*/

int unur_distr_is_discr( UNUR_DISTR *distribution );
/* 
   Test if distribution is a univariate discrete distribution.
*/


/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous distributions.
*/

/* Essential parameters. */

int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, void *pdf );
int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, void *dpdf );
int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, void *cdf );
/* 
   Set respective pointer to the probability density function (pdf),
   the derivative of the probability density function (dpdf) and the
   cumulative distribution function (cdf) of the distribution.
   The type of each of these functions must be either
   double funct(double x) or 
   double funct(double x, UNUR_DISTR *distr).
   The second form is necessary when parameters for the p.d.f. are
   be used.  */


void *unur_distr_cont_get_pdf( UNUR_DISTR *distribution );
void *unur_distr_cont_get_dpdf( UNUR_DISTR *distribution );
void *unur_distr_cont_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the p.d.f., the derivative of the 
   p.d.f. and the c.d.f. of the distribution. The pointer is of type
   double funct(double x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cont_eval_pdf( double x, UNUR_DISTR *distribution );
double unur_distr_cont_eval_dpdf( double x, UNUR_DISTR *distribution );
double unur_distr_cont_eval_cdf( double x, UNUR_DISTR *distribution );
/* 
   Evaluate the p.d.f., derivative of the p.d.f. and the c.d.f.,
   respectively, at x.
   Notice that @code{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/


int unur_distr_cont_set_pdfparams( UNUR_DISTR *distribution, double *params, int n_params );
/* 
   Set array of parameters for distribution. There is an upper limit
   for the number of parameters @code{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in unuran_config.h. (It is set to
   5 but can be changed to any appropriate nonnegative number.)
   If @code{n_params} is negative or exceeds this limit no parameters
   are copied into the distribution object and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_distr_cont_get_pdfparams( UNUR_DISTR *distribution, double **params );
/* 
   Get number of parameters of the p.d.f. and sets pointer
   @code{params} to array of parameters. If no parameters are stored
   in the object, @code{0} is returned and @code{params} is set to
   NULL.
*/

int unur_distr_cont_set_domain( UNUR_DISTR *distribution, double left, double right );
/* 
   Set the left and right borders of the domain of the
   distribution. This can also be used truncate an existing
   distribution. For setting the boundary to +/- infinity use +/-
   @code{UNUR_INFINITY}.
   If @code{right} is not strictly greater than @code{left} no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
   It is allowed to use this call to increase the domain.
*/

int unur_distr_cont_get_domain( UNUR_DISTR *distribution, double *left, double *right );
/* 
   Get the left and right borders of the domain of the
   distribution. If the domain is not set explicitly 
   +/- @code{UNUR_INFINITY} is assumed and returned.
   There is no error reported in this case.
*/


/* 
   Derived parameters.
   
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_cont_set_mode( UNUR_DISTR *distribution, double mode );
/* 
   Set mode of distribution.
*/

int unur_distr_cont_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the distribution. This call only works for
   distribution objects from the (=>) UNURAN library of standard
   distributions when the corresponding function is available.
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/


double unur_distr_cont_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of distribution. If the mode is not marked as known, 
   unur_distr_cont_upd_mode() is called to compute the mode. If this
   is not successful @code{UNUR_INFINITY} is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)
*/


int unur_distr_cont_set_pdfarea( UNUR_DISTR *distribution, double area );
/* 
   Set the area below the p.d.f. If @code{area} is non-positive, no
   area is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 
*/

int unur_distr_cont_upd_pdfarea( UNUR_DISTR *distribution );
/*
   Recompute the area below the p.d.f. of the distribution. 
   In most cases the normalization constant is recompute and thus the
   area is 1. This call only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 
*/

double unur_distr_cont_get_pdfarea( UNUR_DISTR *distribution );
/* 
   Get the area below the p.d.f of the distribution. If this area is
   not known, unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/


/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate discrete distributions.
*/

/* Essential parameters */

int unur_distr_discr_set_prob( UNUR_DISTR *distribution, double *prob, int n_prob );
/* 
   Set finite probability vector for a distribution. It is not
   necessary that the entries in the given probability vector sum to
   1. @code{n_prob} must be positive. However there is no testing
   whether all entries in @code{prob} are non-negative. 
*/


/* 
   Derived parameters.
   
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

/*---------------------------------------------------------------------------*/










