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
  entirely by yourself. For this purpose unur_distr_new() returns an
  empty object of a given type (eg. univariate contiuous) which can be
  filled with the appropriate set calls.

  Notice that there are essential data about a distribution, eg. the
  p.d.f., a list of (shape, scale, location) parameters for the
  distribution, and the domain of (truncated) distributions.
  And there exist parameters that are/can be derived from these,
  eg. the mode of the distribution or the area below the given p.d.f.
  (which need not be normalized for may methods).
  UNURAN keeps track of parameters which are known. Thus if one of the
  essential parameters is changed all derived parameters are marked as
  unknown and must be set again if these are required for a the
  chosen generator method.
  For the objects provided by the ((=>) UNURAN library of standard
  distributions calls for updating these parameters exist (one for
  each parameter to avoid computational overhead since not all
  parameters are required for all generator methods).

  All the following calls only handle distribution object. Since every
  generator object has its own copy of a distribution object, there are 
  calls for a chosen method that change this copy of distribution object.
  NEVER extract the distribution object out of the generator object and 
  run one of the below set calls on it.
  (How should the poor generator object know what you have done?)
*/

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x010u,     /* univariate continuous distribution       */ 
  UNUR_DISTR_CEMP  = 0x011u,     /* empirical univ. cont. distr. (a sample)  */ 
  UNUR_DISTR_CVEC  = 0x110u,     /* mulitvariate continuous distribution     */ 
  UNUR_DISTR_DISCR = 0x020u,     /* univariate discrete distribution         */ 
  UNUR_DISTR_DEMP  = 0x021u,     /* empirical univariate discr. distribution */ 
};

/*---------------------------------------------------------------------------*/

/*
  Parameters common to all distributions.
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


int unur_distr_get_dim( UNUR_DISTR *distribution );
/* 
   Get number of components of random vector (its dimension).
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

int unur_distr_is_cvec( UNUR_DISTR *distribution );
/* 
   Test if distribution is a multivariate continuous distribution.
*/

int unur_distr_is_cemp( UNUR_DISTR *distribution );
/* 
   Test if distribution is an empirical univariate continuous distribution,
   i.e. a sample.
*/

int unur_distr_is_discr( UNUR_DISTR *distribution );
/* 
   Test if distribution is a univariate discrete distribution.
*/

int unur_distr_is_demp( UNUR_DISTR *distribution );
/* 
   Test if distribution is an empirical univariate discrete distribution.
*/


/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous distributions (CONT).
*/

UNUR_DISTR *unur_distr_cont_new( void );
/* 
   Create a new (empty) object for univariate continuous distribution.
*/

/* Essential parameters. */

int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *pdf );
/* See @code{unur_distr_cont_set_cdf}           */

int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dpdf );
/* See @code{unur_distr_cont_set_cdf}           */

int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *cdf );
/* 
   Set respective pointer to the probability density function (pdf),
   the derivative of the probability density function (dpdf) and the
   cumulative distribution function (cdf) of the distribution.
   The type of each of these functions must be of type
   double funct(double x, UNUR_DISTR *distr).

   It is not necessary that the given p.d.f. is normalized, i.e. the
   integral need not be 1 (but see below).
   Nevertheless the area below the given p.d.f. can be provided by a 
   unur_distr_cont_set_pdfarea() call.
   The given derivative must be the derivative of the given pdf of
   course.
   On the other hand the given c.d.f. must not be a multiple of the
   real c.d.f., i.e. lim_{x->infinity} cdf(x) = 1.

   It is important to note that all these functions must return a
   result for all floats @var{x}. Eg., if the domain of a given
   p.d.f. is the interval [-1,1], than the given function must return
   @code{0.0} for all points outside this interval.

   It is not possible to change such a function. Once the p.d.f. or
   c.d.f. is set it cannot be overwritten. A new distribution object
   has to be used instead.
*/

/*
   (If no parameters are used for the function use the NULL pointer
   for the second argument.)
*/

UNUR_FUNCT_CONT *unur_distr_cont_get_pdf( UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_get_cdf}           */

UNUR_FUNCT_CONT *unur_distr_cont_get_dpdf( UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_get_cdf}           */

UNUR_FUNCT_CONT *unur_distr_cont_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the p.d.f., the derivative of the 
   p.d.f. and the c.d.f. of the distribution. The pointer is of type
   double funct(double x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cont_eval_pdf( double x, UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_eval_cdf}           */

double unur_distr_cont_eval_dpdf( double x, UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_eval_cdf}           */

double unur_distr_cont_eval_cdf( double x, UNUR_DISTR *distribution );
/* 
   Evaluate the p.d.f., derivative of the p.d.f. and the c.d.f.,
   respectively, at x.
   Notice that @code{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/


int unur_distr_cont_set_pdfparams(UNUR_DISTR *distribution,double *params,int n_params);
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
   Get number of parameters of the p.d.f. and set pointer
   @code{params} to array of parameters. If no parameters are stored
   in the object, @code{0} is returned and @code{params} is set to
   NULL.
   Warning: Do not change the entries in @var{params}!
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
*/   
/*   
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
   area is set and @code{unur_errno} is set to @*
   @code{UNUR_ERR_DISTR_SET}. 
   
   For a distribution object created by the 
   (=>) UNURAN library of standard distributions you always should use
   the unur_distr_cont_upd_pdfarea(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_cont_upd_pdfarea( UNUR_DISTR *distribution );
/*
   Recompute the area below the p.d.f. of the distribution. 
   It only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   This call sets the normalization constant such that the given
   p.d.f. is the derivative of a given c.f.d., i.e. the area is 1.
   However for truncated distribution the area smaller than 1.
*/

double unur_distr_cont_get_pdfarea( UNUR_DISTR *distribution );
/* 
   Get the area below the p.d.f of the distribution. If this area is
   not known,@* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous order statistics (CORDER).
*/

UNUR_DISTR *unur_distr_corder_new( UNUR_DISTR *distr, int n, int k );
/* 
   Create an object for order statistics of for a sample size
   @var{n} and rank @var{k}.
   @var{distr} must be a pointer to a univariate continuous distribution.
   The area below the p.d.f. is set to the (possibly unknown) area
   given for the underlying distribution.

   The result is of the same type as of unur_distr_cont_new() calls.
   (Except it cannot be used to for making the order statistics of the
   order statistics.)
*/

/* Essential parameters. */

/* 
   The result of unur_distr_corder_new() can be handled with
   unur_distr_cont_... calls with the following exceptions:

   the p.d.f.s and c.d.f. cannot be set or changed.

   unur_distr_cont_set_pdfparams() changes the parameters of the
      underlying distribution.

   unur_distr_cont_upd_mode() does not work.

   Additionally the following routines can be used.
*/
  

UNUR_DISTR *unur_distr_corder_get_distribution( UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/

int unur_distr_corder_set_rank( UNUR_DISTR *distribution, int n, int k );
/* 
   Change sample size $var{n} and rank @var{k} of order statistics.
   In case of invalid data, no parameters are changed and @code{0} is
   returned.
   The area below the p.d.f. can be set to that of the underlying
   distribution by a unur_distr_upd_pdfarea() call.
*/

int unur_distr_corder_get_rank( UNUR_DISTR *distribution, int *n, int *k );
/* 
   Get sample size $var{n} and rank @var{k} of order statistics.
   In case of error @code{0} is returned.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling multivariate continuous distributions (CVEC).
*/

UNUR_DISTR *unur_distr_cvec_new( int dim );
/* 
   Create a new (empty) object for multivariate continuous
   distribution. @var{dim} is the number of components of the random
   vector (i.e. its dimension). It must be at least 2; otherwise
   unur_distr_cont_new() should be used to create an object for a
   univariate distribution.
*/

/* Essential parameters. */

int unur_distr_cvec_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *pdf );
/* 
   Set respective pointer to the probability density function (pdf) of
   the distribution. The type of this function must be of type
   double funct(double *x, UNUR_DISTR *distr),
   where @var{x} must be a pointer to a double array of appropriate
   size (i.e. of the same size as given to the unur_distr_cvec_new()
   call).

   It is not necessary that the given p.d.f. is normalized and can be
   any (positive) multiple of the p.d.f., i.e. the integral need not
   be 1. Nevertheless it can be provided by a
   unur_distr_cvec_set_pdfarea() call.
*/

int unur_distr_cvec_set_dpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dpdf );
/* 
   Set pointer to the gradiant of the probability density function
   (pdf). The type of this function must be
   int funct(double *result, double *x, UNUR_DISTR *distr),
   where @var{result} and @var{x} must be pointer to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).
   The gradiant of the p.d.f. is stored in the array @var{result}.
   The function should return @code{0} in case of an error and must
   return a non-zero value otherwise.

   The given function must be proved the gradiant of the function
   given by a unur_distr_cvec_set_pdf() call.
*/

UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the p.d.f. of the distribution. The
   pointer is of type 
   double funct(double *x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the gradiant of the p.d.f. of the
   distribution. The pointer is of type 
   int double funct(double *result, double *x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cvec_eval_pdf( double *x, UNUR_DISTR *distribution );
/* 
   Evaluate the p.d.f. of the distribution at @var{x}.
   @var{x} must be a pointer to a double arrays of appropriate size
   (i.e. of the same size as given to the unur_distr_cvec_new() call)
   that contains the vector for which the function has to be evaluated.

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_eval_dpdf( double *result, double *x, UNUR_DISTR *distribution );
/* 
   Evaluate the gradiant of the p.d.f. of the distribution at @var{x}.
   The result is stored in the double array @var{result}.
   Both @var{result} and @var{x} must be pointer to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   distribution, @code{0} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA} (@var{result} is left unmodified).
*/

int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, int par, double *params, int n_params );
/* 
   Set parameter with number @var{par}. 
   @var{par} indicates directly which of the parameters is set and
   must be a number between @code{0} and @code{UNUR_DISTR_MAXPARAMS}-1
   (the upper limit of possible parameters defined in
   unuran_config.h; it is set to 5 but can be changed to any 
   appropriate nonnegative number.)

   All parameters are given by the array @var{params} of size
   @var{n_params}. Thus for a single parameter an array of size 1 and
   @var{n_params} has to be used. 
   For a vector @var{n_params} has to be set to the size of the array.
   An (n x m)-matrix has to be stored in an array of length
   @var{n_params} = n m; where the rows of the matrix are stored
   consecutively in this array.

   When more than one type of parameters are used (e.g. the mean
   vector and the covariance matrix) then there for each of these 
   unur_distr_cvec_set_pdfparams() are required.

   Due to great variety of possible parameters for a multivariate
   distribution there is no simpler interface.

   If an error occurs no parameters are copied into the parameter
   object @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_get_pdfparams( UNUR_DISTR *distribution, int par, double **params );
/* 
   Get parameter of the p.d.f. with number @var{par}.
   The pointer to the parameter array is stored in @var{params}, its
   size is returned by the function.
   If the requested parameter is not set, then @code{0} is returned
   and @code{params} is set to NULL.

   Warning: Do not change the entries in @var{params}!
*/

/* 
   Derived parameters.
*/   
/*   
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_cvec_set_mode( UNUR_DISTR *distribution, double *mode );
/* 
   Set mode of distribution.
*/

int unur_distr_cvec_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the distribution. This call only works for
   distribution objects from the (=>) UNURAN library of standard
   distributions when the corresponding function is available.
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

double *unur_distr_cvec_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of distribution. The function return a pointer to an array
   of the size returned by unur_distr_get_dim().
   If the mode is not marked as known, 
   unur_distr_cvec_upd_mode() is called to compute the mode. If this
   is not successful NULL is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)

   Warning: Do not modify the array that holds the mode!
*/

int unur_distr_cvec_set_pdfvol( UNUR_DISTR *distribution, double volume );
/* 
   Set the volume below the p.d.f. If @var{vol} is non-positive, no
   volume is set and @code{unur_errno} is set to @*
   @code{UNUR_ERR_DISTR_SET}. 
*/

int unur_distr_cvec_upd_pdfvol( UNUR_DISTR *distribution );
/*
   Recompute the volume below the p.d.f. of the distribution. 
   It only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 
   
   This call sets the normalization constant such that the volume is
   1.
*/

double unur_distr_cvec_get_pdfvol( UNUR_DISTR *distribution );
/* 
   Get the volume below the p.d.f of the distribution. If this volume is
   not known,@* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical univariate continuous distributions (CEMP).
*/

UNUR_DISTR *unur_distr_cemp_new( void );
/* 
   Create a new (empty) object for empirical univariate continuous distribution.
*/

/* Essential parameters */

int unur_distr_cemp_set_data( UNUR_DISTR *distribution, double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
*/


int unur_distr_cemp_get_data( UNUR_DISTR *distribution, double **sample );
/* 
   Get number of samples and set pointer @code{sample} to array of
   observations. If no sample has beeb given,
   @code{0} is returned and @code{sample} is set to NULL.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate discrete distributions (DISCR).
*/

UNUR_DISTR *unur_distr_discr_new( void );
/* 
   Create a new (empty) object for univariate discrete distribution.
*/


/* Essential parameters. */

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
   c.d.f. is set it cannot be overwritten. A new distribution object
   has to be used instead.

*/

UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( UNUR_DISTR *distribution );
/* See @code{unur_distr_discr_get_cdf}           */

UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the p.m.f. and the c.d.f. of the 
   distribution. The pointer is of type
   double funct(int k, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_discr_eval_pmf( int k, UNUR_DISTR *distribution );
/* See @code{unur_distr_discr_eval_cdf}           */

double unur_distr_discr_eval_cdf( int k, UNUR_DISTR *distribution );
/* 
   Evaluate the p.m.f., and the c.d.f., respectively, at k.
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

/* 
   Derived parameters.
*/   
/* 
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_discr_set_pmfarea(UNUR_DISTR *distribution, double area);
/* 
   Set the area below the p.m.f. If @code{area} is non-positive, no
   area is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 

   For a distribution object created by the 
   (=>) UNURAN library of standard distributions you always should use
   the unur_distr_discr_upd_pmfarea(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_discr_upd_pmfarea( UNUR_DISTR *distribution);
/*
   Recompute the area below the p.m.f. of the distribution. 
   In most cases the normalization constant is recompute and thus the
   area is 1. This call only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 
*/

double unur_distr_discr_get_pmfarea(UNUR_DISTR *distribution);
/* 
   Get the area below the p.m.f of the distribution. If this area is
   not known, @* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical univariate discrete distributions (DEMP).
*/

UNUR_DISTR *unur_distr_demp_new( void );
/* 
   Create a new (empty) object for empirical univariate discrete distribution.
*/

/* Essential parameters */

int unur_distr_demp_set_prob( UNUR_DISTR *distribution, double *prob, int n_prob );
/* 
   Set finite probability vector for a distribution. It is not
   necessary that the entries in the given probability vector sum to
   1. @code{n_prob} must be positive. However there is no testing
   whether all entries in @code{prob} are non-negative. 
*/


int unur_distr_demp_get_prob( UNUR_DISTR *distribution, double **prob );
/* 
   Get length of probability vector of the distribution and set pointer
   @code{prob} to array of probabilities. If no probability vector is given,
   @code{0} is returned and @code{prob} is set to NULL.
*/

/*---------------------------------------------------------------------------*/
