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
  entirely from scratch. For this purpose there exists an
  unur_distr_<type>_new() call for each object type that returns an
  empty object of this type (eg. univariate contiuous) which can be
  filled with the appropriate set calls.

  Notice that there are essential data about a distribution, eg. the
  PDF, a list of (shape, scale, location) parameters for the
  distribution, and the domain of (truncated) distributions.
  And there exist parameters that are/can be derived from these,
  eg. the mode of the distribution or the area below the given PDF
  (which need not be normalized for may methods).
  UNURAN keeps track of parameters which are known. Thus if one of the
  essential parameters is changed all derived parameters are marked as
  unknown and must be set again if these are required for a the
  chosen generator method.

  The library can handle truncated distributions, that is,
  distribution that are derived from (standard) distribution by simply
  restrict its domain to a subset. There is a subtle difference
  between changing the domain by a unur_distr_cont_set_domain() call
  and changing the (truncated) domain for an existing generator object
  by a unur_<method>_chg_truncated() call (if available).
  The domain of the given distribuiton is used to create the generator
  object. This is always handled as the domain of a non-truncated
  distribution (although it really was derived from one UNURAN
  standard distributions by resetting the domain). This domain can
  then be restricted to a subset for the generator object. 
  Generator methods that require a recreation of the generator object
  Notice: when the domain is changed have a unur_<method>_chg_domain()
  instead. For this call there are of course no restrictions on the
  given domain (i.e., it is possible to increase the domain of the
  distribution).

  For the objects provided by the ((=>) UNURAN library of standard
  distributions calls for updating these parameters exist (one for
  each parameter to avoid computational overhead since not all
  parameters are required for all generator methods).

  All the following calls only handle distribution object. Since every
  generator object has its own copy of a distribution object, there are 
  calls for a chosen method that change this copy of distribution object.
  NEVER extract the distribution object out of the generator object and 
  run one of the below set calls on it.
  (How should the poor generator object know what has happend?)
*/

/*---------------------------------------------------------------------------*/
/* types of distribtuions                                                    */

enum {
  UNUR_DISTR_CONT  = 0x010u,     /* univariate continuous distribution       */ 
  UNUR_DISTR_CEMP  = 0x011u,     /* empirical univ. cont. distr. (a sample)  */ 
  UNUR_DISTR_CVEC  = 0x110u,     /* mulitvariate continuous distribution     */ 
  UNUR_DISTR_CVEMP = 0x111u,     /* empirical multiv. cont. distr. (sample)  */ 
  UNUR_DISTR_DISCR = 0x020u,     /* univariate discrete distribution         */ 
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

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous distributions (CONT).
*/

UNUR_DISTR *unur_distr_cont_new( void );
/* 
   Create a new (empty) object for univariate continuous distribution.
*/

/* 
   Essential parameters.
*/

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

   Due to the fact that some of the methods do not require a
   normalized PDF the following is important:

   (*) the given CDF must be the cumulative distribution function of
   the (non-truncated) distribution. If a distribution from the 
   (=>) UNURAN library of standard distributions is truncated,
   there is no need to change the CDF.

   (*) If both the CDF and PDF are used (for a method or for order
   statistics), the PDF must be the derivative of the CDF.
   If a truncated distribution for one of the standard distributions
   from the (=>) UNURAN library of standard distributions is used,
   there is no need to change the PDF.

   (*) If the area below the PDF is required for a given distribution
   it must be given by the unur_distr_cont_set_pdfarea() call.
   For a truncated distribution this must be of course the integral of
   the PDF in the given truncated domain.
   For distributions from the (=>) UNURAN library of standard
   distributions this is done automatically by the
   unur_distr_cont_upd_pdfarea() call.

   It is important to note that all these functions must return a
   result for all floats @var{x}. Eg., if the domain of a given
   PDF is the interval [-1,1], than the given function must return
   @code{0.0} for all points outside this interval.

   It is not possible to change such a function. Once the PDF or
   CDF is set it cannot be overwritten. A new distribution object
   has to be used instead.
*/

UNUR_FUNCT_CONT *unur_distr_cont_get_pdf( UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_get_cdf}           */

UNUR_FUNCT_CONT *unur_distr_cont_get_dpdf( UNUR_DISTR *distribution );
/* See @code{unur_distr_cont_get_cdf}           */

UNUR_FUNCT_CONT *unur_distr_cont_get_cdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PDF, the derivative of the 
   PDF and the CDF of the distribution. The pointer is of type
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
   Evaluate the PDF, derivative of the PDF and the CDF, respectively,
   at x. 
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
   Get number of parameters of the PDF and set pointer
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

int unur_distr_cont_get_truncated( UNUR_DISTR *distribution, double *left, double *right );
/* 
   Get the left and right borders of the (truncated) domain of the
   distribution. For non-truncated distribution this call is
   equivalent to the unur_distr_cont_get_domain() call.
   If the (truncated) domain is not set explicitly 
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
   Set the area below the PDF. If @code{area} is non-positive, no
   area is set and @code{unur_errno} is set to @*
   @code{UNUR_ERR_DISTR_SET}. 
   
   For a distribution object created by the 
   (=>) UNURAN library of standard distributions you always should use
   the unur_distr_cont_upd_pdfarea(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_cont_upd_pdfarea( UNUR_DISTR *distribution );
/*
   Recompute the area below the PDF of the distribution. 
   It only works for distribution objects from the
   (=>) UNURAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   This call sets the normalization constant such that the given
   PDF is the derivative of a given c.f.d., i.e. the area is 1.
   However for truncated distribution the area smaller than 1.
*/

double unur_distr_cont_get_pdfarea( UNUR_DISTR *distribution );
/* 
   Get the area below the PDF of the distribution. If this area is
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
   @var{distr} must be a pointer to a univariate continuous
   distribution. 
   The result is of the same type as of unur_distr_cont_new() calls.
   (However it cannot be used to make an order statistics out of an
   order statistics.)

   To have a PDF for the order statistics, the given distribution
   obejct must contain a CDF and a PDF. Moreover it is assumed that
   the given PDF is the derivative of the given CDF. Otherwise the
   area below the PDF of the order statistics is not computed correctly.

   Important: There is no warning when the computed area below the PDF
   of the order statistics is wrong.
*/


UNUR_DISTR *unur_distr_corder_get_distribution( UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/



/* 
   Essential parameters.
*/

int unur_distr_corder_set_rank( UNUR_DISTR *distribution, int n, int k );
/* 
   Change sample size $var{n} and rank @var{k} of order statistics.
   In case of invalid data, no parameters are changed and @code{0} is
   returned.
   The area below the PDF can be set to that of the underlying
   distribution by a unur_distr_upd_pdfarea() call.
*/

int unur_distr_corder_get_rank( UNUR_DISTR *distribution, int *n, int *k );
/* 
   Get sample size $var{n} and rank @var{k} of order statistics.
   In case of error @code{0} is returned.
*/

/* 
   Additionally most of the set and get calls for continuous
   univariate distributions work. The most important exceptions are
   that the PDF and CDF cannot be changed and
   unur_distr_cont_upd_mode() uses in any way a (slow) numerical
   method that might fail.
*/


#define unur_distr_corder_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_pdf( UNUR_DISTR *distribution ); */

#define unur_distr_corder_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_dpdf( UNUR_DISTR *distribution ); */

#define unur_distr_corder_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_cdf( UNUR_DISTR *distribution ); */
/* 
   Get the respective pointer to the PDF, the derivative of the 
   PDF and the CDF of the distribution, respectively. The pointer is of type
   double funct(double x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
   See also unur_distr_cont_get_pdf().
   (Macro)
*/

#define unur_distr_corder_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
/*  double unur_distr_corder_eval_pdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_corder_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
/*  double unur_distr_corder_eval_dpdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_corder_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
/*  double unur_distr_corder_eval_cdf( double x, UNUR_DISTR *distribution ); */
/* 
   Evaluate the PDF, derivative of the PDF. and the CDF,
   respectively, at x.
   Notice that @code{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
   See also unur_distr_cont_eval_pdf().
   (Macro)
*/


#define unur_distr_corder_set_pdfparams(distr,params,n)  \
  unur_distr_cont_set_pdfparams((distr),(params),(n))
/*  int unur_distr_corder_set_pdfparams(UNUR_DISTR *distribution,double *params,int n_params); */
/* 
   Set array of parameters for underlying distribution.
   See unur_distr_cont_set_pdfparams() for details.
   (Macro)
*/

#define unur_distr_corder_get_pdfparams(distr,params) \
   unur_distr_cont_get_pdfparams((distr),(params))
/*  int unur_distr_corder_get_pdfparams( UNUR_DISTR *distribution, double **params ); */
/* 
   Get number of parameters of the PDF of the underlying distribution
   and set pointer @code{params} to array of parameters. 
   See unur_distr_cont_get_pdfparams() for details.
   (Macro)
*/


#define unur_distr_corder_set_domain(distr,left,right) \
   unur_distr_cont_set_domain((distr),(left),(right))
/*  int unur_distr_corder_set_domain( UNUR_DISTR *distribution, double left, double right ); */
/* 
   Set the left and right borders of the domain of the
   distribution. 
   See unur_distr_cont_set_domain() for details.
   (Macro)
*/

#define unur_distr_corder_get_domain(distr,left,right) \
   unur_distr_cont_get_domain((distr),(left),(right))
/*  int unur_distr_corder_get_domain( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the domain of the
   distribution. 
   See unur_distr_cont_get_domain() for details.
   (Macro)
*/


#define unur_distr_corder_get_truncated(distr,left,right) \
   unur_distr_cont_get_truncated((distr),(left),(right))
/*  int unur_distr_corder_get_truncated( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the (truncated) domain of the
   distribution.
   See unur_distr_cont_get_truncated() for details.
   (Macro)
*/

/* 
   Derived parameters.
*/   
/*   
   The following paramters MUST be set whenever one of the essential
   parameters have been set or changed (and the parameter is required
   for the chosen method).
*/

#define unur_distr_corder_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
/*  int unur_distr_corder_set_mode( UNUR_DISTR *distribution, double mode ); */
/* 
   Set mode of distribution. 
   See also unur_distr_corder_set_mode().
   (Macro)
*/

#define unur_distr_corder_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
/*  double unur_distr_corder_upd_mode( UNUR_DISTR *distribution ); */
/* 
   Recompute the mode of the distribution numerically. Notice that
   this routine is slow might not work properly in every case.
   See also unur_distr_cont_upd_mode() for further details.
   (Macro)
*/

#define unur_distr_corder_get_mode(distr)   unur_distr_cont_get_mode((distr))
/*  double unur_distr_corder_get_mode( UNUR_DISTR *distribution ); */
/* 
   Get mode of distribution.
   See unur_distr_cont_get_mode() for details.
   (Macro)
*/


#define unur_distr_corder_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
/*  int unur_distr_corder_set_pdfarea( UNUR_DISTR *distribution, double area ); */
/* 
   Set the area below the PDF.
   See unur_distr_cont_set_pdfarea() for details.
   (Macro)
*/


#define unur_distr_corder_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
/*  double unur_distr_corder_upd_pdfarea( UNUR_DISTR *distribution ); */
/*
   Recompute the area below the PDF of the distribution. 
   It only works for order statistics for distribution objects from
   the (=>) UNURAN library of standard distributions when the
   corresponding function is available.
   unur_disr_cont_upd_pdfares() assumes that the PDF of the underlying
   distribution is normalized, i.e. it is the derivative its CDF.
   Otherwise the computed area is wrong and there is NO warning
   about this failure.
   See unur_distr_cont_upd_pdfarea() for further details.
   (Macro)
*/

#define unur_distr_corder_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
/*  double unur_distr_corder_get_pdfarea( UNUR_DISTR *distribution ); */
/* 
   Get the area below the PDF of the distribution.
   See unur_distr_cont_get_pdfarea() for details.
   (Macro)
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

   It is not necessary that the given PDF is normalized and can be
   any (positive) multiple of the PDF, i.e. the integral need not
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
   The gradiant of the PDF is stored in the array @var{result}.
   The function should return @code{0} in case of an error and must
   return a non-zero value otherwise.

   The given function must be proved the gradiant of the function
   given by a unur_distr_cvec_set_pdf() call.
*/

UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PDF of the distribution. The
   pointer is of type 
   double funct(double *x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the gradiant of the PDF of the
   distribution. The pointer is of type 
   int double funct(double *result, double *x, UNUR_DISTR *distr).
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cvec_eval_pdf( double *x, UNUR_DISTR *distribution );
/* 
   Evaluate the PDF of the distribution at @var{x}.
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
   Evaluate the gradiant of the PDF of the distribution at @var{x}.
   The result is stored in the double array @var{result}.
   Both @var{result} and @var{x} must be pointer to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   distribution, @code{0} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA} (@var{result} is left unmodified).
*/

int unur_distr_cvec_set_mean( UNUR_DISTR *distribution, double *mean );
/* 
   Set mean vector for multivariate distribution.
   @var{mean} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension returned by unur_distr_get_dim().
   A @code{NULL} pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
*/

double *unur_distr_cvec_get_mean( UNUR_DISTR *distribution );
/* 
   Get the mean vector of the distribution. The function returns a
   pointer to an array of size @code{dim}.
   If the mean vector is not marked as known the @code{NULL} pointer is
   returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   (However note that the @code{NULL} pointer also indicates the zero
   vector to avoid unnecessary computations. But then
   @code{unur_errno} is not set.)
   
   Warning: Do not modify the array that holds the mean vector!
*/

int unur_distr_cvec_set_covar( UNUR_DISTR *distribution, double *covar );
/* 
   Set covariance matrix for multivariate distribution.
   @var{covar} must be a pointer to an array of size
   @code{dim}x@code{dim}, where @code{dim} is the dimension returned
   by unur_distr_get_dim(). The rows of the matrix have to be stored
   consecutively in this array.

   @var{covar} must be a variance-covariance matrix of the
   distribution, i.e. it must be symmetric and positive definite and
   its diagonal entries (i.e. the variance of the components of the
   random vector) must be positive.
   There is no check for the positive definitness yet.

   A @code{NULL} pointer for @var{covar} is interpreted as the
   identity matrix.
*/

double *unur_distr_cvec_get_covar( UNUR_DISTR *distribution );
/* 
   Get covariance matrix of distribution. The function returns a
   pointer to an array of size @code{dim}x@code{dim}.
   The rows of the matrix have to be stored consecutively in this
   array.
   If the covariance matrix is not marked as known the @code{NULL}
   pointer is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   (However note that the @code{NULL} pointer also indicates the
   identity matrix to avoid unnecessary computations. But then
   @code{unur_errno} is not set.)

   Warning: Do not modify the array that holds the covariance matrix!
*/

int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, int par, double *params, int n_params );
/* 
   This function provides an interface for additional parameters for a
   multivariate distribution besides mean vector and covariance matrix
   which have their own calls.

   It sets the parameter with number @var{par}. 
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
   Get parameter of the PDF with number @var{par}.
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
   Set mode of distribution. @var{mode} must be a pointer to an array
   of the size returned by unur_distr_get_dim().
*/

double *unur_distr_cvec_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of distribution. The function returns a pointer to an array
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
   Set the volume below the PDF. If @var{vol} is non-positive, no
   volume is set and @code{unur_errno} is set to @*
   @code{UNUR_ERR_DISTR_SET}. 
*/

double unur_distr_cvec_get_pdfvol( UNUR_DISTR *distribution );
/* 
   Get the volume below the PDF of the distribution. If this volume is
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
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given,
   @code{0} is returned and @code{sample} is set to NULL.
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling empirical multivariate continuous distributions (VCEMP).
*/

UNUR_DISTR *unur_distr_cvemp_new( int dim ); 
/* 
   Create a new (empty) object for empirical multivariate continuous
   distribution. @var{dim} is the number of components of the random
   vector (i.e. its dimension). It must be at least 2; otherwise
   unur_distr_cemp_new() should be used to create an object for an
   empirical univariate distribution.
*/

/* Essential parameters */

int unur_distr_cvemp_set_data( UNUR_DISTR *distribution, double *sample, int n_sample );
/* 
   Set observed sample for empirical distribution.
   @var{sample} is an array of double arrays of size 
   @code{dim}x@var{n_sample}, where
   @code{dim} is the dimension of the distribution returned by
   unur_distr_get_dim(). 
   The data points must be stored consecutively in @var{sample}.
*/


int unur_distr_cvemp_get_data( UNUR_DISTR *distribution, double **sample );
/* 
   Get number of samples and set pointer @var{sample} to array of
   observations. If no sample has been given,
   @code{0} is returned and @var{sample} is set to NULL.
   If successful @var{sample} points to an array of length
   @code{dim}x@code{n_sample}, where
   @code{dim} is the dimension of the distribution returned by
   unur_distr_get_dim() and @code{n_sample} the return value of the
   function.
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
   Default is [0,INT_MAX].
*/

int unur_distr_discr_get_domain( UNUR_DISTR *distribution, int *left, int *right );
/* 
   Get the left and right borders of the domain of the
   distribution. If the domain is not set explicitly 
   the interval [0,INT_MAX] is assumed and returned.
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
