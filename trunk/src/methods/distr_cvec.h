/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_cvec.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CVEC  (continuous multivariate distribution)                *
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
   =NODE   CVEC   Continuous multivariate distributions

   =UP Distribution_objects [30]

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling multivariate continuous distributions (CVEC).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cvec_new( int dim );
/* 
   Create a new (empty) object for multivariate continuous
   distribution. @var{dim} is the number of components of the random
   vector (i.e. its dimension). It must be at least 2; otherwise
   unur_distr_cont_new() should be used to create an object for a
   univariate distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cvec_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *pdf );
/* 
   Set respective pointer to the PDF of the @var{distribution}.
   The type of this function must be of type 
   @code{double funct(double *x, UNUR_DISTR *distr)},
   where @var{x} must be a pointer to a double array of appropriate
   size (i.e. of the same size as given to the unur_distr_cvec_new()
   call).

   It is not necessary that the given PDF is normalized, i.e. the
   integral need not be 1. 
   Nevertheless the volume below the PDF can be provided by a
   unur_distr_cvec_set_pdfvol() call.
*/

int unur_distr_cvec_set_dpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dpdf );
/* 
   Set pointer to the gradient of the PDF. The type of this function must be
   @code{int funct(double *result, double *x, UNUR_DISTR *distr)},
   where @var{result} and @var{x} must be pointers to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).
   The gradient of the PDF is stored in the array @var{result}.
   The function should return @code{0} in case of an error and must
   return a non-zero value otherwise.

   The given function must be proved the gradient of the function
   given by a unur_distr_cvec_set_pdf() call.
*/

UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the PDF of the @var{distribution}. The
   pointer is of type 
   @code{double funct(double *x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the
   @var{distribution}, the NULL pointer is returned.
*/

UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the gradient of the PDF of the
   @var{distribution}. The pointer is of type 
   @code{int double funct(double *result, double *x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the
   @var{distribution}, the NULL pointer is returned.
*/

double unur_distr_cvec_eval_pdf( const double *x, const UNUR_DISTR *distribution );
/* 
   Evaluate the PDF of the @var{distribution} at @var{x}.
   @var{x} must be a pointers to a double arrays of appropriate size
   (i.e. of the same size as given to the unur_distr_cvec_new() call)
   that contains the vector for which the function has to be evaluated.

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   @var{distribution}, @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_eval_dpdf( double *result, const double *x, const UNUR_DISTR *distribution );
/* 
   Evaluate the gradient of the PDF of the @var{distribution} at
   @var{x}. 
   The result is stored in the double array @var{result}.
   Both @var{result} and @var{x} must be pointer to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   @var{distribution}, @code{0} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA} (@var{result} is left unmodified).
*/

int unur_distr_cvec_set_mean( UNUR_DISTR *distribution, const double *mean );
/* 
   Set mean vector for multivariate @var{distribution}.
   @var{mean} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension returned by unur_distr_get_dim().
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).
*/

const double *unur_distr_cvec_get_mean( const UNUR_DISTR *distribution );
/* 
   Get the mean vector of the @var{distribution}. The function returns a
   pointer to an array of size @code{dim}.
   If the mean vector is not marked as known the NULL pointer is
   returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   (However note that the NULL pointer also indicates the zero
   vector to avoid unnecessary computations. But then
   @code{unur_errno} is not set.)
   
   @emph{Important:} Do @strong{not} modify the array that holds the
   mean vector!
*/

int unur_distr_cvec_set_covar( UNUR_DISTR *distribution, const double *covar );
/* 
   Set covariance matrix for multivariate @var{distribution}.
   @var{covar} must be a pointer to an array of size
   @code{dim} x @code{dim}, where @code{dim} is the dimension returned
   by unur_distr_get_dim(). The rows of the matrix have to be stored
   consecutively in this array.

   @var{covar} must be a variance-covariance matrix of the
   @var{distribution}, i.e. it must be symmetric and positive definite and
   its diagonal entries (i.e. the variance of the components of the
   random vector) must be positive.
   There is no check for the positive definitnes yet.

   A NULL pointer for @var{covar} is interpreted as the
   identity matrix.
*/

const double *unur_distr_cvec_get_covar( const UNUR_DISTR *distribution );
/* 
   Get covariance matrix of @var{distribution}. The function returns a
   pointer to an array of size @code{dim} x @code{dim}.
   The rows of the matrix have to be stored consecutively in this
   array.
   If the covariance matrix is not marked as known the NULL
   pointer is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   (However note that the NULL pointer also indicates the
   identity matrix to avoid unnecessary computations. But then
   @code{unur_errno} is not set.)

   @emph{Important:} Do @strong{not} modify the array that holds the
   covariance matrix!
*/

int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, int par, const double *params, int n_params );
/* 
   This function provides an interface for additional parameters for a
   multivariate @var{distribution} besides mean vector and covariance
   matrix which have their own calls.

   It sets the parameter with number @var{par}. 
   @var{par} indicates directly which of the parameters is set and
   must be a number between @code{0} and @code{UNUR_DISTR_MAXPARAMS}-1
   (the upper limit of possible parameters defined in
   @file{unuran_config.h}; it is set to 5 but can be changed to any 
   appropriate nonnegative number.)

   The entries of a this parameter are given by the array @var{params}
   of size @var{n_params}. Notice that using this interface an
   An (@i{n} x @i{m})-matrix has to be stored in an array of length
   @var{n_params} = @i{n} times @i{m}; where the rows of the matrix
   are stored consecutively in this array.

   Due to great variety of possible parameters for a multivariate
   @var{distribution} there is no simpler interface.

   If an error occurs no parameters are copied into the parameter
   object @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_get_pdfparams( const UNUR_DISTR *distribution, int par, const double **params );
/* 
   Get parameter of the PDF with number @var{par}.
   The pointer to the parameter array is stored in @var{params}, its
   size is returned by the function.
   If the requested parameter is not set, then @code{0} is returned
   and @code{params} is set to NULL.

   @emph{Important:} Do @strong{not} change the entries in @var{params}!
*/

/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the
   essential parameters has been set or changed (and the parameter is
   required for the chosen method).
*/

int unur_distr_cvec_set_mode( UNUR_DISTR *distribution, const double *mode );
/* 
   Set mode of @var{distribution}. @var{mode} must be a pointer to an
   array of the size returned by unur_distr_get_dim().
*/

const double *unur_distr_cvec_get_mode( const UNUR_DISTR *distribution );
/* 
   Get mode of @var{distribution}. The function returns a pointer to
   an array of the size returned by unur_distr_get_dim().
   If the mode is not marked as known the NULL pointer is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}. 
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the @var{distribution} at all.)

   @emph{Important:} Do @strong{not} modify the array that holds the mode!
*/

int unur_distr_cvec_set_pdfvol( UNUR_DISTR *distribution, double volume );
/* 
   Set the volume below the PDF. If @var{vol} is non-positive, no
   volume is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 
*/

double unur_distr_cvec_get_pdfvol( const UNUR_DISTR *distribution );
/* 
   Get the volume below the PDF of the @var{distribution}. If this volume is
   not known,@* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/


