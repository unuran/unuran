/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_cont.h                                                      *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CONT  (continuous univariate distribution)                  *
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

