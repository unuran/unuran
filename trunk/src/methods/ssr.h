/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ssr.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method SSR                                *
 *         (Simple Setup, Rejection with universal bounds)                   *
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
   =METHOD  SSR   Simple Setup Rejection

   =TYPE
      CONT  continuous univariate

   =DESCRIPTION
      SSR is an acceptance/rejection method that uses universal
      inequalities for constructing (universal) hats and squeezes.
      It works for all T-concave distributions with T(x) = -1/sqrt(x).
      
      It requires the PDF, the (exact) location of the mode and the
      area below the given PDF. The rejection constant is 4 for all
      T-concave distributions with unbounded domain and is less than 4
      when the domain is bounded.  Optionally the CDF at the mode
      can be given to increase the performance of the algorithm by
      means of the unur_ssr_set_cdfatmode() call.  Then the rejection
      constant is reduced by one half and even a universal squeeze can
      (but need not be) used. However using squeezes is not
      recommended unless the evaluation of the PDF is rather expensive.
      
      (The mirror principle is not implemented.)
      
      If the exact location of the mode is not known, then use the
      approximate location and provide the (exact) value of the PDF at
      the mode by means of the unur_ssr_set_pdfatmode() call.  But then
      unur_ssr_set_cdfatmode() must not be used.  Notice if no mode is given
      at all, a (slow) numerical mode finder will be used.
      
      If the (exact) area below the PDF is not known, then an upper
      bound can be used instead (which of course increases the rejection
      constant).  But then the squeeze flag must not be set and
      unur_ssr_set_cdfatmode() must not be used.

      It is even possible to give an upper bound for the PDF only.
      However then the (upper bound for the) area below the PDF has to be
      multiplied by the ratio between the upper bound and the lower bound of
      the PDF at the mode.  Again seting the squeeze flag and using
      unur_ssr_set_cdfatmode() is not allowed.
      
      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object using the
      unur_ssr_chg_pdfparams() and unur_ssr_chg_domain() call, respectively.
      But then unur_ssr_chg_pdfarea(), unur_ssr_chg_mode() and 
      unur_ssr_chg_cdfatmode() have to used to reset the corresponding figures 
      whenever they have changed.
      If the PDF at the mode has been provided by a 
      unur_ssr_set_pdfatmode() call, additionally unur_ssr_chg_pdfatmode() must 
      be used (otherwise this call is not necessary since then this figure is
      computed directly from the PDF).
      If any of mode, PDF or CDF at the mode, or the area below the mode
      has been changed, then unur_ssr_reinit() must be executed.
      (Otherwise the generator produces garbage).

      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_ssr_set_verify() and unur_ssr_chg_verify(), respectively.
      Notice however that sampling is (a little bit) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_ssr_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_ssr_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions has been changed (see below).
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{1} is returned,
   in case of a failure @code{0} is returned.
*/

int unur_ssr_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* 
   Set CDF at mode. 
   When set the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   (=>) unur_ssr_chg_cdfatmode() has to be used to update this value.
*/

int unur_ssr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* 
   Set pdf at mode. if set the PDF at the mode is never changed.          
   This is to avoid additional computations, when the PDF does not
   change when parameters of the distributions vary. 
   It is only useful when the PDF at the mode does not change with
   changing parameters for the distribution.
*/

int unur_ssr_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

int unur_ssr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

int unur_ssr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* 
   Set flag for using universal squeeze (default: off).
   using squeezes is only useful when the evaluation of the PDF is 
   (extremely) expensive.
   Using squeezes is automatically disabled when the CDF at the mode
   is not given (then no universal squeezes exist).
*/

/*...........................................................................*/

int unur_ssr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.
   IMPORTANT: The given parameters are not checked against domain errors;
   in opposition to the (=>) unur_<distr>_new() call.
*/

int unur_ssr_chg_domain( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the distribution.  
   If the mode changes when the domain of the distribution is 
   changed, then a correspondig unur_ssr_chg_mode() is required.
   (There is no domain checking as in the unur_init() call.)
*/

int unur_ssr_chg_mode( UNUR_GEN *generator, double mode );
/* 
   Change mode of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. This call only works when
   a distribution object from the (=>) UNURAN library of standard
   distributions is used.
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.

   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
/* 
   Change CDF at mode of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
/* 
   Change PDF at mode of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_chg_pdfarea( UNUR_GEN *generator, double area );
/* 
   Change area below PDF of distribution.
   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_ssr_upd_pdfarea( UNUR_GEN *generator );
/*
   Recompute the area below the PDF of the distribution. 
   It only works when a distribution objects from the
   (=>) UNURAN library of standard distributions is used. 
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/

