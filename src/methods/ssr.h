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
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, mode, area

   =SPEED Set-up: fast, Sampling: slow

   =REF  [LJa01] [HLD04: Sect.6.3.3, Alg.6.6]

   =DESCRIPTION
      SSR is an acceptance/rejection method that uses universal
      inequalities for constructing (universal) hats and squeezes
      (@pxref{Rejection}).
      It works for all @i{T}-concave distributions with 
      @unurmath{T(x) = -1/\sqrt{x}.}

      It requires the PDF, the (exact) location of the mode and the
      area below the given PDF. The rejection constant is 4 for all 
      @i{T}-concave distributions with unbounded domain and is less
      than 4 when the domain is bounded. Optionally the CDF at the
      mode can be given to increase the performance of the algorithm.
      Then the rejection constant is at most 2 and a universal squeeze
      can (but need not be) used. However, using squeezes is not
      recommended unless the evaluation of the PDF is expensive.

      The exact location of the mode and/or the area below the PDF can
      be replace by appropriate bounds. Then the algorithm still works
      but has larger rejection constants.

   =HOWTOUSE
      SSR works for any continuous univariate distribution object with
      given @i{T}-concave PDF (with @unurmath{T(x) = -1/\sqrt{x},)}
      mode and area below PDF. Optional the CDF at the mode
      can be given to increase the performance of the algorithm by
      means of the unur_ssr_set_cdfatmode() call. Additionally
      squeezes can be used and switched on via
      unur_ssr_set_usesqueeze().

      If the (exact) area below the PDF is not known, then an upper
      bound can be used instead (which of course increases the rejection
      constant).  But then the squeeze flag must not be set and
      unur_ssr_set_cdfatmode() must not be used.

      If the exact location of the mode is not known, then use the
      approximate location and provide the (exact) value of the PDF at
      the mode by means of the unur_ssr_set_pdfatmode() call. But then
      unur_ssr_set_cdfatmode() must not be used. Notice, that a (slow)
      numerical mode finder will be used if no mode is given at all.
      It is even possible to give an upper bound for the PDF only.
      However, then the (upper bound for the) area below the PDF has to be
      multiplied by the ratio between the upper bound and the lower bound of
      the PDF at the mode.  Again setting the squeeze flag and using
      unur_ssr_set_cdfatmode() is not allowed.
      
      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object using the
      unur_ssr_chg_pdfparams() and unur_ssr_chg_domain() call, respectively.
      Then unur_ssr_chg_pdfarea(), unur_ssr_chg_mode() and 
      unur_ssr_chg_cdfatmode() have to be used to reset the corresponding figures 
      whenever they have changed.
      If the PDF at the mode has been provided by a
      unur_ssr_set_pdfatmode() call, additionally
      unur_ssr_chg_pdfatmode() must be used (otherwise this call is
      not necessary since then this figure is computed directly from the PDF).

      @emph{Important:}
      If any of mode, PDF or CDF at the mode, or the area below the mode
      has been changed, then unur_ssr_reinit() must be executed.
      (Otherwise the generator produces garbage).

      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not while sampling. It can be
      switched on/off by calling unur_ssr_set_verify() and
      unur_ssr_chg_verify(), respectively.
      Notice, however, that sampling is (a little bit) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_ssr_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_ssr_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distribution has been changed (see below).
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_ssr_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* 
   Set CDF at mode. 
   When set, the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   unur_ssr_chg_cdfatmode() has to be used to update this value.

   Default: not set.
*/

int unur_ssr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* 
   Set pdf at mode.
   When set, the PDF at the mode is never changed.          
   This is to avoid additional computations, when the PDF does not
   change when parameters of the distributions vary. 
   It is only useful when the PDF at the mode does not change with
   changing parameters for the distribution.

   Default: not set.
*/

int unur_ssr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* 
   Set flag for using universal squeeze (default: off).
   Using squeezes is only useful when the evaluation of the PDF is 
   (extremely) expensive.
   Using squeezes is automatically disabled when the CDF at the mode
   is not given (then no universal squeezes exist).

   Default is FALSE.
*/

int unur_ssr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_ssr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/*...........................................................................*/

int unur_ssr_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then an error code is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then an error code is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
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
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

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
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_ssr_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/

