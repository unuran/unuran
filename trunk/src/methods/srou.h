/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: srou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method SROU                               *
 *         (Simple universal generator, ratio-of-uniforms method)            *
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
   =METHOD  SROU   Simple Ratio-Of-Uniforms method

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, mode, area

   =SPEED Set-up: fast, Sampling: slow

   =REF  [LJa01] [LJa02] [HLD04: Sect.6.3.1, Sect.6.4.1, Alg.6.4, Alg.6.7]

   =DESCRIPTION
      SROU is based on the ratio-of-uniforms method
      (@pxref{Ratio-of-Uniforms}) that uses universal inequalities for
      constructing a (universal) bounding rectangle. 
      It works for all @i{T}-concave distributions, including
      log-concave and @i{T}-concave distributions with 
      @unurmath{T(x) = -1/\sqrt{x}.}

      Moreover an (optional) parameter @code{r} can be given, to
      adjust the generator to the given distribution. This parameter
      is strongly related to the parameter @code{c} for transformed
      density rejection (@pxref{TDR}) via the formula
      @i{c = -r/(r+1)}. The rejection constant increases with higher
      values for @code{r}. On the other hand, the given density must
      be @unurmath{T_c}-concave for the corresponding @i{c}.
      The default setting for @code{r} is 1 which results in a very
      simple code. (For other settings, sampling uniformly from the
      acceptance region is more complicated.)

      Optionally the CDF at the mode can be given to increase the
      performance of the algorithm. Then the rejection constant is
      reduced by 1/2 and (if @code{r=1}) even a universal squeeze can
      (but need not be) used. 
      A way to increase the performance of the algorithm when the
      CDF at the mode is not provided is the usage of the mirror
      principle (only if @code{r=1}). However using squeezes and using
      the mirror principle is only recommended when the PDF is
      expensive to compute.

      The exact location of the mode and/or the area below the PDF can
      be replace by appropriate bounds. Then the algorithm still works
      but has larger rejection constants.

   =HOWTOUSE
      SSR works for any continuous univariate distribution object with
      given @unurmath{T_c}-concave PDF with @unurmath{c<1,})
      mode and area below PDF. Optional the CDF at the mode
      can be given to increase the performance of the algorithm by
      means of the unur_ssr_set_cdfatmode() call. Additionally
      squeezes can be used and switched on via
      unur_srou_set_usesqueeze() (only if @code{r=1}).
      A way to increase the performance of the algorithm when the
      CDF at the mode is not provided is the usage of the mirror
      principle which can be swithced on by means of a
      unur_srou_set_usemirror() call (only if @code{r=1}) .
      However using squeezes and using
      the mirror principle is only recommended when the PDF is
      expensive to compute.

      The parameter @code{r} can be given, to adjust the generator to
      the given distribution. This parameter is strongly related
      parameter @code{c} for transformed density rejection via the
      formula @i{c = -r/(r+1)}. 
      The parameter @code{r} can be any value larger than or equal to
      1. Values less then 1 are automatically set to 1.
      The rejection constant depends on the chosen parameter
      @code{r} but not on the particular distribution. It is 4 for
      @code{r} equal to 1 and higher for higher values of @code{r}.
      It is important to note that different algorithms for different
      values of @code{r}: If @code{r} equal to 1 this is much faster
      than the algorithm for @code{r} greater than 1.
      The default setting for @code{r} is 1.

      If the (exact) area below the PDF is not known, then an upper
      bound can be used instead (which of course increases the rejection
      constant).  But then the squeeze flag must not be set and
      unur_srou_set_cdfatmode() must not be used.

      If the exact location of the mode is not known, then use the
      approximate location and provide the (exact) value of the PDF at
      the mode by means of the unur_srou_set_pdfatmode() call. But then
      unur_srou_set_cdfatmode() must not be used. Notice, that a (slow)
      numerical mode finder will be used if no mode is given at all.
      It is even possible to give an upper bound for the PDF only.
      However, then the (upper bound for the) area below the PDF has to be
      multiplied by the ratio between the upper bound and the lower bound of
      the PDF at the mode.  Again setting the squeeze flag and using
      unur_srou_set_cdfatmode() is not allowed.
      
      
      It is possible to change the parameters and the domain of the
      chosen distribution without building a new generator object
      using the unur_srou_chg_pdfparams() and unur_srou_chg_domain()
      call, respectively. But then unur_srou_chg_pdfarea(),
      unur_srou_chg_mode() and unur_srou_chg_cdfatmode() have to be
      used to reset the corresponding figures whenever they have
      changed. If the PDF at the mode has been provided by a 
      unur_srou_set_pdfatmode() call, additionally
      unur_srou_chg_pdfatmode() must be used (otherwise this call is
      not necessary since then this figure is computed directly from
      the PDF). 

      @emph{Important:}
      If any of mode, PDF or CDF at the mode, or the area below the
      mode has been changed, then unur_srou_reinit() must be executed.
      (Otherwise the generator produces garbage).
      
      There exists a test mode that verifies whether the conditions
      for the method are satisfied or not while sampling. It can be
      switched on by calling unur_srou_set_verify() and
      unur_srou_chg_verify(), respectively. Notice however that
      sampling is (a little bit) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_srou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_srou_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions have been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_srou_set_r( UNUR_PAR *parameters, double r );
/* 
   Set parameter @var{r} for transformation.
   Only values greater than or equal to 1 are allowed.
   The performance of the generator decreases when @var{r} is
   increased. On the other hand @var{r} must not be set to small,
   since the given density must be T_c-concave for 
   @i{c = -r/(r+1)}.

   @emph{Notice:} If @var{r} is set to @code{1} a simpler and much
   faster algorithm is used then for @var{r} greater than one.
   
   For computational reasons values of @var{r} that are greater than
   @code{1} but less than @code{1.01} are always set to @code{1.01}.

   Default is @code{1}.
*/

int unur_srou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* 
   Set CDF at mode. 
   When set, the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   unur_srou_chg_cdfatmode() has to be used to update this value.

   Default: not set.
*/

int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* 
   Set pdf at mode. 
   When set, the PDF at the mode is never changed.          
   This is to avoid additional computations, when the PDF does not
   change when parameters of the distributions vary. 
   It is only useful when the PDF at the mode does not change with
   changing parameters of the distribution.

   @emph{IMPORTANT:}
   This call has to be executed after a possible call of 
   unur_srou_set_r().

   Default: not set.
*/

int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* 
   Set flag for using universal squeeze (default: off).
   Using squeezes is only useful when the evaluation of the PDF is 
   (extremely) expensive.
   Using squeezes is automatically disabled when the CDF at the mode
   is not given (then no universal squeezes exist).

   Default is FALSE.
*/

int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
/* 
   Set flag for using mirror principle (default: off).
   Using the mirror principle is only useful when the CDF at the
   mode is not known and the evaluation of the PDF is rather cheap compared
   to the marginal generation time of the underlying uniform random
   number generator.
   It is automatically disabled when the CDF at the mode is given.
   (Then there is no necessity to use the mirror principle. However disabling
   is only done during the initialization step but not at a re-initialization
   step.)

   Default is FALSE.
*/

int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_srou_chg_verify( UNUR_GEN *generator, int verify );
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

int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
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

int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
/* 
   Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_srou_chg_mode() is required.
   (There is no checking whether the domain is set or not as in the
   unur_init() call.)
*/

int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
/* 
   Change mode of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
/* 
   Change CDF at mode of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_chg_pdfatmode( UNUR_GEN *generator, double fmode );
/* 
   Change PDF at mode of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
/* 
   Change area below PDF of distribution.
   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_srou_upd_pdfarea( UNUR_GEN *generator );
/*
   Recompute the area below the PDF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_srou_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/


