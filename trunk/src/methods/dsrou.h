/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dsrou.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DSROU                              *
 *         (Discrete, Simple universal generator, Ratio-Of-Uniforms method)  *
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
   =METHOD  DSROU   Discrete Simple Ratio-Of-Uniforms method

   =UP  Methods_for_DISCR

   =REQUIRED T-concave PMF, mode, sum over PMF

   =SPEED Set-up: fast, Sampling: slow

   =REF  [LJa01]

   =DESCRIPTION
      DSROU is based on the ratio-of-uniforms method but uses universal 
      inequalities for constructing a (universal) bounding rectangle.
      It works for all T-concave distributions with T(x) = -1/sqrt(x).
      
      It requires the PMF, the (exact) location of the mode and the
      sum over the given PDF. The rejection constant is 4 for all
      T-concave distributions. Optionally the CDF at mode-1 can
      be given to increase the performance of the algorithm by means
      of the unur_dsrou_set_cdfbeforemode() call. Then the rejection
      constant is reduced to 2.
      
      If the (exact) sum over the PMF is not known, then an upper
      bound can be used instead (which of course increases the
      rejection constant). But then unur_dsrou_set_cdfbeforemode()
      must not be called.
      
      It is possible to change the parameters and the domain of the
      chosen distribution without building a new generator object
      using the unur_dsrou_chg_pmfparams() and unur_dsrou_chg_domain()
      call, respectively. But then unur_dsrou_chg_pmfsum(),
      unur_dsrou_chg_mode() and unur_dsrou_chg_cdfbeforemode() have to
      be used to reset the corresponding figures whenever they have
      changed. 

      If any of mode, CDF at mode-1, or the sum over the PMF has been
      changed, then unur_dsrou_reinit() must be executed. 
      (Otherwise the generator produces garbage).

      There exists a test mode that verifies whether the conditions
      for the method are satisfied or not while sampling. It can be
      switched on or off by calling unur_dsrou_set_verify() and
      unur_dsrou_chg_verify(), respectively.
      Notice however that sampling is (a little bit) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dsrou_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_dsrou_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distribution have been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{1} is returned,
   in case of a failure @code{0} is returned.
*/

int unur_dsrou_set_cdfbeforemode( UNUR_PAR *parameters, double Fbmode );
/* 
   Set CDF at mode-1. 
   When set, the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   unur_dsrou_chg_cdfbeforemode() has to be used to update this value.

   Default: not set.
*/

int unur_dsrou_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_dsrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PMF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/*...........................................................................*/

int unur_dsrou_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. 

   For standard distributions from the UNURAN library the parameters
   are checked. It these are invalid, then @code{0} is
   returned. Moreover the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   Notice that optional parameters are (re-)set to their default
   values if not given for UNURAN standard distributions.

   For other distributions @var{params} is simply copied into to
   distribution object. It is only checked that @var{n_params} does
   not exceed the maximum number of parameters allowed.
   Then @code{0} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.
*/

int unur_dsrou_chg_domain( UNUR_GEN *generator, int left, int right );
/* 
   Change left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_dsrou_chg_mode() is required.
   (There is no checking whether the domain is set or not as in the
   unur_init() call.)
*/

int unur_dsrou_chg_mode( UNUR_GEN *generator, int mode );
/* 
   Change mode of distribution.
   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. 
   See unur_distr_cont_upd_mode() for more details.

   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_chg_cdfbeforemode( UNUR_GEN *generator, double Fbmode );
/* 
   Change CDF at mode-1 of distribution.
   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/


int unur_dsrou_chg_pmfsum( UNUR_GEN *generator, double sum );
/* 
   Change sum over PMF of distribution.
   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dsrou_upd_pmfsum( UNUR_GEN *generator );
/*
   Recompute the sum over the the PMF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,Standard distributions,Standard distributions}).
   Otherwise @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   unur_dsrou_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/


