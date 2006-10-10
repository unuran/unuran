/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dari.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DARI                               *
 *         (discrete automatic rejection inversion)                          *
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
   =METHOD   DARI   discrete automatic rejection inversion

   =UP  Methods_for_DISCR

   =REQUIRED T-concave PMF, mode, approximate area

   =SPEED Set-up: moderate, Sampling: fast

   =REF  [HDa96] [HLD04: Sect.10.2, Alg.10.4]

   =DESCRIPTION
      DARI is based on rejection inversion, which can be seen as an
      adaptation of transformed density rejection to discrete
      distributions. The used transformation is @unurmath{-1/\sqrt{x}}.

      DARI uses three almost optimal points for constructing the
      (continuous) hat. Rejection is then done in horizontal
      direction. Rejection inversion uses only one uniform random
      variate per trial. 

      DARI has moderate set-up times (the PMF is evaluated nine
      times), and good marginal speed, especially if an auxiliary
      array is used to store values during generation.
      
      DARI works for all @unurmath{T_{-1/2}}-concave distributions. It requires the PMF
      and the location of the mode. Moreover the approximate sum over the PMF
      is used. (If no sum is given for the distribution the algorithm
      assumes that it is approximately 1.)
      The rejection constant is bounded from above by 4 for all @i{T}-concave
      distributions.

   =HOWTOUSE
      DARI works for discrete distribution object with given PMF.   
      The sum over probabilities should be approximately
      one. Otherwise it must be set by a unur_distr_discr_set_pmfsum()
      call to its (approximate) value.

      The size of an auxiliary table can be set by unur_dari_set_tablesize().
      The expected number of evaluations can be reduced by switching
      the use of squeezes by means of unur_dari_set_squeeze().

      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object by using the
      unur_dari_chg_pmfparams() and unur_dari_chg_domain() call, respectively.
      But then unur_dari_chg_mode() and unur_dari_chg_pmfsum() have to be used
      to reset the corresponding figures whenever they were changed.
      Before sampling from the distribution again, unur_dari_reinit() must be 
      executed. (Otherwise the generator might produce garbage).
      
      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_dari_set_verify() and unur_dari_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

UNUR_PAR *unur_dari_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_dari_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters or the domain
   of the distributions has been changed (see below).
   It is faster than destroying the existing object and building
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is returned,
   in case of a failure an error code is returned.
*/

int unur_dari_set_squeeze( UNUR_PAR *parameters, int squeeze );
/* 
   Turn utilization of the squeeze of the algorithm on/off.
   This squeeze does not resamble the squeeze of the continuous TDR
   method. It was especially designed for rejection inversion.

   The squeeze is not necessary if the size of the auxiliary table is
   big enough (for the given distribution). 
   Using a squeeze is suggested to speed up the algorithm if the
   domain of the distribution is very big or if only small samples are
   produced.  

   Default: no squeeze.
*/

int unur_dari_set_tablesize( UNUR_PAR *parameters, int size );
/* 
   Set the size for the auxiliary table, that stores constants
   computed during generation. 
   If @var{size} is set to @code{0} no table is used.
   The speed-up can be impressive if the PMF is expensive to
   evaluate and the ``main part of the distribution'' is concentrated
   in an interval shorter than the size of the table.

   Default is @code{100}.
*/

int unur_dari_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
/* 
   Set factor for position of the left and right construction point,
   resp. 
   The @var{cp_factor} is used to find almost optimal construction
   points for the hat function.
   The @var{cp_factor} must be positive and should not exceed 2.
   There is no need to change this factor in almost all situations.

   Default is @code{0.664}.
*/

int unur_dari_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_dari_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition is violated for some @i{x} then @code{unur_errno}
   is set to @code{UNUR_ERR_GEN_CONDITION}. However notice that this
   might happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/

/**********************
 *  Deprecated calls  *
 **********************/

int unur_dari_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of the distribution in a given generator
   object. Notice that this call simply copies the parameters into
   the generator object. Thus if fewer parameters are provided then
   the remaining parameters are left unchanged.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.

   @emph{Important:} The given parameters are not checked against
   domain errors; in opposition to the @command{unur_<distr>_new} calls.
*/

int unur_dari_chg_domain( UNUR_GEN *generator, int left, int right );
/* 
   Change the left and right border of the domain of the 
   (truncated) distribution.  
   If the mode changes when the domain of the (truncated) distribution is 
   changed, then a correspondig unur_dari_chg_mode() call is required.
   (There is no domain checking as in the unur_init() call.)
   Use @code{INT_MIN} and @code{INT_MAX} for (minus) infinity.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_chg_mode( UNUR_GEN *generator, int mode );
/* 
   Change mode of distribution.
   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/


int unur_dari_upd_mode( UNUR_GEN *generator );
/* 
   Recompute the mode of the distribution. This call only works well
   when a distribution object from the UNURAN library of standard
   distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise a (slow) numerical mode finder is called.
   If no mode can be found, then an error code is returnded and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_chg_pmfsum( UNUR_GEN *generator, double sum );
/* 
   Change sum over the PMF of distribution.
   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

int unur_dari_upd_pmfsum( UNUR_GEN *generator );
/* 
   Recompute sum over the PMF of the distribution. 
   It only works when a distribution objects from the
   UNURAN library of standard distributions is used
   (@pxref{Stddist,,Standard distributions}).
   Otherwise an error code is returned and @code{unur_errno} 
   is set to @code{UNUR_ERR_DISTR_DATA}.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/



