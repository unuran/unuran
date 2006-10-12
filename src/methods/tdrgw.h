/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdrgw.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TDRGW                              *
 *         (Transformed Density Rejection - Gilks & Wild variant)            *
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
   =METHOD  TDRGW   Transformed Density Rejection - Gilks&Wild variant

   =UP  Methods_for_CONT

   =REQUIRED concave logPDF, derivative of logPDF

   =OPTIONAL mode

   =SPEED Set-up: fast, Sampling: slow

   =REINIT not implemented

   =REF  [GWa92] [HLD04: Cha.4]

   =DESCRIPTION
      TDRGW is an acceptance/rejection method that uses the concavity
      of the log-density function to construct hat function and
      squeezes automatically. 
      It is very similar to method TDR (@pxref{TDR}) with variant GW,
      parameter @code{c = 0}, and DARS switched off. 
      Moreover, method TDRGW requires the logPDF and its derivative
      dlogPDF to run. On the other hand, it is designed to draw only a
      (very) small samples and it is much more robust against
      densities with very large or small areas below the PDF as 
      it occurs, for example, in conditional distributions of
      (high dimensional) multivariate distributions. 
      Additionally, it can be re-initialized when the underlying
      distribution has been modified. 
      Thus it is well suited for Gibbs sampling. 

      Notice, that method TDRGW is a restricted version of TDR. If the
      full functionally of Transformed Density Rejection is needed use
      method @ref{TDR}.

   =HOWTOUSE
      Method TDRGW is designed for distributions with log-concave
      densities. To use this method you need a distribution object
      with the logarithm of the PDF and its derivative given.

      The number of construction points as well as a set of such
      points can be provided using unur_tdrgw_set_cpoints().
      Notice that addition construction points are added by means of
      adaptive rejection sampling until the maximal number of
      intervals given by unur_tdrgw_set_max_intervals() is reached.
      
      A generated distribution object can be re-initialized using the
      unur_tdrgw_reinit() call. However, this only works if the
      generator does not use a private copy but a pointer to the
      distribution object. This feature can be switched on by means of 
      a unur_set_use_distr_privatecopy() call. Of course the
      distribution object given as argument to the unur_tdrgw_new()
      must not be freed and unur_tdrgw_reinit() must be called
      whenever this distribution object has been changed.

      When unur_tdrgw_reinit() is called construction points for the
      new generator are necessary. There are two options:
      Either the same construction points as for the initial generator 
      (given by a unur_tdrgw_set_cpoints() call) are used (this is the
      default), or percentiles of the old hat function can be used.
      This can be set or changed using unur_tdrgw_set_reinit_percentiles()
      and unur_tdrgw_chg_reinit_percentiles().
      This feature is usefull when the underlying distribution object
      is only moderately changed. (An example is Gibbs sampling with
      small correlations.)

      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_tdrgw_set_verify() and unur_tdrgw_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tdrgw_new( const UNUR_DISTR* distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_tdrgw_reinit( UNUR_GEN *generator );
/* 
   Update an existing generator object after the distribution has been
   modified. It must be executed whenever the parameters of the
   distribution has been changed. 
   It is faster than destroying the existing object and build
   a new one from scratch.
   If reinitialization has been successful @code{UNUR_SUCCESS} is
   returned. In case of a failure an error code is returned. Then
   @var{generator} cannot be used before another successful reinit
   (with proper parameters for the underlying distribution).

   Notice that currently reinit only makes sense when the
   @var{generator} object does not use a private copy of the given
   distribution object. This behavior can be switched on by means of a
   unur_set_use_distr_privatecopy() call.
   The it is possible to manipulate the (external) distribution and
   (immediately after doing this) reinitialize the generator object.
   Obviously using this procedure must be done with extreme care.
*/

/*...........................................................................*/

int unur_tdrgw_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added after the setup when the number of
   intervals suceeds @var{max_ivs}. 
   It is increased automatically to twice the number of construction
   points if this is larger.

   Default is @code{200}.
*/

int unur_tdrgw_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
/* 
   Set construction points for the hat function. If @var{cpoints} is
   NULL then a heuristic rule of thumb is used to get @var{n_cpoints}
   construction points. This is the default behavior. 
   @var{n_cpoints} should be at least @code{2}, otherwise defaults are used.

   The default number of construction points is 2.
*/

int unur_tdrgw_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
/* */ 

int unur_tdrgw_chg_reinit_percentiles( UNUR_GEN *generator, int n_percentiles, const double *percentiles );
/* 
   By default, when the @var{generator} object is reinitialized, it
   used the same construction points as for the initialization
   procedure.
   Often the underlying distribution object has been changed only
   moderately. For example, the full conditional distribution of a
   multivariate distribution. 
   In this case it might be more appropriate to use
   percentilesm of the hat function for the last (unchanged)
   distribution. @var{percentiles} must then be a pointer to an
   ordered array of numbers between @code{0.01} and @code{0.99}.
   If @var{percentiles} is NULL, then a heuristic rule of thumb is
   used to get @var{n_percentiles} values for these percentiles.
   Notice that @var{n_percentiles} must be at least @code{2},
   otherwise defaults are used.
   (Then the first and third quartiles are used by default.) 
*/

int unur_tdrgw_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
/* */ 

int unur_tdrgw_chg_reinit_ncpoints( UNUR_GEN *generator, int ncpoints );
/* 
   When reinit fails with the given construction points or the percentiles 
   of the old hat function, another trial is undertaken with @var{ncpoints}
   construction points. @var{ncpoints} must be at least @code{10}.

   Default: @code{30}
 */

int unur_tdrgw_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_tdrgw_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_tdrgw_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   PDF is not log-concave. 

   With @var{pedantic} being TRUE, the
   sampling routine is exchanged by a routine that simply returns
   @code{UNUR_INFINITY}. Otherwise the new point is not added to the
   list of construction points. At least the hat function remains
   log-concave.

   Setting @var{pedantic} to FALSE allows sampling from a
   distribution which is ``almost'' log-concave and small errors are
   tolerated. However it might happen that the hat function cannot be
   improved significantly. When the hat functions that has been
   constructed by the unur_init() call is extremely large then it
   might happen that the generation times are extremely high
   (even hours are possible in extremely rare cases).

   Default is FALSE.
*/

double unur_tdrgw_get_loghatarea( const UNUR_GEN *generator );
/* 
   Get the logarithm of area below the hat for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_tdrgw_eval_invcdfhat( const UNUR_GEN *generator, double u );
/* 
   Evaluate the inverse of the CDF of the hat distribution at @var{u}.

   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).
*/

/* =END */
/*---------------------------------------------------------------------------*/
