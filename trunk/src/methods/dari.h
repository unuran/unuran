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
   =METHOD   DARI   discrete automatic rejection inversion

   =UP  Methods_for_DISCR

   =DESCRIPTION
      DARI is based on rejection inversion, which can be seen as an
      adaptation of transformed density rejection to discrete
      distributions. The used transformation is  -1/sqrt(x).

      DARI uses three almost optimal points for constructing the
      (continuous) hat. Rejection is then done in horizontal
      direction. Rejection inversion uses only one uniform random
      variate per trial. 

      DARI has moderate set-up times (the PMF is evaluated nine
      times), and good marginal speed, especially if an auxilliary
      array is used to store values during generation.
      
      DARI works for all T-(-1/2)-concave distributions. It requires the PMF
      and the location of the mode. Moreover the approximate sum over the PMF
      is used. (If no sum is given for the distribution the algorithm
      assumes that it is approximately 1.)
      The rejection constant is bounded from above by 4 for all T-concave
      distributions.

      It is possible to change the parameters and the domain of the chosen 
      distribution without building a new generator object by using the
      unur_dari_chg_pmfparams() and unur_dari_chg_domain() call, respectively.
      But then unur_dari_chg_mode() and unur_dari_chg_pmfsum() have to be used
      to reset the corresponding figures whenever these have changed.
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

UNUR_PAR *unur_dari_new( UNUR_DISTR *distribution );
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
   If reinitialization has been successful @code{1} is returned,
   in case of a failure @code{0} is returned.
*/

int unur_dari_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_dari_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

int unur_dari_set_squeeze( UNUR_PAR *parameters, char squeeze );
/* 
   Turn utilization of the squeeze of the algorithm on/off.
   This squeeze does not resamble the squeeze of the continuous TDR
   method. It was especially designed for rejection inversion.

   The squeeze is not necessary if the size of the auxiliary table is
   big enough (for the given distribution). 
   Using a squeeze is suggested to speed up the algorithm if the
   domain of the distribution is very big or if only small samples are
   produced.  

   Default is off.
*/

int unur_dari_set_size( UNUR_PAR *parameters, int size );
/* 
   Set the size for the auxiliary table, that stores constants
   computed during generation. 
   If @var{size} is set to @code{0} no table is used.
   The speed-up can be impressive if the PMF is expensive to
   evaluate and the ``main part of the distribution'' is concentrated
   in an interval shorter than the size of the table.

   Default is 100.
*/

int unur_dari_set_cfactor( UNUR_PAR *parameters, double cfactor );
/* 
   Set factor for position of the left and right construction point,
   resp. 
   The c_factor is used to find almost optimal construction points for
   the hat function.
   There is no need to change this factor in almost all situations.

   Default is ??.
*/

/*...........................................................................*/

int unur_dari_chg_pmfparams( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.

   @emph{IMPORTANT:} The given parameters are not checked against
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
   for a distribution object from the UNURAN library of standard
   distributions is used
   (@pxref{Stddist,Standard distributions,Standard distributions}).
   Otherwise a (slow) numerical mode finder is used.
   If no mode can be found, then @code{0} is returnded and
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
   (@pxref{Stddist,Standard distributions,Standard distributions}).
   Otherwise @code{0} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.

   unur_dari_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/


