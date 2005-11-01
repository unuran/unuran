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
   =METHOD  TDRGW   Transformed Density Rejection - Gild&Wild variant

   =UP  Methods_for_CONT

   =REQUIRED concave logPDF, derivative of logPDF

   =OPTIONAL mode

   =SPEED Set-up: fast, Sampling: slow

   =REF  [GWa92] [HLD04: Cha.4]

   =DESCRIPTION
      TDRGW is an acceptance/rejection method that uses the concavity
      of the log-density function to construct hat function and
      squeezes automatically. 
      It is very similar to method TDR (@pxref{TDR}) with variant GW,
      parameter @code{c = 0}, and DARS switched off. 
      Moreover, method TDRGW requires the logPDF and its derivative
      dlogPDF to run. On the other hand, it is much more robust
      against densities with very large or small areas below the PDF as
      it occurs for example in conditional distributions of
      (high dimensional) multivariate distributions. Thus it is well
      suited for Gibbs sampling. 

      Notice, that method TDRGW is a restricted version of TDR. If the
      full functionally of Transformed Density Rejection is needed use
      method @ref{TDR}.

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

double unur_tdrgw_get_hatarea( const UNUR_GEN *generator );
/* 
   Get the area below the hat for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

int unur_tdrgw_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added after the setup when the number of
   intervals suceeds @var{max_ivs}. 
   It is increased automatically to twice the number of construction
   points if this is larger.

   Default is @code{100}.
*/

int unur_tdrgw_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
/* 
   Set construction points for the hat function. If @var{cpoints} is
   NULL than a heuristic rule of thumb is used to get @var{n_cpoints}
   construction points. This is the default behavior. 

   The default number of construction points is 2.
*/

int unur_tdrgw_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is 2.
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

/* =END */
/*---------------------------------------------------------------------------*/
