/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdr.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TDR                                *
 *         (Transformed Density Rejection)                                   *
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
   =METHOD  TDR   Transformed Density Rejection

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, dPDF

   =OPTIONAL mode

   =SPEED Set-up: slow, Sampling: fast

   =REF  [GWa92]  [HWa95]

   =DESCRIPTION
      TDR is an acceptance/rejection method that uses the concavity of a
      transformed density to construct hat function and squeezes
      automatically. Such PDFs are called T-concave. Currently the
      following transformations are implemented and can be selected by
      setting their @code{c}-values by a unur_tdr_set_c() call:

      @table @code
      @item c = 0
      T(x) = log(x)
      @item c = -0.5
      T(x) = -1/sqrt(x) @ @ @ @ @ (Default)
      @end table
   
      In future releases the transformations T(x) = -(x)^c will be
      available for any c with 0 > c > -1.
      Notice that if a PDF is T-concave for a c then it also T-concave
      for every c'<c. However the performance decreases when c' is
      smaller than c. For computational reasons we suggest the usage of 
      c = -0.5 (this is the default). 
      For c <= -1 is not bounded any more if the domain of the PDF is
      unbounded. But in the case of a bounded domain using method TABL is
      preferred to a TDR with c < -1 (except in a few special cases).
      
      We offer three variants of the algorithm. 

      @table @code
      @item GW
      squeezes between construction points
      @item PS
      squeezes proportional to hat function  @ @ @ @ @ (Default)
      @item IA
      same as variant PS but uses a compositon method with
      ``immediate acceptance'' in the region below the squeeze.
      @end table

      @code{GW} has a slightly faster setup but higher marginal generation
      times.
      @code{PS} is faster than @code{GW}. @code{IA} uses less uniform
      random numbers and is therefore faster than @code{PS}.
      
      There are lots of parameters for these methods, see below.
      
      It is possible to use this method for correlation induction by
      setting an auxilliary uniform random number generator via the
      unur_set_urng_aux() call. (Notice that this must be done after a
      possible unur_set_urng() call.)
      When an auxilliary generator is used then the number of
      uniform random numbers from the first URNG that are used for one
      generated random variate is constant and given in the following table:

      @table @code
      @item GW ... 2
      @item PS ... 2
      @item IA ... 1
      @end table
      
      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_tdr_set_verify() and unur_tdr_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.

      For densities with modes not close to 0 it is suggested either
      to set the mode of the distribution or to use the
      unur_tdr_set_center() call for provide some information about
      the main part of the PDF to avoid numerical problems.

      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_tdr_chg_truncated() call.

      @emph{Important:} The ratio between the area below the hat and
      the area below the squeeze changes when the sampling region is
      restricted. Especially it becomes (very) small when sampling
      from the (far) tail of the distribution. Then it is better to
      create a new generator object for the tail of the distribution
      only.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tdr_new( UNUR_DISTR* distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_tdr_set_c( UNUR_PAR *parameters, double c );
/* 
   Set parameter @var{c} for transformation T. 
   Currently only values between 0 and -0.5 are allowed.
   If @code{c} is between 0 and -0.5 it is set to -0.5.

   Default is @code{-0.5}.
*/


int unur_tdr_set_variant_gw( UNUR_PAR *parameters );
/* 
   Use original version with squeezes between construction points as
   proposed by Gilks & Wild  (1992).
*/

int unur_tdr_set_variant_ps( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function. This is faster than
   the original version.
   This is the default.
*/

int unur_tdr_set_variant_ia( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function together with a 
   composition method that required less uniform random numbers.
*/

int unur_tdr_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Change the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call. The hat function will not be changed.

   @emph{Important:}
   The ratio between the area below the hat and the area below the
   squeeze changes when the sampling region is restricted. Especially
   it becomes (very) small when sampling from the (far) tail of the
   distribution. Then it is better to create a generator object for the
   tail of distribution only.

   @emph{Important:}
   This call does not work for variant @code{IA} (immediate
   acceptance). In this case it switches to variant @code{PS}.

   @emph{Important:}
   It is not a good idea to use adaptave rejection sampling while 
   sampling from a domain that is a strict subset of the domain that
   has been used to construct the hat.
   For that reason adaptive adding of construction points is
   automatically disabled by this call.

   @emph{Important:} If the CDF of the hat is (almost) the same 
   for @var{left} and @var{right} and (almost) equal to @code{0} or
   @code{1}, then the truncated domain is not chanced and the call
   returns @code{0}.
*/


int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area below squeeze) / (area below hat).
   It must be a number between 0 and 1.
   When the ratio exceeds the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.

   Default is @code{0.99}.
*/

double unur_tdr_get_sqhratio( UNUR_GEN *generator );
/* 
   Get the current ratio (area below squeeze) / (area below hat)
   for the generator. (In case of an error 0 is returned.)
*/

int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added after the setup when the number of
   intervals suceeds @var{max_ivs}.

   Default is @code{100}.
*/

int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, double *stp );
/* 
   Set construction points for the hat function. If @var{stp} is NULL
   than a heuristic rule of thumb is used to get @var{n_stp}
   construction points. This is the default behavior. 

   The default number of construction points is 30.
*/

int unur_tdr_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (approximate mode) of the PDF.
   It is used to find construction points by means of a heuristical
   rule of thumb. If the mode is given the center is set equal to the
   mode.

   It is suggested to use this call to provide some information about
   the main part of the PDF to avoid numerical problems.

   By default the mode is used as center if available. 
   Otherwise @code{0} is used.
*/

int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* 
   Use the center as construction point. Default is TRUE.
*/

int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
/* 
   Use the (exact!) mode as construction point.
   Notice that the behavior of the algorithm is different to simply
   adding the mode in the list of construction points via a
   unur_tdr_set_cpoints() call. In the latter case the mode is treated
   just like any other point. However when @code{usemode} is TRUE, the
   tangent in the mode is always set to 0. Then the hat of the
   transformed density can never cut the x-axis which must never
   happen if c < 0, since otherwise the hat would not be bounded.

   Default is TRUE.
*/

int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is 2.
*/

int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_tdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_tdr_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   PDF is not T-concave. 

   With @var{pedantic} being TRUE, the
   sampling routine is exchanged by a routine that simply returns
   @code{UNUR_INFINITY}. Otherwise the new point is not added to the
   list of construction points. At least the hat function remains
   T-concave.

   Setting @var{pedantic} to FALSE allows sampling from a
   distribution which is ``almost'' T-concave and small errors are
   tolerated. However it might happen that the hat function cannot be
   improved significantly. When the hat functions that has been
   constructed by the unur_init() call is extremely large then it
   might happen that the generation times are extremely high
   (even hours are possible in extremely rare cases).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/





