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
   =METHOD  TDR   Transformed Density Rection

   TDR is an acceptance/rejection method that uses the concavity of a
   transformed density to construct hat function and squeezes
   automatically. Such p.d.f.s are called T-concave. Currently the
   following transformations are implemented and can be selected by
   setting their @code{c}-values by a unur_tdr_set_c() call:

      c = 0     ... T(x) = log(x)
      c = -0.5  ... T(x) = -1/sqrt(x)
   
   In future releases the transformations T(x) = -(x)^c will be
   available for any c with 0 > c > -1.
   Notice that if a p.d.f. is T-concave for a c then it also T-concave
   for every c'<c. However the performance decreases when c' is
   smaller than c. For computational reasons we suggest the usage of 
   c = -0.5 (this is the default). 
   For c <= -1 is not bounded any more if the domain of the p.d.f. is
   unbounded. But in the case of a bounded domain using method TABL is
   preferred to a TDR with c < -1 (except in a few special cases).

   We offer three variant of the algorithm. 

      GW  ... squeezes between construction points
      PS  ... squeezes proportional to hat function
      IA  ... same as variant PS but uses a compositon method with
              "immediate acceptance" in the region below the squeeze.

   GW has a slightly faster setup but higher marginal generation
   times.
   PS is faster than GW. IA uses less uniform random numbers is faster
   than PS.

   There are lots of parameters for this methods, see below.

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tdr_new( UNUR_DISTR* distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_tdr_set_c( UNUR_PAR *parameters, double c );
/* 
   Set parameter c for transformation T. 
   Currently only values between 0 and -0.5 are allowed.
   If @code{c} is between 0 and -0.5 it is set to -0.5.
   Default is -0.5.
*/


int unur_tdr_set_variant_gw( UNUR_PAR *parameters );
/* 
   Use original version with squeezes between construction points as
   proposed by Gilks & Wild  (1992).
   This is the default.
*/

int unur_tdr_set_variant_ps( UNUR_PAR *parameters );
/*
  Use squeezes proportional to the hat function. This is faster than the
  original version.
*/

int unur_tdr_set_variant_ia( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function together with a 
   composition method that required less uniform random numbers.
*/

int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area below squeeze) / (area below hat).
   It must be a number between 0 and 1.
   When the ratio exceed the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.
   Default is ??.
*/

int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals (default is ??).
   No construction points are added after the setup when the number of
   intervals suceeds @code{max_ivs}.
*/

int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, double *stp );
/* 
   Set construction points for the hat function. If @code{stp} is NULL
   than a heuristical rule of thumb is used to get @code{n_stp}
   construction points. This is the default behavior. The default
   number of construction points is ??.
*/

int unur_tdr_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (approximate mode) of the p.d.f.
   It is used to find construction points by means of a heuristical
   rule of thumb. If the mode is given the center is set equal to the
   mode.
*/

int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* 
   Use the center as construction point. Default is TRUE.
*/

int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
/* 
   Use the (exact!) mode as construction point. Default is TRUE.
   Notice that the behavior of the algorithm is different to simply
   adding the mode in the list of construction points via a
   unur_tdr_set_cpoints() call. In the latter case the mode is treated
   just like any other point. However when @code{usemode} is TRUE, the
   tangent in the mode is always set to 0. Then the hat of the
   transformed density can never cut the x-axis which must never
   happen if c < 0, since otherwise the hat would not be bounded.
*/

int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT). It must be greater than or equal to 0.
   If it is set to 0, then sequential search is used.
   Default is ??.
*/

int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

int unur_tdr_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   p.d.f. is not T-concave. With @code{pedantic} being TRUE, the
   sampling routine is exchanged by a routine that simply returns
   UNUR_INFINITY. Otherwise the new point is not added to the list of
   construction points. At least the hat function remains T-concave.
   Setting @code{pedantic} to FALSE allows sampling from a
   distribution which is "almost" T-concave and small errors are
   acceptable. However it might happen that the hat function cannot be
   improved significantly. Then when the hat functions that has been
   constructed by the unur_init() call is extremely large and the
   generation times is extremely high (in theory even hours to get one
   random number are possible).
   Default is TRUE.
*/

/* =END */
/*---------------------------------------------------------------------------*/



