/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: arou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method AROU                               *
 *         (Adaptive Ratio-Of-Uniforms)                                      *
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
   =METHOD  AROU   Automatic Ratio-Of-Uniforms method

   =UP  Methods_for_CONT

   =DESCRIPTION
      AROU is a variant of the ratio-of-uniforms method that uses the
      fact that the transformed region is convex for many distributions.
      It works for all T-concave distributions with T(x) = -1/sqrt(x).
      
      There are lots of parameters for this methods, see below.
      
      It is possible to use this method for correlation induction by
      setting an auxilliary uniform random number generator via the
      unur_set_urng_aux() call. (Notice that this must be done after a
      possible unur_set_urng() call.)
      When an auxilliary generator is used then the number of used
      uniform random numbers that is used up for one generated random
      variate is constant and equal to 1.
      
      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_arou_set_verify() and unur_arou_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_arou_new( UNUR_DISTR *distribution );
/*    Get default parameters for generator. */

/*...........................................................................*/

int unur_arou_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area inside squeeze) / (area inside envelope).
   It must be a number between 0 and 1.
   When the ratio exceed the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.
   Default is ??.
*/

double unur_arou_get_sqhratio( UNUR_GEN *generator );
/* 
   Get the current ratio (area inside squeeze) / (area inside envelope)
   for the generator. (In case of error 0 is returned.)
*/


int unur_arou_set_max_segments( UNUR_PAR *parameters, int max_segs );
/* 
   Set maximum number of segements (default is ??).
   No construction points are added after the setup when the number of
   intervals suceeds @code{max_segs}.
*/


int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, double *stp );
/* 
   Set construction points for enveloping polygon. If @code{stp} is NULL
   than a heuristical rule of thumb is used to get @code{n_stp}
   construction points. This is the default behavior. The default
   number of construction points is ??.
*/


int unur_arou_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (approximate mode) of the PDF.
   It is used to find construction points by means of a heuristical
   rule of thumb. If the mode is given the center is set equal to the
   mode.
*/


int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* 
   Use the center as construction point. Default is TRUE.
*/


int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT). It must be greater than or equal to 0.
   If it is set to 0, then sequential search is used.
   Default is ??.
*/


int unur_arou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/

int unur_arou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
*/


int unur_arou_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   PDF is not T-concave. With @code{pedantic} being TRUE, the
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



