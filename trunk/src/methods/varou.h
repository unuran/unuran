/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: varou.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VAROU                              *
 *         (Vector Adaptive Ratio-Of-Uniforms method)                           *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 
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
   =METHOD  VAROU   Multivariate Adaptive Ratio-Of-Uniforms method

   =UP  Methods_for_CVEC

   =REQUIRED PDF, dPDF, ...

   =OPTIONAL ...

   =REF  [WGS91]

   =SPEED Set-up: slow, Sampling: slow

   =DESCRIPTION
      VAROU is an implementation of the multivariate ratio-of-uniforms
      method which uses a piecewice continuous hut consisting of a family
      of conii with tops at the origin of the @unurmath{(U_1,U_2,...,U_d,V)}
      coordinate system and spanning vectors pointing towards an initial 
      triangulation of the upper half unit-sphere.
      
      The base of such conii is given as the tangent-plane to the
      surface :
     
      @display
      @unurmath{ \{ (U_1, ..., U_d, V) : V=PDF(U_1/V,...,U_d/V)^{1/(d+1)} \} }
      @end display
     
      through one of the pierce-points of the spanninge vectors.
      In particular, we choose the point 'nearest' to 
      @unurmath{(0,...,0,V_{max})} as the construction point for the
      tangent-plane.
      
      VAROU is based on the rejection method (@pxref{Rejection}).
      And it is important to note that the acceptance probability
      decreases exponentially with dimension. Thus even for moderately
      many dimensions (e.g. 5) the number of repetitions to get one
      random vector can be prohibitively large and the algorithm seems
      to stay in an infinite loop.

   =HOWTOUSE
      For using the VAROU method UNURAN needs the PDF of the
      distribution. 
      
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_varou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_varou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_varou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
