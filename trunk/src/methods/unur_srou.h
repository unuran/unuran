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

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

UNUR_PAR *unur_srou_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_srou_set_Fmode( UNUR_PAR *parameters, double Fmode );
/* Set cdf at mode                                                           */

int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
/* Turn verifying of algorithm while sampling on/off                         */

int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* Set flag for using universal squeeze (default: off)                       */

int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
/* Set flag for using mirror principle (default: off)                        */

/*...........................................................................*/

int unur_srou_chg_pdfparam( UNUR_GEN *generator, double *params, int n_params );
/* 
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of 
   distribution parameters.
   IMPORTANT: The given parameters are not checked against domain errors;
   in opposition to the (=>) unur_<distr>_new() call.
   (=>) unur_reinit() must be called afterwards.
*/

int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
/* Change mode of distribution                                               */

int unur_srou_chg_Fmode( UNUR_GEN *generator, double Fmode );
/* Change c.d.f. at mode of distribution                                     */

int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
/* Change left and right border of the domain of the 
   (truncated) distribution.                                                 */

int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
/* Change area below p.d.f. of distribution                                  */

/*---------------------------------------------------------------------------*/
