/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ninv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NINV                               *
 *         (Numerical INVersion of cumulative distribution function)         *
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
  =METHOD NINV Numerical INVersion

NINV is the implementation of numerical inversion.
For finding the root it is possible to choose between
Newton's method and the regula falsi.
*/

/*
 =REQUIRED NINV

 cdf, (in addition pdf only in case of Newton's method) 
*/


/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_ninv_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *_unur_ninv_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

double _unur_ninv_sample_regula( UNUR_GEN *generator );
double _unur_ninv_sample_newton( UNUR_GEN *generator );
/* sample from generator                                                     */

void _unur_ninv_free( UNUR_GEN *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_ninv_use_newton( UNUR_PAR *parameters );
/* Use Newton's method                                                       */

int unur_ninv_use_regula( UNUR_PAR *parameters );
/* Use regula falsi                                                          */

int unur_ninv_set_max_iter( UNUR_PAR *parameters, int max_iter );
/* Set number of maximal iterations                                          */

int unur_ninv_set_x_resolution( UNUR_PAR *parameters, double x_resolution);
/* Set maximal relative error in x                                           */

int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
/*  Set starting points.
    If not set, suitable values are chosen automatically.                   */
/*   Newton:        left:           starting point                          */
/*   Regula falsi: left, right:       boundary of starting interval         */
/*   If left == right, UNURAN sets starting values as follows:              */
/*   Newton: left: CDF(left) = 0.5                                          */
/*   Regula falsi: left, right: CDF(left) = 0.1, CDF(right) = 0.9           */

int unur_ninv_use_table(UNUR_PAR *parameters);
/* Generates a table (100 points) with suitable starting values
   for the iteration. Depending on the uniform random number (used by
   NINV) an adequate starting value is choosen from the table.      
*/   
/*
   The table itself is generated by applying NINV to
   equidistant points (100) between 0 and 1.
 */

int unur_ninv_chg_max_iter(UNUR_GEN *gen, int max_iter);
/* 
   Change the maximum number of iterations of on inversion step.
*/

int unur_ninv_chg_x_resolution(UNUR_GEN *gen, double x_resolution);
/*
  Change the maximal relative error in x
*/

int unur_ninv_chg_start(UNUR_GEN *gen, double left, double right);
/* Change the starting points for numerical inversion. */
/* If left==right, UNURAN chooses the starting points (see the
   function @code{unur_ninv_set_start()}*/

int unur_ninv_chg_table(UNUR_GEN *gen);
/*
   Regenerates a table as described in @code{unur_ninv_use_table()}
   and uses it for further random number generations. 
*/

int unur_ninv_table_onoff(UNUR_GEN *gen, int onoff);
/*   If onoff = 1, the table will be used (a table must already exist!)*/
/*   If onoff = 0, the table won't be used. */

int unur_ninv_chg_domain(UNUR_GEN *gen, double left, double right);
/*
   Change the borders of the (truncated) distribution. Notice that
   neither the starting point(s) nor the table will be changed!
*/

int unur_ninv_chg_pdfparams(UNUR_GEN *generator, double *params, int n_params);
/*
   Change array of parameters of distribution in given generator object.
   Notice that it is not possible to change the number of parameters.
   This function only copies the given arguments into the array of
   distribution parameters.
   IMPORTANT: The given parameters are not checked against domain errors;
   in opposition to the (=>) unur_<distr>_new() call.
*/ 


/* =END */
/*---------------------------------------------------------------------------*/
