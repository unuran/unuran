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

int unur_ninv_set_start( UNUR_PAR *parameters, double s1, double s2, double s3 );
/*  Set starting points.
    If not set, suitable values are chosen automatically.                   */
/*   Newton:        s1:           starting point                            */
/*   Regula falsi: s1, s2:       boundary of starting interval              */
/*   If s1 == s2, UNURAN sets starting values as follows:                   */
/*   Newton: CDF(s1) = 0.5                                                  */
/*   Regula falsi: CDF(s1) = 0.1, CDF(s2) = .9                              */

/* =END */
/*---------------------------------------------------------------------------*/










