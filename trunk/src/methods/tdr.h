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

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

UNUR_PAR *unur_tdr_new( UNUR_DISTR* distribution );
/* get default parameters for generator                                      */


/* sample from generator                                                     */

/*...........................................................................*/

int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, double *stp );
/* set construction points for envelope and/or its number for initialization */

int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
/* set factor for relative size of guide table                               */

int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* set bound for ratio A(squeeze) / A(hat)                                   */

int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* set maximum number of intervals                                           */

int unur_tdr_set_center( UNUR_PAR *parameters, double center );
/* set center (approximate mode) of p.d.f.                                   */

int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* set flag for using center as construction point                           */

int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
/* set flag for using (exact) mode as construction point                     */

int unur_tdr_set_version_orig( UNUR_PAR *parameters );
/* 
   Use original version with squeezes as proposed by Gilks & Wild 
   (but with generalized transformations T).
   This is the default.
*/

int unur_tdr_set_version_ps( UNUR_PAR *parameters );
/*
  Use squeezes proportional to the hat function. This is faster than the
  original version.
*/

int unur_tdr_set_version_ia( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function together with a 
   composition method that required less uniform random numbers.
*/

int unur_tdr_set_c( UNUR_PAR *parameters, double c );
/* set parameter c for transformation T_c                                    */

int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
/* turn verifying of algorithm while sampling on/off                         */

/*---------------------------------------------------------------------------*/

