/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tabl.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TABL                               *
 *         (Ahren's TABLe method: piecewise constant hat)                    *
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

UNUR_PAR *unur_tabl_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *unur_tabl_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

double unur_tabl_sample( UNUR_GEN *generator );
double unur_tabl_sample_check( UNUR_GEN *generator );
/* sample from generator                                                     */

void unur_tabl_free( UNUR_GEN *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_tabl_set_nstp( UNUR_PAR *parameters, int n_stp );
/* set number of construction points for hat at initialization               */

int unur_tabl_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* set bound for ratio A(squeeze) / A(hat)                                   */

int unur_tabl_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* set maximum number of intervals                                           */

int unur_tabl_set_areafraction( UNUR_PAR *parameters, double fraction );
/* set parameter for equal area rule                                         */

int unur_tabl_set_slopes( UNUR_PAR *parameters, double *slopes, int n_slopes );
/* set slopes of p.d.f.                                                      */

int unur_tabl_set_boundary( UNUR_PAR *parameters, double left, double right );
/* set left and right boundary of computation interval                       */

int unur_tabl_set_variant( UNUR_PAR *parameters, unsigned variant );
/* set variant for method                                                    */

int unur_tabl_set_guidefactor( UNUR_PAR *parameters, double factor );
/* set factor for relative size of guide table                               */

int unur_tabl_set_verify( UNUR_PAR *parameters, int verify );
/* turn verifying of algorithm while sampling on/off                         */

#define unur_dis_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/



