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

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

UNUR_PAR *unur_arou_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *_unur_arou_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

double _unur_arou_sample( UNUR_GEN *generator );
double _unur_arou_sample_check( UNUR_GEN *generator );
/* sample from generator                                                     */

void _unur_arou_free( UNUR_GEN *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, double *stp );
/* set construction points for envelope and/or its number for initialization */

int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
/* set factor for relative size of guide table                               */

int unur_arou_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* set bound for ratio A(squeeze) / A(hat)                                   */

int unur_arou_set_max_segments( UNUR_PAR *parameters, int max_segs );
/* set maximum number of segments                                            */

int unur_arou_set_center( UNUR_PAR *parameters, double center );
/* set center (approximate mode) of p.d.f.                                   */

int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* set flag for using center as construction point                           */

int unur_arou_set_verify( UNUR_PAR *parameters, int verify );
/* turn verifying of algorithm while sampling on/off                         */

/*---------------------------------------------------------------------------*/



