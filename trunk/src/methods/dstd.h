/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dstd.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DSTD                               *
 *         (wrapper for special generators for                               *
 *         Discrete STanDard distributions)                                  *
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

UNUR_PAR *unur_dstd_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *unur_dstd_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

/** 
    double unur_dstd_sample( UNUR_GEN *gen );
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

void unur_dstd_free( UNUR_GEN *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_dstd_set_variant( UNUR_PAR *parameters, unsigned variant );
/* set variant of method                                                     */

int unur_dstd_chg_param( UNUR_GEN *gen, double *params, int n_params );
/* change array of parameters for distribution                               */

#define unur_dstd_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

