/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: stdr.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method STDR                               *
 *         (Simple Transformed Density Rejection with universal bounds)      *
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

UNUR_PAR *unur_stdr_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *unur_stdr_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

double unur_stdr_sample( UNUR_GEN *generator );
double unur_stdr_sample_check( UNUR_GEN *generator );
/* sample from generator                                                     */

void unur_stdr_free( UNUR_GEN *generator);
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_stdr_set_Fmode( UNUR_PAR *parameters, double Fmode );
/* set cdf at mode                                                           */

int unur_stdr_set_verify( UNUR_PAR *parameters, int verify );
/* turn verifying of algorithm while sampling on/off                         */

int unur_stdr_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
/* set flag for using universal squeeze (default: off)                       */

#define unur_stdr_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

