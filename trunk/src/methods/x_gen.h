/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_gen.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for handling               *
 *         generator objects.                                                *
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
/* (Re-) Initialize generators                                               */

UNUR_GEN *unur_init( UNUR_PAR *parameters );
int    unur_reinit( UNUR_GEN *generator );

/*---------------------------------------------------------------------------*/
/* Sample from generator                                                     */

int    unur_sample_discr(UNUR_GEN *generator);
double unur_sample_cont(UNUR_GEN *generator);
void   unur_sample_vec(UNUR_GEN *generator, double *vector);

/*---------------------------------------------------------------------------*/
/* Destroy (free) generator object                                           */

void   unur_free( UNUR_GEN *gen );

/*---------------------------------------------------------------------------*/
/* Get dimension of generator for (multivariate) distribution                */

int unur_get_dimension( UNUR_GEN *generator );

/*---------------------------------------------------------------------------*/
