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
/*
  Initialize a generator object. All necessary information must be
  stored in the parameter object.
  Warning: If an error has occurred a NULL pointer is return. This
  must not be used for the sampling routines (this causes a
  segmentation fault). Thus always check if the call was successful.
*/

int unur_reinit( UNUR_GEN *generator );
/*
  Update an existing generator object after the distribution has been
  changed. It is faster than destroying the existing object and build
  a new one from scratch.
  However it is not implemented all methods yet.
  (was passiert wenn ein fehler auftritt??) return value = 0 ??
*/

/*---------------------------------------------------------------------------*/
/* Sample from generator                                                     */

int    unur_sample_discr(UNUR_GEN *generator);
double unur_sample_cont(UNUR_GEN *generator);
void   unur_sample_vec(UNUR_GEN *generator, double *vector);
/*
  Sample from generator object. The three routines depend on the type
  of the generator object (discrete or continuous univariate
  distribution, or multivariate distribution).
  Warning: These routines do not check if generator is an invalid NULL
  pointer.
*/


/*---------------------------------------------------------------------------*/
/* Destroy (free) generator object                                           */

void   unur_free( UNUR_GEN *gen );
/*
  Destroy (free) the given generator object.
*/

/*---------------------------------------------------------------------------*/
/* Get dimension of generator for (multivariate) distribution                */

int unur_get_dimension( UNUR_GEN *generator );
/*
  Get the number of dimension of a (multivariate) distribution.
  For a univariate distribution 1 is return.
*/

/*---------------------------------------------------------------------------*/

const char *unur_get_genid( UNUR_GEN *generator );
/*
  Get identifier string for generator.
  If UNUR_ENABLE_GENID is not defined in unuran_config.h then only the method
  used for the generator is return.
*/

/*---------------------------------------------------------------------------*/

UNUR_DISTR *unur_get_distr( UNUR_GEN *generator );
/* 
   Get pointer to distribution object from generator object.
*/

/*---------------------------------------------------------------------------*/





