/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dau.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DAU                                *
 *         ((Discrete) Alias-Urn)                                            *
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

UNUR_PAR *unur_dau_new( UNUR_DISTR *distribution );
/* get default parameters for generator                                      */

UNUR_GEN *_unur_dau_init( UNUR_PAR *parameters );
/* initialize new generator                                                  */

int _unur_dau_sample( UNUR_GEN *generator );
/* sample from generator                                                     */

void _unur_dau_free( UNUR_GEN *generator );
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_dau_set_urnfactor( UNUR_PAR *parameters, double factor );
/* set factor for relative size of urn                                       */

/*---------------------------------------------------------------------------*/



