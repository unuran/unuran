/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_dau.h                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method DAU        *
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
/* Information for constructing the generator                                */

struct unur_dau_par { 
  double  urn_factor;  /* relative length of table for alias-urn method      */
                       /*    (DEFAULT = 1 --> alias method)                  */
                       /*   length of table = urn_factor * len               */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_dau_gen { 
  int     len;         /* length of probability vector                       */
  int     urn_size;    /* size of table for alias-urn method                 */
  double *qx;          /* pointer to cut points for strips                   */
  int    *jx;          /* pointer to donor                                   */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_dau_new(  struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_dau_init( struct unur_par *parameters );
/* initialize new generator                                                  */

int unur_dau_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_dau_free( struct unur_gen *generator );
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_dau_set_urnfactor( struct unur_par *par, double factor );
/* set factor for relative size of urn                                       */

#define unur_dau_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/



