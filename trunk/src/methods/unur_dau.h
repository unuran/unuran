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
 *         only used in unur_methods.h                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_dau_par { 
  double *prob;        /* pointer to probability vector    (DEFAULT = NULL)  */
  int     len;         /* length of probability vector     (DEFAULT = 0)     */
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

  double *prob;        /* pointer to probability vector    (DEFAULT = NULL)  */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_dau_new( double *probabilities, int plength );
/* get default parameters for generator                                      */

struct unur_gen *unur_dau_init( struct unur_par *parameters );
/* initialize new generator                                                  */

int unur_dau_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_dau_free( struct unur_gen *generator );
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/



