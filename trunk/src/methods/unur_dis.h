/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_dis.h                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures and function prototypes for method DIS        *
 *         ((Discrete) Alias-Urn)                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only used in unur_methods.h                                       *
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

struct unur_dis_par { 
  double *prob;         /* pointer to probability vector    (DEFAULT = NULL) */
  int     len;          /* length of probability vector     (DEFAULT = 0)    */
  double  guide_factor; /* relative length of guide table.  (DEFAULT = 1)    */
                        /*   length of guide table = guide_factor * len      */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_dis_gen { 
  double  sum;          /* sum of all probabilities = cumprob[len-1]         */
  double *cumprob;      /* pointer to the vector of cumulated probabilities  */
  int    *guide_table;  /* pointer to guide table                            */
  int     guide_size;   /* length of guide table                             */

  double *prob;         /* pointer to probability vector    (DEFAULT = NULL) */
  int     len;          /* length of probability vector     (DEFAULT = 0)    */
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_dis_new( double *probabilities, int plength );
/* get default parameters for generator                                      */

struct unur_gen *unur_dis_init( struct unur_par *parameter );
/* initialize new generator                                                  */

int unur_dis_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_dis_free( struct unur_gen *generator );
/* destroy generator object                                                  */

/*---------------------------------------------------------------------------*/

