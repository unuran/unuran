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
};

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

struct unur_par *unur_dis_new( struct unur_distr *distr );
/* get default parameters for generator                                      */

struct unur_gen *unur_dis_init( struct unur_par *parameter );
/* initialize new generator                                                  */

int unur_dis_sample( struct unur_gen *generator );
/* sample from generator                                                     */

void unur_dis_free( struct unur_gen *generator );
/* destroy generator object                                                  */

/*...........................................................................*/

int unur_dis_set_variant( struct unur_par *par, unsigned variant );
/* set variant of method                                                     */

int unur_dis_set_guidefactor( struct unur_par *par, double factor );
/* set factor for relative size of guide table                               */

#define unur_dis_set_debug(par,debugflags)  unur_set_debug((par),(debugflags))
/* set debuging flags                                                        */

/*---------------------------------------------------------------------------*/

