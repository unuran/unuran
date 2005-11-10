/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mcgibbs_struct.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method MCGIBBS                            *
 *         (Markov Chain - GIBBS sampler)                                    *
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

struct unur_mcgibbs_par { 
  int thinning;              /* thinning factor for generated chain          */
  const double *x0;          /* starting point of chain                      */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_mcgibbs_gen {
  int dim;                   /* dimension of distribution                    */
  int thinning;              /* thinning factor for generated chain          */

  double *state;             /* current point                                */

  struct unur_distr *distr_condi; /* conditional distribution                */
  struct unur_gen **gen_condi; /* generators for conditional distributions   */

  int    coord;             /* current coordinate used for GIBBS chain       */
  double *direction;        /* working array for random direction            */

/*   double *direction;        /\* random direction                              *\/ */
/*   double *tdr_points;       /\* starting points for the tdr method in each dim*\/ */
/*   int    coordinate;        /\* current coordinate used for stepping          *\/ */
/*   long   pdfcount;          /\* counting the number of PDF calls              *\/ */
/*   struct unur_par   *par;   /\* pointer to parameter object for gibbs sampler *\/ */
/*   struct unur_par   *par_conditional;   /\* pointer to conditional par.   object *\/ */
/*   struct unur_distr *distr_conditional; /\* pointer to conditional distr. object *\/ */
/*   struct unur_gen   *gen_conditional;   /\* pointer to conditional gen.   object *\/ */


/* replaced */
/*   double *point_current;    --> state */

};

/*---------------------------------------------------------------------------*/

