/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hitro_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method HITRO                              *
 *         (Markov Chain - HIT-and-run sampler with Ratio-Of-uniforms)       *
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

struct unur_hitro_par { 
  double r;                  /* r-parameter of HITRO method                  */
  int thinning;              /* thinning factor for generated chain          */
  int burnin;                /* length of burn-in for chain                  */
  double vmax;               /* bounding rectangle v-coordinate              */
  const double *umin, *umax; /* bounding rectangle u-coordinates             */
  const double *x0;          /* starting point of chain                      */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_hitro_gen {
  int dim;                   /* dimension of distribution                    */
  int thinning;              /* thinning factor for generated chain          */
  double r;                  /* r-parameter of HITRO method                  */

  double *state;             /* state of chain / current point.
				the state is a point in the acceptance region
				of the RoU-method in vu-hyperplane.
				the coordinates are stored in the following order:
				state = {v, u_1, u_2, ... u_n}               */
        
  double vmax;               /* bounding rectangle v-coordinate              */
  double *umin, *umax;       /* bounding rectangle u-coordinates             */
  const double *center;      /* center of distribution                       */

  int    coord;              /* current coordinate used for HITRO chain      */
  double *direction;         /* working array for random direction           */

  int burnin;                /* length of burn-in for chain                  */
  double *x0;                /* starting point of chain                      */
};

/*---------------------------------------------------------------------------*/

/* struct unur_hitrou_gen { */
/*   double *umin, *umax;      /\* boundary rectangle u-coordinates              *\/ */
/*   double vmax;              /\* boundary rectangle/strip v-coordinate         *\/ */
/*   const double *center;     /\* center of distribution                        *\/ */
/*   double *direction;        /\* random direction vector                       *\/ */
/*   double *point_current;    /\* current point inside the shape                *\/ */
/*   double *point_random;     /\* random point, can be inside shape or not      *\/ */
/*   double *x;                /\* working point in the (xy)-coordinate system   *\/ */
/*   int coordinate;           /\* current coordinate used for stepping          *\/ */
/* }; */
