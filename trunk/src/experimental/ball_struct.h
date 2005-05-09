/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ball_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method BALL                               *
 *         (Multivariate HIT and run Ratio Of Uniforms)                      *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
 *                                                                           *
 *****************************************************************************

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

struct unur_ball_par {
  int    dim;               /* dimension of distribution                     */
  long skip;                /* skip-parameter of the ball method             */
  double ball_radius;       /* (initial) ball radius                         */
  int adaptive_ball;        /* adaptive flag for ball sampler radius         */
  double adaptive_factor;   /* factor of adaptive radius increase/decrease   */
  double r;                 /* r-parameter for the RoU variant               */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ball_gen {
  int    dim;               /* dimension of distribution                     */
  long   skip;              /* skip-parameter of the ball method             */
  const double *center;     /* center of distribution                        */
  double *direction;        /* random direction vector                       */
  double *point_current;    /* current point inside the shape                */
  double *point_random;     /* random point, can be inside shape or not      */
  double *x;                /* working point in the (xy)-coordinate system   */
  double ball_radius;       /* (initial) ball radius                         */
  int adaptive_ball;        /* adaptive flag for ball sampler radius         */
  double adaptive_factor;   /* factor of adaptive radius increase/decrease   */
  double r;                 /* r-parameter for the RoU variant               */
};

/*---------------------------------------------------------------------------*/

