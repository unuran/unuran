/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hitrou_struct.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method HITROU                             *
 *         (Multivariate HIT and run Ratio Of Uniforms)                      *
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

struct unur_hitrou_par {
  int    dim;               /* dimension of distribution                     */
  double r;                 /* r-parameter of the hitrou method              */
  long skip;                /* skip-parameter of the hitrou method           */
  int bounding_rectangle;   /* if we calculate and use the bounding rect     */ 
  double *umin, *umax;      /* boundary rectangle u-coordinates              */
  double vmax;              /* boundary rectangle v-coordinate               */
  int adaptive_points;      /* adaptive reusability of outside points (0/1)  */
  int adaptive_strip;       /* adaptive reusability of strip position : vmax */
  int unidirectional_flag;  /* should we sample uni- or bi-directionally     */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_hitrou_gen {
  int    dim;               /* dimension of distribution                     */
  double r;                 /* r-parameter of the hitrou method              */
  long skip;                /* skip-parameter of the hitrou method           */
  int bounding_rectangle;   /* if we calculate and use the bounding rect     */ 
  double *umin, *umax;      /* boundary rectangle u-coordinates              */
  double vmax;              /* boundary rectangle/strip v-coordinate         */
  const double *center;     /* center of distribution                        */
  double *direction;        /* random direction vector                       */
  double *point_current;    /* current point inside the shape                */
  double *point_random;     /* random point, can be inside shape or not      */
  double *x;                /* working point in the (xy)-coordinate system   */
  int coordinate;           /* current coordinate used for stepping          */
  int adaptive_points;      /* adaptive reusability of outside points (0/1)  */
  int adaptive_strip;       /* adaptive reusability of strip position : vmax */
  int unidirectional_flag;  /* should we sample uni- or bi-directionally     */
};

/*---------------------------------------------------------------------------*/

