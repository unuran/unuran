/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: varou_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method VAROU                              *
 *         (Vector Adaptive Ratio Of Uniforms)                               *
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
/* Cone structure */
struct unur_varou_cone {
  long   *index;   /* indices of the spanning cone vectors (unit verteces)   */
  double *length;    /* lengths of the spanning vectors                      */
  double *spoint;  /* touching point of surface tangential plane             */
  double *normal;  /* normal vector to tangential plane through spoint       */
  double unit_volume; /* volume of cone spanned by unit vectors              */
  double volume;   /* volume of cone : spanning vectors + tangential surface */
};

/*---------------------------------------------------------------------------*/
/* Parameter object                                                          */

struct unur_varou_par { 
  int    dim;               /* dimension of distribution                     */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_varou_gen { 
  int    dim;               /* dimension of distribution */
  double *umin, *umax;      /* boundary rectangle u-coordinates */
  double vmax;              /* boundary rectangle v-coordinate  */
  const double *center;     /* center of distribution */
  int aux_dim;              /* parameter used my auxiliary functions */
  long n_vertex;            /* number of verteces on upper half-unit-sphere */
  double **vertex_list;     /* list of verteces on upper half-unit-sphere */
  long n_cone;              /* number of constructed cones */
  struct unur_varou_cone **cone_list; /* list of constructed cones */
};

/*---------------------------------------------------------------------------*/
