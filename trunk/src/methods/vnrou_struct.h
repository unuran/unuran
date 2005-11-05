/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vnrou_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method VNROU                              *
 *         (Vector Naive Ratio Of Uniforms)                                  *
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

struct unur_vnrou_par { 
  int    dim;               /* dimension of distribution                     */
  double r;		    /* r-parameter of the vnrou method 	             */
  double *umin, *umax;      /* boundary rectangle u-coordinates              */
  double vmax;              /* boundary rectangle v-coordinate               */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_vnrou_gen { 
  int    dim;               /* dimension of distribution                     */
  double r;		    /* r-parameter of the vnrou method 	             */
  double *umin, *umax;      /* boundary rectangle u-coordinates              */
  double vmax;              /* boundary rectangle v-coordinate               */
  const double *center;     /* center of distribution                        */  
};

/*---------------------------------------------------------------------------*/

