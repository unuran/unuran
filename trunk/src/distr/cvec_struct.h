/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cvec_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines structures and function prototypes for                    *
 *         handling multivariate distributions.                              *
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
#ifndef UNUR_CVEC_STRUCT_H_SEEN
#define UNUR_CVEC_STRUCT_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* In the following we describe a hyper-rectangle :
   [xlow[0], xhigh[0]] x ... x [xlow[dim-1], xhigh[dim-1]]
   where xlow[] and xhigh[] contain the (finite or unbounded) coordinates of 
   the two vertices (with lowest/highest coordinates) that span the rectangle.
*/


/*---------------------------------------------------------------------------*/

struct unur_rectangle {
  double *xlow; 
  double *xhigh;
};

/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
#endif   /* UNUR_CVEC_STRUCT_H_SEEN */
/*---------------------------------------------------------------------------*/
