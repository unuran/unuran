/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: utdr_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method UTDR                               *
 *         (Universal Transformed Density Rejection; 3-point method)         *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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

struct unur_utdr_par { 
  double  c_factor;     /* constant for choosing the design points           */
  double  delta_factor; /* constant for choosing delta to replace the tangent*/
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_utdr_gen { 

  double  il;                   /* left border of the domain                 */
  double  ir;                   /* right border of the domain                */

  double  vollc,volcompl,voll,fm,hm,
    al,ar,col,cor,sal,sar,bl,br,tlx,trx,
    brblvolc,drar,dlal,ooar2,ooal2;/* constants of the hat and for generation*/
};

/*---------------------------------------------------------------------------*/
