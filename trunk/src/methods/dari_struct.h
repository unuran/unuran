/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dari_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method DARI                               *
 *         ((Discrete) Alias-Urn)                                            *
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

struct unur_dari_par { 
  int    squeeze;         /* should the squeeze be used  
			     0.. no squeeze,  1..squeeze                     */
  int    size;            /* size of table for speeding up generation        */
  double c_factor;        /* constant for choosing the design points         */
};

/*---------------------------------------------------------------------------*/
/* The generator object     Das kommentiere ich, wenn ich es programmiere    */

struct unur_dari_gen { 
  double  vt, vc, vcr, xsq[2], y[2], ys[2], ac[2], pm, Hat[2];
  long    m, n[2], x[2], sign[2],s[2];
  double *hp; /*pointer to double array of length size*/
  char   *hb; /*pointer to boolean array of length size*/
};

/*---------------------------------------------------------------------------*/


