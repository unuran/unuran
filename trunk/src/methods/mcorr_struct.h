/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mcorr_struct.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method MCORR                              *
 *         (Matrix -- CORRelation matrix)                                    *
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

struct unur_mcorr_par { 
  int    dim;          /* dimension (number of rows and columns) of matrix */
  double *eigenvalues; /* optional eigenvalues of the correlation matrix   */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_mcorr_gen { 
  int    dim;          /* dimension (number of rows and columns) of matrix */
  double *H;           /* working array                                    */
  double *eigenvalues; /* optional eigenvalues of the correlation matrix   */
};

/*---------------------------------------------------------------------------*/

