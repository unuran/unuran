/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ninv_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method NINV                               *
 *         (Numerical INVersion of cumulative distribution function)         *
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

struct unur_ninv_par { 
  int     max_iter;          /* maximal number of iterations                 */
  double  rel_x_resolution;  /* maximal relative error in x                  */
  double  s[2];              /* interval boundaries at start (left/right)    */
  int     table_on;          /* if TRUE a table for starting points is used  */
  int     table_size;        /* size of table                                */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ninv_gen { 
  int     max_iter;          /* maximal number of iterations                 */
  double  rel_x_resolution;  /* maximal relative error in x                  */
  double *table;             /* table with possible starting values for NINV */
  int     table_on;          /* if TRUE a table for starting points is used  */
  int     table_size;        /* size of table                                */
  double  Umin, Umax;        /* bounds for iid random variable in respect to
                                the given (truncated) domain of the distr.   */
  double  CDFmin, CDFmax;    /* CDF-bounds of the table                      */
  double  s[2];              /* interval boundaries at start (left/right) ...*/
  double  CDFs[2];           /* ... and their CDF-values                     */
};

/*---------------------------------------------------------------------------*/
























