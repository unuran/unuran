/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: matrix_source.h                                                   *
 *                                                                           *
 *   Routines for computations with square matrices.                         *
 *   Matrices are stored as double array (i.e.: double *matrix).             *
 *   The rows of the matrix have to be stored consecutively in this array,   *
 *   i.e. the entry with index [i,j] can be entered via (i*dim+j), where     *
 *   dim is the dimension of the dim x dim - matrix.                         *
 *                                                                           *
 *   Routines are mainly adapted from the GSL (GNU Scientifiy Library)       *
 *                                                                           *
 *****************************************************************************
    $Id$
 *****************************************************************************
 *                                                                           *
 *   adapted by Wolfgang Hoermann and Josef Leydold                          *
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
				                                                                                    
/*--------------------------------------------------------------------------*/
int _unur_matrix_cholesky_decomposition (int dim, const double *S, double *L );
/* The Colesky factor L of a variance-covariance matrix S is computed: S = LL' */

int _unur_matrix_invert_matrix (int dim, double *A, double detmin, double *Ainv, double *det);
/* Calculates the inverse matrix (by means of LU decomposition).             */
/* Calculates the inverse matrix when |det(A)| > detmin */

void _unur_matrix_debug(int dim, double *M, const char *info, const char *genid);
/* Writes the matrix-elements to the log file */
