/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vector_source.h                                                   *
 *                                                                           *
 *   Routines for computations with vectors (arrays).                        *
 *                                                                           *
 *****************************************************************************
 
 *****************************************************************************
 *                                                                           *
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

/* vector structure */
typedef struct unur_vector
{
  int dim;                   /* dimension of vector */
  double *x;                 /* coordinates of vector */
} UNUR_VECTOR;
    
/*---------------------------------------------------------------------------*/

UNUR_VECTOR *unur_vector_new(int dim);
/* allocate memory for new vector structure */

void unur_vector_free(UNUR_VECTOR *v);
/* free allocated memory used by vector structure */

double unur_vector_norm(UNUR_VECTOR *v);
/* calculation of vector norm */

double unur_vector_scalar_product(UNUR_VECTOR *v1, UNUR_VECTOR *v2);
/* calculation of scalar product */


/*---------------------------------------------------------------------------*/
