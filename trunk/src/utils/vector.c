/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vector.c                                                          *
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

#include <unur_source.h>
#include "vector_source.h"

/*---------------------------------------------------------------------------*/

/* allocate memory for new vector structure */
UNUR_VECTOR *
unur_vector_new(int dim)
{
  int i;
 
  UNUR_VECTOR *v;
  v = _unur_xmalloc(sizeof(UNUR_VECTOR));

  v->dim = dim;
  v->x = _unur_xmalloc(dim*sizeof(double));

  /* setting all coordinates to 0 */
  for (i=0; i<dim; i++) v->x[i]=0.;

  return v;
}

/*--------------------------------------------------------------------------*/

/* free allocated memory used by vector structure */
void 
unur_vector_free(UNUR_VECTOR *v)
{
  if (v) {
     if (v->x) free(v->x);
     free(v);
     v=NULL;
  }
}

/*--------------------------------------------------------------------------*/

/* calculation of vector norm */
double 
unur_vector_norm(UNUR_VECTOR *v)
{
  int i;
  double norm=0.;

  /* TODO: checking if v or v->x are NULL */

  for (i=0; i<v->dim; i++) {
    norm += v->x[i]*v->x[i];
  }
  norm = sqrt(norm);

  return norm;
}

/*--------------------------------------------------------------------------*/

/* calculation of scalar product */
double 
unur_vector_scalar_product(UNUR_VECTOR *v1, UNUR_VECTOR *v2)
{
  int i;
  double scalar_product=0.;
  
  /* TODO: checking if vi or vi->x are NULL or have unequal dimensions */
  
  for (i=0; i<v1->dim; i++) {
    scalar_product += v1->x[i]*v2->x[i];
  }

  return scalar_product;
}

/*--------------------------------------------------------------------------*/

