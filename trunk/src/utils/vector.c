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

/* allocate memory for new vector */
double *
_unur_vector_new(int dim)
{
  int i;
  double *v;

  v = _unur_xmalloc(dim*sizeof(double));

  /* setting all coordinates to 0 */
  for (i=0; i<dim; i++) v[i]=0.;

  return v;
} /* end of _unur_vector_new() */

/*--------------------------------------------------------------------------*/

/* free allocated memory used by vector structure */
void 
_unur_vector_free(double *v)
{
  if (v) {
     free(v);
     v=NULL;
  }
} /* end of _unur_vector_free() */

/*--------------------------------------------------------------------------*/

/* calculation of vector norm */
double 
_unur_vector_norm(int dim, double *v)
{
  int i;
  double norm=0.;
  double vmax=0;
  double p;

  /* checking if v is NULL */
  /* TODO: warning ? */
  if (v==NULL) return 0; 

  /* determining the largest element (absolute values) */
  for (i=0; i<dim; i++) {
    if (vmax < fabs(v[i])) vmax = fabs(v[i]); 
  }
  
  if (vmax<=UNUR_EPSILON) return 0;
  
  /* it's nummerically more stable to calculate the norm this way */
  for (i=0; i<dim; i++) {
    p=v[i]/vmax;
    norm += p*p;
  }
  norm = vmax * sqrt(norm);

  return norm;
} /* end of _unur_vector_norm() */

/*--------------------------------------------------------------------------*/

/* calculation of scalar product */
double 
_unur_vector_scalar_product(int dim, double *v1, double *v2)
{
  int i;
  double scalar_product=0.;
  
  /* checking if v1 or v2 are NULL */
  /* TODO: warning ? */
  
  for (i=0; i<dim; i++) {
    scalar_product += v1[i]*v2[i];
  }

  return scalar_product;
} /* end of _unur_vector_scalar_product() */

/*--------------------------------------------------------------------------*/

