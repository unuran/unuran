/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      printsample.c                                                *
 *                                                                           *
 *   print a sample of random numbers                                        *
 *                                                                           *
 *****************************************************************************
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

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/
static char test_name[] = "Sample";
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

void
unur_test_printsample( struct unur_gen *gen, int n_rows, int n_cols )
     /*----------------------------------------------------------------------*/
     /* print a sample of generator output in small (n_rows x n_cols) table  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   n_rows ... number of rows                                          */
     /*   n_cols ... number of columns (for univariate case only)            */
     /*----------------------------------------------------------------------*/
{
  int i,j,k;

  /* check arguments */
  _unur_check_NULL(test_name,gen,/*void*/);

  printf("\nSAMPLE: ");              

  switch (gen->method & UNUR_MASK_TYPE) {
  
  case UNUR_METH_DISCR:
    for( j=0; j<n_rows; j++ ) {
      for(i=0; i<n_cols; i++)
	printf("%04d ",_unur_sample_discr(gen));
      printf("\n        "); 
    }
    break;

  case UNUR_METH_CONT:
    for( j=0; j<n_rows; j++ ) {
      for(i=0; i<n_cols; i++)
	printf("%8.5f ",_unur_sample_cont(gen));
      printf("\n        "); 
    }
    break;

  case UNUR_METH_VEC:
    { /* we need an array for the vector */
      double *vec;
      int dim;
      dim = unur_get_dimension(gen);
      vec = _unur_malloc( dim * sizeof(double) );
	
      for( j=0; j<n_rows; j++ ) {
	_unur_sample_vec(gen,vec);
	printf("( %8.5f",vec[0]);
	for (k=1; k<dim; k++)
	  printf(", %8.5f",vec[k]);
	printf(" )\n        ");
      }
      free(vec);
    }
    break;
  default: /* unknown ! */
    _unur_error(test_name,UNUR_ERR_GENERIC,"method unknown!");
    return;
  }

  printf("\n");

} /* end of unur_test_printsample() */

/*---------------------------------------------------------------------------*/






