/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: matrix.c                                                          *
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

#include <unur_source.h>
#include "matrix_source.h"

/*---------------------------------------------------------------------------*/

static int _unur_matrix_swap_rows (int dim, double *A, int i, int j);
/* Swap rows i and j in square matrix A of dimension dim.                    */

static int _unur_matrix_permutation_swap (int dim, int *p, int i, int j);
/* Swap entries i and j of the integer array p.                              */

static int _unur_matrix_LU_decomp (int dim, double *A, int *P, int *signum);
/* Factorise a general dim x dim matrix A into P A = L U.                    */

static int _unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X);
/* ????? */

static int _unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X);
/* ????? */

static int _unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse);
/* Compute inverse of matrix with given LU decomposition.                    */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_swap_rows (int dim, double *A, int i, int j)
     /*----------------------------------------------------------------------*/
     /* Swap rows i and j in square matrix A of dimension dim.               */
     /* The resulting matrix is again stored in A.                           */ 
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   A   ... dim x dim -matrix					     */
     /*   i,j ... row numbers that should be swapped			     */
     /*                                                                      */
     /* output:								     */
     /*   the given matrix A contains the new matrix with swapped rows       */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);

  if (i != j)
  {
    double *row1 = A + i * dim;
    double *row2 = A + j * dim;
      
    int k;
      
    for (k = 0; k < dim; k++)
    {
       double tmp = row1[k] ;
       row1[k] = row2[k] ;
       row2[k] = tmp ;
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_matrix_swap_rows() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_permutation_swap (int dim, int *p, int i, int j)
     /*----------------------------------------------------------------------*/
     /* Swap entries i and j of the integer array p.                         */
     /*                                                                      */
     /* input:							             */
     /*   dim ... number of elements				             */
     /*   p   ... permutation vector of length dim		             */
     /*           (contains all integers between 0 and dim-1) 	             */
     /*   i,j ... indexes that should be swapped		             */
     /*                                                                      */
     /* output:							             */
     /*   permutation p with swapped elements  			             */
     /*----------------------------------------------------------------------*/
{
  if (i != j)
    {
      int tmp = p[i];
      p[i] = p[j];
      p[j] = tmp;
    }
  
  return UNUR_SUCCESS;
} /* end of _unur_matrix_permutation_swap() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_LU_decomp (int dim, double *A, int *p, int *signum)
     /*----------------------------------------------------------------------*/
     /* Factorise a general dim x dim matrix A into,			     */
     /*									     */
     /*   P A = L U							     */
     /*									     */
     /* where P is a permutation matrix, L is unit lower triangular and U    */
     /* is upper triangular.						     */
     /*									     */
     /* L is stored in the strict lower triangular part of the input	     */
     /* matrix. The diagonal elements of L are unity and are not stored.     */
     /*									     */
     /* U is stored in the diagonal and upper triangular part of the	     */
     /* input matrix.  							     */
     /* 								     */
     /* P is stored in the permutation p. Column j of P is column k of the   */
     /* identity matrix, where k = permutation->data[j]			     */
     /*									     */
     /* signum gives the sign of the permutation, (-1)^n, where n is the     */
     /* number of interchanges in the permutation. 			     */
     /*									     */
     /* See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss    */
     /* Elimination with Partial Pivoting).				     */
     /*									     */
     /* input:								     */
     /*   dim    ... number of columns and rows of			     */
     /*   A      ... dim x dim matrix					     */
     /*   p      ... pointer to array where permutation should be stored     */
     /*   signum ... pointer where sign of permutation should be stored      */
     /*									     */
     /* output:								     */
     /*   A      ... L (strict lower triangular part) and U (upper part)     */
     /*   p      ... permutation vector of length dim			     */
     /*              (contains all integers between 0 and dim-1) 	     */
     /*   signum ... sign of the permutation (1 or -1)			     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+b)
  int i, j, k;

  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(signum,UNUR_ERR_NULL);

  /* initialize permutation vector */
  *signum = 1;
  for(i=0;i<dim;i++) p[i]=i;

  for (j = 0; j < dim - 1; j++){
    /* Find maximum in the j-th column */

    double ajj, max = fabs (A[idx(j,j)]);
    int i_pivot = j;
    for (i = j + 1; i < dim; i++){
      double aij = fabs (A[idx(i,j)]);
      if (aij > max){
	max = aij;
	i_pivot = i;
      }
    }
    if (i_pivot != j){
      _unur_matrix_swap_rows (dim, A, j, i_pivot);
      _unur_matrix_permutation_swap (dim, p, j, i_pivot);
      *signum = -(*signum);
    }

    ajj = A[idx(j,j)];

    if (ajj != 0.0){
      for (i = j + 1; i < dim; i++){
	double aij = A[idx(i,j)] / ajj;
	A[idx(i,j)] = aij;
	for (k = j + 1; k < dim; k++){
	  double aik = A[idx(i,k)];
	  double ajk = A[idx(j,k)];
	  A[idx(i,k)]= aik - aij * ajk;
	}
      }
    }
  }
      
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_LU_decomp() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_backsubstitution_dtrsv(int dim, double *LU, double *X)
     /*----------------------------------------------------------------------*/
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix that contains LU decomposition of matrix */
     /*   X   ... vector 						     */
     /*                                                                      */
     /* output:								     */
     /*   X 								     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+b)

  int ix,jx,i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);

  /* backsubstitution */
  ix = (dim - 1);
   
  X[ix] = X[ix] / LU[idx(ix,ix)];
 
  ix--;
  for (i = dim - 1; i > 0 && i--;) {
    double tmp = X[ix];
    jx = ix + 1;
    for (j = i + 1; j < dim; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx ++;
    }
      
    X[ix] = tmp / LU[idx(i,i)];
    ix --;
  }

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_backsubstitution_dtrsv() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_forwardsubstitution_dtrsv(int dim, double *LU, double *X)
     /*----------------------------------------------------------------------*/
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix that contains LU decomposition of matrix */
     /*   X   ... vector 						     */
     /*                                                                      */
     /* output:								     */
     /*   X 								     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+b)

  int ix,jx,i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(X,UNUR_ERR_NULL);

  /* forward substitution */
  ix = 0;
  ix++;
  for (i = 1; i < dim; i++) {
    double tmp = X[ix];
    jx = 0;
    for (j = 0; j < i; j++) {
      tmp -= LU[idx(i,j)] * X[jx];
      jx += 1;
    }
    X[ix] = tmp;
    ix ++;
  }
  
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_forwardsubstitution_dtrsv() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_LU_invert (int dim, double *LU, int *p, double *inverse)
     /*----------------------------------------------------------------------*/
     /* Compute inverse of matrix with given LU decomposition.               */
     /*                                                                      */
     /* input:								     */
     /*   dim ... number of columns and rows of				     */
     /*   LU  ... dim x dim -matrix, LU decomposition			     */
     /*   p   ... permutation vector of length dim			     */
     /*           (contains all integers between 0 and dim-1)     	     */
     /*   inverse ... pointer to array where inverse should be stored        */
     /*                                                                      */
     /* output								     */
     /*   inverse ... the inverse matrix 				     */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   error code      otherwise                                          */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+b)

  double *vector;
  int i,j;

  /* check arguments */
  CHECK_NULL(LU,UNUR_ERR_NULL);
  CHECK_NULL(p,UNUR_ERR_NULL);
  CHECK_NULL(inverse,UNUR_ERR_NULL);

  /* allocate working array */
  vector = _unur_malloc(dim*sizeof(double));

  for (i = 0; i < dim; i++){
    for(j=0;j<dim;j++) vector[j] = 0.;
    vector[i] = 1.;
    /* Solve for c using forward-substitution, L c = P b */
    _unur_matrix_forwardsubstitution_dtrsv (dim, LU, vector);
    /* Perform back-substitution, U x = c */
    _unur_matrix_backsubstitution_dtrsv (dim, LU, vector);

    for(j=0;j<dim;j++){
      inverse[idx(j,p[i])] = vector[j]; /* set column vector of inverse matrix */
    }
  }

  /* free working arrays */
  free(vector);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_LU_invert() */

/*---------------------------------------------------------------------------*/

int 
_unur_matrix_invert_matrix(int dim, double *A, double detmin, double *Ainv, double *det)
     /*-------------------------------------------------------------------------*/
     /* Calculates the inverse matrix (by means of LU decomposition)         	*/
     /* If |det(A)| <= detmin a message is printed 				*/
     /* the array Ainv is computed whenever |det(A)|/dim >= UNUR_EPSILON        */
     /*	or whenever norm(A)*dim/|det(A)| <= DBL_MAX/2	 			*/
     /*										*/
     /* input:                                                                  */
     /*   dim    ... dimension of the square matrix A                        	*/
     /*   A      ... dim x dim -matrix                                       	*/
     /*   detmin ... threshold value for |det(A)| for illconditioned matrix  	*/
     /*   Ainv   ... pointer to array where inverse matrix should be stored  	*/
     /*   det    ... pointer where det(A) should be stored                   	*/
     /*                                                                         */
     /* output:                                                                 */
     /*   Ainv   ... inverse matrix of A	                                */
     /*   det    ... determinant of A                                           */
     /*									     	*/
     /* return:								     	*/
     /*   UNUR_SUCCESS on success                                            	*/
     /*   UNUR_FAILURE when matrix is ill-conditioned, i.e. when		*/
     /*                |det(A)|/dim >= UNUR_EPSILON  or                         */
     /*		       norm(A)*dim/|det(A)| <= DBL_MAX/2			*/
     /*                (array Ainv remains unchanged in this case)              */
     /*   other error code, otherwise                                           */
     /*-------------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+b)

  int *p, s, i, j;
  double *LU;
  double norm, halfnorm;
  
  /* check arguments */
  CHECK_NULL(A,UNUR_ERR_NULL);
  CHECK_NULL(Ainv,UNUR_ERR_NULL);
  CHECK_NULL(det,UNUR_ERR_NULL);
  if (dim<2) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension <= 1");
    return UNUR_ERR_GENERIC;
  }

  /* allocate working space */
  p = _unur_malloc(dim*sizeof(int));
  LU = _unur_malloc(dim*dim*sizeof(double));

  /* compute LU decomposition */
  memcpy(LU, A, dim*dim*sizeof(double));
  _unur_matrix_LU_decomp(dim, LU, p, &s);

  /* compute determinat */
  *det = s;
  for(i=0;i<dim;i++)
    *det *= LU[idx(i,i)];

  /* check for small determinant */
  if (fabs(*det) <= detmin) {
    _unur_warning("matrix",UNUR_ERR_GENERIC,"det(A) < detmin");
  }

  /* calculate matrix norm */
  norm=0.;
  halfnorm=0.;
  for(i=0;i<dim;i++) {
    for(j=0;j<dim;j++) {
      halfnorm += fabs(A[idx(i,j)]);
    }
    if (halfnorm > norm) norm=halfnorm;
  }

  if ( fabs(*det)/dim < UNUR_EPSILON ) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"|det(A)| < dim*UNUR_EPSILON"); 
    return UNUR_FAILURE; 
  }
 
  /* check for ill-conditioned matrix */
  if ( norm * dim / fabs(*det) > DBL_MAX / 2 ) { 
    _unur_error("matrix",UNUR_ERR_GENERIC,"matrix not computationally stable");  
    return UNUR_FAILURE; 
  } 
  
  /* compute inverse by means of LU factors */
  _unur_matrix_LU_invert(dim, LU, p, Ainv);   
  
  /* free working space */
  free(LU);
  free(p);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_invert_matrix() */

/*---------------------------------------------------------------------------*/

double 
_unur_matrix_qf (int dim, double *x, double *A)
     /*----------------------------------------------------------------------*/
     /* Compute quadratic form x'Ax                                          */
     /*									     */
     /* input:                                                               */
     /*   dim ... number of columns and rows of                              */
     /*   x   ... vector                                                     */
     /*   A   ... dim x dim matrix                                           */
     /*									     */
     /* return:								     */
     /*   returns the result of x'Ax                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+b)

  int i,j;
  double sum,outersum;
  
  /* check arguments */
  CHECK_NULL(x,INFINITY);
  CHECK_NULL(A,INFINITY);
  if (dim<2) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension <= 1");
    return INFINITY;
  }

  outersum=0.;
  for(i=0;i<dim;i++){
    sum=0.;
    for(j=0;j<dim;j++)
      sum+=A[idx(i,j)]*x[j];
    outersum+=sum*x[i];
  }

  return outersum;

#undef idx
} /* end of _unur_matrix_qf() */

/*---------------------------------------------------------------------------*/

int
_unur_matrix_cholesky_decomposition (int dim, const double *S, double *L )
     /*----------------------------------------------------------------------*/
     /* The Colesky factor L of a variance-covariance matrix S is computed:  */
     /*    S = LL'                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... dimension of covariance matrix S                           */
     /*   S   ... variance-covariance matrix                                 */
     /*   L   ... pointer to array where cholesky factor should be stored    */
     /*                                                                      */
     /* return:								     */
     /*   UNUR_SUCCESS on success                                            */
     /*   UNUR_FAILURE when matrix is not positive definite                  */
     /*   other error code, otherwise                                        */
     /*----------------------------------------------------------------------*/
{ 
#define idx(a,b) ((a)*dim+b)

  int i,j,k;
  double sum1,sum2;

  /* check arguments */
  CHECK_NULL(S,UNUR_ERR_NULL);
  if (dim<2) {
    _unur_error("matrix",UNUR_ERR_GENERIC,"dimension <= 1");
    return UNUR_ERR_GENERIC;
  }

  L[idx(0,0)] = sqrt( S[idx(0,0)] );

  for(j=1; j<dim; j++) {
    L[idx(j,0)] = S[idx(j,0)] / L[idx(0,0)];

    sum1 = L[idx(j,0)] * L[idx(j,0)];
    for(k=1; k<j; k++) {
      sum2 = 0.;
      for(i=0; i<k; i++)
	sum2 += L[idx(j,i)] * L[idx(k,i)];
      
      L[idx(j,k)] = (S[idx(j,k)] - sum2) / L[idx(k,k)];
      sum1 += L[idx(j,k)] * L[idx(j,k)];
    }

    if (S[idx(j,j)] <= sum1) {
      /* covariance matrix not positive definite */
      return UNUR_FAILURE;
    }

    L[idx(j,j)] = sqrt( S[idx(j,j)] - sum1 );
  }

  /* although not necessary upper triangular of L - matrix is set to 0 */
  for(j=0; j<dim; j++)
    for(k=j+1; k<dim; k++)
      L[idx(j,k)]=0.;

  /* o.k. */
  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_matrix_cholesky_decomposition() */


/*--------------------------------------------------------------------------*/

void 
_unur_matrix_debug ( int dim, const double *M, const char *info, const char *genid )
     /*----------------------------------------------------------------------*/
     /* The elemets of the rectangular dim x dim matrix M                    */
     /* are written row-wise into the logfile.                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim   ... dimension                                                */
     /*   M     ... rectangular matrix with dim rows and columns             */
     /*   info  ... additional info-string to be printed                     */
     /*   genid ... id string                                                */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+b)

  FILE *log;
  int i,j;
 
  log = unur_get_stream();

  fprintf(log,"%s: %s\n", genid, info); 

  if (M==NULL) {
    fprintf(log,"%s: NULL pointer\n", genid);
  }

  else {
    for (i=0; i<dim; i++) {
      fprintf(log, "%s: ", genid); 
      for (j=0; j<dim; j++) {
        M[idx(i,j)]<0 ? 
        fprintf(log, " %e",  M[idx(i,j)]): 
        fprintf(log, "  %e", M[idx(i,j)]); 
      }
      fprintf(log, "\n");
    }
  }

  fprintf(log,"%s:\n", genid); 

#undef idx
} /* end of _unur_matrix_debug() */  

/*---------------------------------------------------------------------------*/

