/*
 * matrix.c
 *
 * implements the inversion of matrices using LU decomposition
 * and the product x'A x 
 *
 * Code by WH 21.1.04 using and changing code of the GSL library
 * Modifications and unuran adaption by RK
 *
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <unur_source.h>
#include "unur_errno.h"

int _unur_matrix_swap_rows(int dim, double *A, int i, int j)
/*
 * input:
 *   dim ... number of columns and rows of
 *   A ... dim x dim -matrix
 *   i,j ... row numbers that should be swapped
 * output:
 *   matrix A with swapped rows
 *
 *******************************/
{

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
}

/*********************************************/

int _unur_matrix_permutation_swap (int dim, int * p, int i, int j)
/*
 * input:
 *   dim ... number of elements
 *   p ... permutation vector (length dim, contains all integers between 0 and dim-1) 
 *   i,j ... indexes that should be swapped
 * output:
 *   permutation p with swapped elements  
 */

{
  if (i != j)
    {
      int tmp = p[i];
      p[i] = p[j];
      p[j] = tmp;
    }
  
  return UNUR_SUCCESS;
}


/*********************************************/


/* Factorise a general dim x dim matrix A into,
 *
 *   P A = L U
 *
 * where P is a permutation matrix, L is unit lower triangular and U
 * is upper triangular.
 *
 * L is stored in the strict lower triangular part of the input
 * matrix. The diagonal elements of L are unity and are not stored.
 *
 * U is stored in the diagonal and upper triangular part of the
 * input matrix.  
 * 
 * P is stored in the permutation p. Column j of P is column k of the
 * identity matrix, where k = permutation->data[j]
 *
 * signum gives the sign of the permutation, (-1)^n, where n is the
 * number of interchanges in the permutation. 
 *
 * See Golub & Van Loan, Matrix Computations, Algorithm 3.4.1 (Gauss
 * Elimination with Partial Pivoting).
 */

int _unur_matrix_LU_decomp (int dim, double *A, int *p, int *signum)
/*
 * input:
 *   dim ... number of columns and rows of
 *   A ... dim x dim -matrix
 * output:
 *   A
 *   p ... permutation vector (length dim, contains all integers between 0 and dim-1) 
 *   signum ... sign of the permuttion (1 or -1)
 *
 */

{
  int i, j, k;

  *signum = 1;
  for(i=0;i<dim;i++) p[i]=i;

  for (j = 0; j < dim - 1; j++){
    /* Find maximum in the j-th column */

    double ajj, max = fabs (A[dim*j+j]);
    int i_pivot = j;
    for (i = j + 1; i < dim; i++){
      double aij = fabs (A[dim*i+j]);
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

    ajj = A[dim*j+j];

    if (ajj != 0.0){
      for (i = j + 1; i < dim; i++){
	double aij = A[dim*i+j] / ajj;
	A[dim*i+j] = aij;
	for (k = j + 1; k < dim; k++){
	  double aik = A[dim*i+k];
	  double ajk = A[dim*j+k];
	  A[dim*i+k]= aik - aij * ajk;
	}
      }
    }
  }
      
  return UNUR_SUCCESS;
}

/**********************************************/

int _unur_matrix_backsubstitution_dtrsv(int dim, double *A, double *X)
/*                                (NonUnit)
 * input:
 *   dim ... number of columns and rows of
 *   A ... dim x dim -matrix
 *   X ... vector 
 * output:
 *   X 
 *
 */
{
  int ix,jx,i,j;

  /* backsubstitution */
  ix = (dim - 1);
   
  X[ix] = X[ix] / A[dim * (dim - 1) + (dim - 1)];/*as nonunit*/
 
  ix --;
  for (i = dim - 1; i > 0 && i--;) {
    double tmp = X[ix];
    jx = ix + 1;
    for (j = i + 1; j < dim; j++) {
      double Aij = A[dim * i + j];
      tmp -= Aij * X[jx];
      jx ++;
    }
      
    X[ix] = tmp / A[dim * i + i];/*as nonunit*/
    ix --;
  }

  return UNUR_SUCCESS;
}

/*********************************************/

int _unur_matrix_forwardsubstitution_dtrsv(int dim, double *A, double *X)
/*                                (Unit)
 * input:
 *   dim ... number of columns and rows of
 *   A ... dim x dim matrix
 *   X ... vector 
 * output:
 *   X 
 *
 */
{ 
  int ix,jx,i,j;

  /* forward substitution */
  ix = 0;
  ix ++;
  for (i = 1; i < dim; i++) {
    double tmp = X[ix];
    jx = 0;
    for (j = 0; j < i; j++) {
      double Aij = A[dim * i + j];
      tmp -= Aij * X[jx];
      jx += 1;
    }
    X[ix] = tmp;
    ix ++;
  }
  
  return UNUR_SUCCESS;
}

/*********************************************/

int _unur_matrix_LU_invert (int dim, double * LU, int * p, double * inverse)
/*
 * input:
 *   dim ... number of columns and rows of
 *   LU  ... dim x dim -matrix, LU decomposition
 *   p   ... permutation vector (length dim, contains all integers between 0 and dim-1) 
 * output
 *   inverse ... the inverse matrix 
 *
 */
{ 
  double *vector;
  int i,j;

  vector = malloc(dim*sizeof(double));


  for (i = 0; i < dim; i++){
    for(j=0;j<dim;j++) vector[j] = 0.;
    vector[i] = 1.;
    /* Solve for c using forward-substitution, L c = P b */
    _unur_matrix_forwardsubstitution_dtrsv (dim, LU, vector);
    /* Perform back-substitution, U x = c */
    _unur_matrix_backsubstitution_dtrsv (dim, LU, vector);


    for(j=0;j<dim;j++){
      inverse[dim*j+p[i]] = vector[j];/*set column vector of inverse matrix*/
    }
  }
  free(vector);

  return UNUR_SUCCESS;
}

/********************************************************/

int _unur_matrix_invert_matrix(int dim, double *A, double *res, double *det)
/*
 * input:
 *   dim ... number of columns and rows of
 *   A ... dim x dim -matrix
 * output:
 *   res... inverse matrix of A
 *   det... determinant of A
 *
 * calculates the inverse matrix. If succesfull UNUR_SUCCESS is returned
 * If |det(A)| < 1.e-10 a message is printed and UNUR_FAIL is returned.
 * the array res remains unchanged in this case
 */
{ 
  FILE *log;
  
  int *p, s,i;
  double *lu;

  log = unur_get_stream();

  p = malloc(dim*sizeof(int));
  lu = malloc(dim*dim*sizeof(double));
  for(i=0;i<dim*dim;i++)
    lu[i] = A[i];

  _unur_matrix_LU_decomp(dim, lu, p, &s);
  *det = s;
  for(i=0;i<dim;i++)
    *det *=lu[dim*i+i];
  if(fabs(*det)<1.e-10) {
#ifdef UNUR_ENABLE_LOGGING
    fprintf(log,"Determinant of matrix is close to 0\n");  
    fprintf(log,"Inverse could not be calculated\n");  
    /* _unur_error(); */
#endif
    return UNUR_FAILURE;
  }

  _unur_matrix_LU_invert(dim, lu, p, res);   
  
  free(lu);
  free(p);

  return UNUR_SUCCESS;
}

/*****************************************************/

double _unur_matrix_quadratic_form(int dim, double *x, double *A)
/*
 * input:
 *   dim ... number of columns and rows of
 *   x ... vector
 *   A ... dim x dim -matrix
 * output:
 *   returns the result of x'Ax
 */
{
  int i,j;
  double sum,outersum;
  
  outersum=0.;
  for(i=0;i<dim;i++){
    sum=0.;
    for(j=0;j<dim;j++)
      sum+=A[dim*i+j]*x[j];
    outersum+=sum*x[i];
  }

  return outersum;
}


/*****************************************************/


double *
_unur_matrix_cholesky_decomposition( double *S, int dim )
/*----------------------------------------------------------------------*/
/* the Colesky factor L of a variance-covariance matrix S is computed   */
/* (the necessary array is allocated)                                   */
/*                                                                      */
/* parameters:                                                          */
/*   S   ... variance-covariance matrix                                 */
/*   dim ... dimension (S is a dim x dim matrix)                        */
/*                                                                      */
/* return:                                                              */
/*   pointer to cholesky factor                                         */
/*                                                                      */
/* error:                                                               */
/*   return NULL                                                        */
/*----------------------------------------------------------------------*/
{ 
#define idx(a,b) (a*dim+b)

  double *L;
  int i,j,k;
  double sum1,sum2;

  /* allocate memory for cholesky factor */
  L = _unur_malloc( dim * dim * sizeof(double) );

  /* run cholesky decomposition */
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
      free(L); return NULL;
    }

    L[idx(j,j)] = sqrt( S[idx(j,j)] - sum1 );
  }

  /* although not necessary upper triangular of L - matrix is set to 0 */
  for(j=0; j<dim; j++)
    for(k=j+1; k<dim; k++)
      L[idx(j,k)]=0.;

  /* return (pointer to) cholesky factor */
  return L;

#undef idx
} /* end of _unur_matrix_cholesky_decomposition() */

void _unur_matrix_debug( double *M, int dim, const char *info )
/*----------------------------------------------------------------------*/
/* The elemets of the rectangular dim x dim matrix M                    */
/* are written row-wise into the logfile.                               */
/*                                                                      */
/* parameters:                                                          */
/*   M    ... rectangular matrix with dim rows and columns              */
/*   dim  ... dimension                                                 */
/*   info ... additional info-string to be printed                      */
/*                                                                      */
/*----------------------------------------------------------------------*/
#define idx(a,b) (a*dim+b)
{
  FILE *log;
  int i,j;
  
  log = unur_get_stream();

  fprintf(log,"matrix debug : %s\n", info); 

  if (M==NULL) {
    fprintf(log,"Matrix -> NULL pointer\n");
  }

  else {
    for (i=0; i<dim; i++) {
      fprintf(log, "     (%e", M[idx(i,0)]); 
      for (j=1; j<dim; j++) {
        fprintf(log, ",%e", M[idx(i,j)]); 
      }
      fprintf(log, ")\n");
    }
  }

  fprintf(log,"\n"); 

#undef idx
} /* end of _unur_matrix_debug */  
