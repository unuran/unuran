/*
 * matrix_source.h
 * header file for matrix.c
 *
 */


double *_unur_matrix_cholesky_decomposition( double *S, int dim );
/*
   The Cholesky factor L of a variance-covariance matrix S is computed   
   (the necessary array is allocated)                                   
                                                                        
   parameters:                                                          
     S   ... variance-covariance matrix                                 
     dim ... dimension (S is a dim x dim matrixes)                      
                                                                        
   return:                                                              
     pointer to Cholesky factor                                         
                                                                        
   error:                                                               
     return NULL                                                        
*/

/********************************************************/
int _unur_matrix_invert_matrix(int dim, double *A, double *res, double *det);
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

/*****************************************************/
double _unur_matrix_quadratic_form(int dim, double *x, double *A);
/*
 * input:
 *   dim ... number of columns and rows of
 *   x ... vector
 *   A ... dim x dim -matrix
 * output:
 *   returns the result of x'Ax
 */


void _unur_matrix_debug(double *M, int dim, const char *info); 

