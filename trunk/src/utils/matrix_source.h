/********************************/
/* matrix_source.h		*/
/* header file for matrix.c	*/
/********************************/

double *_unur_matrix_cholesky_decomposition( int dim, const double *S );
/* The Cholesky factor L of a variance-covariance matrix S is computed */

int _unur_matrix_invert_matrix(int dim, double *A, double detmin, double *Ainv, double *det);
/* Calculates the inverse matrix when |det(A)| > detmin */

double _unur_matrix_quadratic_form(int dim, double *x, double *A);
/* Calculates the result of x'Ax */

void _unur_matrix_debug(int dim, double *M, char *info, char *genid); 
/* Writes the matrix-elements to the log file */
