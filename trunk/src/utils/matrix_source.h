/********************************/
/* matrix_source.h		*/
/* header file for matrix.c	*/
/********************************/

int _unur_matrix_cholesky_decomposition (int dim, const double *S, double *L );
/* The Colesky factor L of a variance-covariance matrix S is computed: S = LL' */

int _unur_matrix_invert_matrix (int dim, double *A, double detmin, double *Ainv, double *det);
/* Calculates the inverse matrix (by means of LU decomposition).             */
/* Calculates the inverse matrix when |det(A)| > detmin */

void _unur_matrix_debug(int dim, double *M, const char *info, const char *genid);
/* Writes the matrix-elements to the log file */
