/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: test_special_matrix.c                                             *
 *                                                                           *
 *   Tests for matrix and eigenvalues/vectors calculations                   *
 *                                                                           *
 *****************************************************************************/

#include <unur_source.h>
#include <time.h>
#include "../src/utils/matrix_source.h"
void _unur_test_set_matrix_1(int dim, double *M);
void _unur_test_set_matrix_2(int dim, double *M);

/* number of failed tests */
static int n_failed = 0;

/* file handles */
static FILE *UNURANLOG;       /* unuran log file   */
static FILE *TESTLOG;         /* test log file     */

/*-------------------------------------------------------------------------*/

void _unur_test_set_matrix_1(int dim, double *M) {
#define idx1(a,b) ((a-1)*dim+(b-1))
  int i,j,k;
  double p;
  
  /* original derflinger test-matrix */
  p=3.141592653589793/(dim+1);
  for (i=1; i<=dim; i++) {
  for (k=1; k<=i; k++) {
    M[idx1(i,k)]=0.;
    for (j=1; j<=dim; j++) {
      M[idx1(i,k)] += sin(p*i*j)*j*sin(p*j*k);
    }
    M[idx1(k,i)] = M[idx1(i,k)]; 
  }}

#undef idx1
}

/*-------------------------------------------------------------------------*/

void _unur_test_set_matrix_2(int dim, double *M) {
#define idx1(a,b) ((a-1)*dim+(b-1))
  int i,j,k;

  /* another test-matrix */
  for (i=1; i<=dim; i++) {
  for (k=1; k<=i; k++) {
    M[idx1(i,k)]=0.;
    for (j=1; j<=dim; j++) {
      M[idx1(i,k)] = 1./(i+k);
    }
    M[idx1(k,i)] = M[idx1(i,k)]; 
  }}

#undef idx1
}

/*-------------------------------------------------------------------------*/

void _unur_test_matrix() {
#define idx1(a,b) ((a-1)*dim+(b-1))

  int i, dim, ret;

  double *M;
  double *values;
  double *vectors;

  int error_counter;
  
  dim=50; 
  M=malloc(dim*dim*sizeof(double));
  values=malloc(dim*sizeof(double));
  vectors=malloc(dim*dim*sizeof(double));

  _unur_test_set_matrix_1(dim, M);
  
  error_counter = 0;
  ret = _unur_matrix_eigensystem(dim, M, values, vectors);
  if (ret == UNUR_SUCCESS) {
    for (i=0; i<dim; i++) {
      fprintf(TESTLOG, "Eigenvalue #%02d : %e  (scaled %e)\n", i+1, values[i], 
      2.*values[i]/(dim+1));
      if ( _unur_FP_approx(2.*values[i]/(dim+1), i+1) ) {
        /* eigenvalue is ok */
      } 
      else {
        /* eigenvalue is not ok */
	error_counter++;
      }
    }
    if (error_counter==0) {printf("[Eigensystem --> ok]");}
    else printf("[Eigensystem --> failed at %d eigenvalues]", error_counter);
  }
  else {
    printf("[Eigensystem --> failed]");
  }
 
  /* TODO ? additional tests using _unur_test_set_matrix_2 */
  
#undef idx1
}


/*-------------------------------------------------------------------------*/


int main()
{

  /* open log file for unuran and set output stream for unuran messages */
  if ( (UNURANLOG = fopen( "test_special_unuran.log","w" )) == NULL )
    exit (EXIT_FAILURE);
  unur_set_stream( UNURANLOG );
                                                                                                           
  /* open log file for testing */
  if ( (TESTLOG = fopen( "test_special_matrix.log","w" )) == NULL )
    exit (EXIT_FAILURE);
                                                                                                           
  /* write header into log file */
  {
    time_t started;
    fprintf(TESTLOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");
    if (time( &started ) != -1)
      fprintf(TESTLOG,"%s",ctime(&started));
    fprintf(TESTLOG,"\n======================================================\n\n");
  }
                                                                                                           
  /* run tests */
  _unur_test_matrix();

  /* close files */
  fclose (UNURANLOG);
  fclose (TESTLOG);
                                                                                                           
  /* end */
  if (n_failed > 0) {
    printf("[Matrix --> failed]\n");
    exit (EXIT_FAILURE);
  }
  else {
    printf("[Matrix --> ok]\n");
    exit (EXIT_SUCCESS);
  }
  return 0;
}

