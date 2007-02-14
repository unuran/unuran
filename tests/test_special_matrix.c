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
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

#include <unur_source.h>
#include <time.h>
#include "../src/utils/matrix_source.h"

/*-------------------------------------------------------------------------*/

/* function prototypes */
void _unur_test_set_matrix_1(int dim, double *M);
void _unur_test_set_matrix_2(int dim, double *M);
void _unur_test_matrix(void);

/* dimension of test matrix */
static const int dim = 50;

/* max tolerable absolute error */
static const double error_absolute_max = 1.e-12;

/* number of failed tests */
static int n_failed = 0;

/* file handles */
static FILE *TESTLOG;         /* test log file     */

/*-------------------------------------------------------------------------*/

void _unur_test_set_matrix_1(int dim, double *M) {
#define idx1(a,b) ((a-1)*dim+(b-1))
  int i,j,k;
  double p;
  
  /* original derflinger test-matrix */
  p = M_PI/(dim+1);
  for (i=1; i<=dim; i++) {
  for (k=1; k<=i; k++) {
    M[idx1(i,k)]=0.;
    for (j=1; j<=dim; j++) {
      M[idx1(i,k)] += sin(p*i*j)*j*sin(p*j*k);
    }
    M[idx1(k,i)] = M[idx1(i,k)]; 
  }}

#undef idx1
} /* end of _unur_test_set_matrix_1() */

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
} /* end of _unur_test_set_matrix_2() */

/*-------------------------------------------------------------------------*/

void _unur_test_matrix(void) {
#define idx1(a,b) ((a-1)*dim+(b-1))

  int i, ret;

  double *M;
  double *values;
  double *vectors;

  double error_absolute; 
  double error_relative; 
 
  int error_counter;
  char error_char;

  /* allocate memory */
  M = malloc(dim*dim*sizeof(double));
  values = malloc(dim*sizeof(double));
  vectors = malloc(dim*dim*sizeof(double));

  _unur_test_set_matrix_1(dim, M);
  
  error_counter = 0;
  ret = _unur_matrix_eigensystem(dim, M, values, vectors);
  if (ret == UNUR_SUCCESS) {
    fprintf(TESTLOG, " #   \t   Eigenvalue \t\t   Expected \t\t Abs. Error \t Rel. Error\n"); 
    for (i=0; i<dim; i++) {
      error_absolute = 2.*values[i]/(dim+1) - (i+1) ;  
      error_relative = error_absolute / (i+1.);

      if (error_absolute > error_absolute_max) { 
        /* eigenvalue is not ok */
        error_char = '*';
        error_counter++;
      }
      else
      {
        /* eigenvalue is ok */
        error_char = ' ';
      }

      fprintf(TESTLOG, "%02d : \t%18.15f   \t%18.15f  \t% 8.4e \t% 8.4e  %c\n", 
              i+1,  
              2.*values[i]/(dim+1),  
              (double) (i+1) , 
              error_absolute,
	      error_relative,
	      error_char);
    }
    if (error_counter==0) {printf("[Eigensystem --> ok]");}
    else printf("[Eigensystem --> failed at %d eigenvalues]", error_counter);
  }
  else {
    printf("[Eigensystem --> failed]");
  }

  n_failed = error_counter;

  /* free memory */
  free(M); 
  free(values);
  free(vectors);

#undef idx1
} /* end of _unur_test_matrix() */

/*-------------------------------------------------------------------------*/

int main(void)
{

  /* open log file for testing */
  if ( (TESTLOG = fopen( "t_special_matrix_test.log","w" )) == NULL )
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


  /* write tail into log file and some info to the console */
  fprintf(TESTLOG,"\n======================================================\n\n");
  if (n_failed > 0) {
    printf("[Matrix --> failed]\n");
    fprintf(TESTLOG,"%d tests FAILED\n", n_failed);
  }
  else {
    printf("[Matrix --> ok]\n");
    fprintf(TESTLOG,"All tests PASSED\n");
  }

  /* close files */
  fclose (TESTLOG);

  /* exit */
  exit( (n_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE );

} /* end of main() */

/*-------------------------------------------------------------------------*/
