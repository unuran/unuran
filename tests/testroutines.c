/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      testroutines.c                                               *
 *                                                                           *
 *   Common test routines                                                    *
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

/*---------------------------------------------------------------------------*/

#include "testunuran.h"

/*---------------------------------------------------------------------------*/
/* function prototypes for internal UNURAN functions                         */

int _unur_isnan (const double x);

/*---------------------------------------------------------------------------*/
/* check whether unur_urng_reset() works                                     */

static int 
cannot_compare_sequence ( FILE *LOG )
{
  fprintf (LOG,"\nURNG cannot be reset. Cannot compare sequences. (Skip)\n");
  printf ("URNG cannot be reset. Cannot compare sequences. (Skip)\n");

  return UNUR_SUCCESS; /* indicate as "not failed" for practical reasons */
}

#define reset_sequence(urng) \
   if (unur_urng_reset((urng)) != UNUR_SUCCESS) {return UNUR_SUCCESS;}
#define reset_sequence_warning(urng) \
   if (unur_urng_reset((urng)) != UNUR_SUCCESS) {return cannot_compare_sequence(LOG);}
#define reset_sequence_par(urng,par) \
   if (unur_urng_reset((urng)) != UNUR_SUCCESS) {if (par) free(par); return UNUR_SUCCESS;}
#define reset_sequence_par_warning(urng,par) \
   if (unur_urng_reset((urng)) != UNUR_SUCCESS) {if (par) free(par); return cannot_compare_sequence(LOG);}

/*---------------------------------------------------------------------------*/
/* check for invalid NULL pointer, that should not happen in this program */

void abort_if_NULL( FILE *LOG, int line, const void *ptr )
{
  if (ptr) return; /* o.k. */
  
  /* 
     There must not be a NULL pointer.
     Since we do not expect a NULL pointer something serious has
     happend. So it is better to abort the tests.
  */
  fprintf(LOG,"line %4d: Unexpected NULL pointer. Panik --> Abort tests!!\n\n",line);
  printf(" Panik --> Tests aborted\n"); fflush(stdout);
  
  /* test finished */
  printf("\n");  fflush(stdout);
  
  /* close log files and exit */
  fclose(LOG);
  exit(EXIT_FAILURE);
  
} /* abort_if_NULL() */

/*---------------------------------------------------------------------------*/
/* compare error code */

int check_errorcode( FILE *LOG, int line, int cherrno )
{
  fprintf(LOG,"line %4d: Error code ...\t\t",line);

  if (unur_errno != cherrno) {
    fprintf(LOG," Failed");
    fprintf(LOG," (observed = %#x, expected = %#x)\n",unur_errno,cherrno);
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_errorcode() */

/*---------------------------------------------------------------------------*/
/* check for expected NULL pointer */

int do_check_expected_NULL( FILE *LOG, int line, int is_NULL )
{
  fprintf(LOG,"line %4d: NULL pointer expected ...\t",line);

  if (is_NULL) {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }

  else {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

} /* end of check_expected_NULL() */

/*---------------------------------------------------------------------------*/
/* check for "set failed" */

int check_expected_setfailed( FILE *LOG, int line, int rcode )
{
  fprintf(LOG,"line %4d: `failed' expected ...\t",line);
  if (rcode==UNUR_SUCCESS) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_setfailed() */

/*---------------------------------------------------------------------------*/
/* check for O (zero) */

int check_expected_zero( FILE *LOG, int line, int k )
{
  fprintf(LOG,"line %4d: 0 (zero) expected ...\t",line);

  if (k != 0) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_zero() */

/*---------------------------------------------------------------------------*/
/* check for INFINITY */

int check_expected_INFINITY( FILE *LOG, int line, double x )
{
  fprintf(LOG,"line %4d: INFINITY expected ...\t",line);

  if (x < UNUR_INFINITY) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_INFINITY() */

int check_expected_negINFINITY( FILE *LOG, int line, double x )
{
  fprintf(LOG,"line %4d: -INFINITY expected ...\t",line);

  if (x > -UNUR_INFINITY) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_negINFINITY() */

int check_expected_INTMAX( FILE *LOG, int line, int k )
{
  fprintf(LOG,"line %4d: INT_MAX expected ...\t",line);

  if (k < INT_MAX) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_INTMAX() */

/*---------------------------------------------------------------------------*/
/* check for reinit */

int check_expected_reinit( FILE *LOG, int line, int rcode )
{
  fprintf(LOG,"line %4d: reinit ...\t\t\t",line);

  if (rcode!=UNUR_SUCCESS) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_reinit() */

/*---------------------------------------------------------------------------*/
/* check for non existing reinit */

int check_expected_no_reinit( FILE *LOG, int line, int rcode )
{
  fprintf(LOG,"line %4d: no reinit ...\t\t",line);

  if (rcode==UNUR_SUCCESS) {
    fprintf(LOG," Failed\n");
    fflush(LOG);
    return UNUR_FAILURE;
  }

  else {
    fprintf(LOG," ok\n");
    fflush(LOG);
    return UNUR_SUCCESS;
  }
} /* end of check_expected_no_reinit() */

/*---------------------------------------------------------------------------*/
/* compare double sequences generated by generator */
/* only when when we can reset the uniform RNG     */

static double *double_sequence_A = NULL;

int compare_double_sequence_par_start( FILE *LOG, int line, UNUR_PAR *par, int sample_size )
{
  int i;
  UNUR_GEN *gen;

  /* allocate memory for storing sequence */
  if (double_sequence_A == NULL) {
    double_sequence_A = malloc( sample_size * sizeof(double) );
    abort_if_NULL(LOG,line, double_sequence_A);
  }

  /* generate sequence */
  reset_sequence_par(unur_get_default_urng(),par);
  gen = unur_init( par ); abort_if_NULL(LOG,line,gen);
  
  for (i=0; i<sample_size; i++)
    double_sequence_A[i] = unur_sample_cont(gen);
  
  unur_free(gen); 

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_double_sequence_par_start() */

/*...........................................................................*/

int compare_double_sequence_urng_start( FILE *LOG, int line, int sample_size )
{
  int i;
  UNUR_URNG *urng = unur_get_default_urng();
  
  /* allocate memory for storing sequence */
  if (double_sequence_A == NULL) {
    double_sequence_A = malloc( sample_size * sizeof(double) );
    abort_if_NULL(LOG, line, double_sequence_A);
  }

  /* generate sequence */
  reset_sequence(urng);
  for (i=0; i<sample_size; i++)
    double_sequence_A[i] = unur_urng_sample(urng);

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_double_sequence_urng() */

/*...........................................................................*/

int compare_double_sequence_par( FILE *LOG, int line, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;
  int ok = TRUE;
  double x = 0.;
  int failed = 0;

  /* init generator */
  reset_sequence_par_warning(unur_get_default_urng(),par);
  gen = unur_init( par ); abort_if_NULL(LOG, line,gen);
  
  /* compare sequence */
  for (i=0; i<sample_size; i++) {
    x = unur_sample_cont(gen);	
    if ( !_unur_FP_equal(double_sequence_A[i], x)) {
      ok = FALSE;
      break;
    }
  }

  /* free generator */
  unur_free(gen); 
  
  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
    fprintf(LOG,"\tx[1] = %g, x[2] = %g, diff = %g\n",double_sequence_A[i],x,double_sequence_A[i]-x);
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_double_sequence_par() */

/*...........................................................................*/

int compare_double_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;

  /* check generator object */
  if (gen==NULL) {
    /* error */
    if (double_sequence_A) free (double_sequence_A);
    double_sequence_A = NULL;
    return UNUR_FAILURE;
  }

  /* allocate memory for storing sequence */
  if (double_sequence_A == NULL) {
    double_sequence_A = malloc( sample_size * sizeof(double) );
    abort_if_NULL(LOG,line, double_sequence_A);
  }

  /* reset uniform RNG */
  reset_sequence(unur_get_urng(gen));
  
  /* generate sequence */
  for (i=0; i<sample_size; i++)
    double_sequence_A[i] = unur_sample_cont(gen);
  
  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_double_sequence_gen_start() */

/*...........................................................................*/

int compare_double_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;
  int ok = TRUE;
  double x = 0.;
  int failed = 0;

  /* check generator object and stored sequence */
  if (gen==NULL || double_sequence_A==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* reset uniform RNG */
  reset_sequence_warning(unur_get_urng(gen));
  
  /* compare sequence */
  for (i=0; i<sample_size; i++) {
    x = unur_sample_cont(gen);	
    if ( !_unur_FP_equal(double_sequence_A[i], x)) {
      ok = FALSE;
      break;
    }
  }

  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
    fprintf(LOG,"\tx[1] = %g, x[2] = %g, diff = %g\n",double_sequence_A[i],x,double_sequence_A[i]-x);
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_double_sequence_gen() */

/*---------------------------------------------------------------------------*/
/* compare int sequences generated by generator    */
/* only when when we can reset the uniform RNG     */

static int *int_sequence_A = NULL;

int compare_int_sequence_par_start( FILE *LOG, int line, UNUR_PAR *par, int sample_size )
{
  int i;
  UNUR_GEN *gen;

  /* allocate memory for storing sequence */
  if (int_sequence_A == NULL) {
    int_sequence_A = malloc( sample_size * sizeof(int) );
    abort_if_NULL(LOG, line, int_sequence_A);
  }
  
  /* generate sequence */
  reset_sequence_par(unur_get_default_urng(),par);
  gen = unur_init( par ); abort_if_NULL(LOG, line,gen);

  for (i=0; i<sample_size; i++)
    int_sequence_A[i] = unur_sample_discr(gen);

  unur_free(gen); 

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_int_sequence_par_start() */

/*...........................................................................*/

int compare_int_sequence_par( FILE *LOG, int line, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;
  int ok = TRUE;
  int failed = 0;

  /* init generator */
  reset_sequence_par_warning(unur_get_default_urng(),par);
  gen = unur_init( par ); abort_if_NULL(LOG, line,gen);
  
  /* compare sequence */
  for (i=0; i<sample_size; i++)
    if (int_sequence_A[i] != unur_sample_discr(gen)) {
      ok = FALSE;
      break;
    }
  
  /* free generator */
  unur_free(gen); 
  
  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_int_sequence_par() */

/*...........................................................................*/

int compare_int_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;

  /* check generator object */
  if (gen==NULL) {
    /* error */
    if (int_sequence_A) free (int_sequence_A);
    int_sequence_A = NULL;
    return UNUR_FAILURE;
  }

  /* allocate memory for storing sequence */
  if (int_sequence_A == NULL) {
    int_sequence_A = malloc( sample_size * sizeof(int) );
    abort_if_NULL(LOG, line, int_sequence_A);
  }
  
  /* generate sequence */
  reset_sequence(unur_get_urng(gen));
  for (i=0; i<sample_size; i++)
    int_sequence_A[i] = unur_sample_discr(gen);

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_int_sequence_gen_start() */

/*...........................................................................*/

int compare_int_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;
  int ok = TRUE;
  int failed = 0;

  /* check generator object and stored sequence */
  if (gen==NULL || int_sequence_A==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* reset uniform RNG */
  reset_sequence_warning(unur_get_urng(gen));
  
  /* compare sequence */
  for (i=0; i<sample_size; i++)
    if (int_sequence_A[i] != unur_sample_discr(gen)) {
      ok = FALSE;
      break;
    }
  
  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_int_sequence_gen() */

/*---------------------------------------------------------------------------*/
/* compare sequences of double vectors generated by generator */
/* only when when we can reset the uniform RNG                */

static double *cvec_sequence_A = NULL;

int compare_cvec_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;
  int dim;

  /* free sequence */
  if (cvec_sequence_A) {
    free (cvec_sequence_A);
    cvec_sequence_A = NULL;
  }

  /* check generator object */
  if (gen==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* get dimension */
  dim = unur_get_dimension (gen);

  /* allocate memory for storing sequence */
  if (cvec_sequence_A == NULL) {
    cvec_sequence_A = malloc( dim * sample_size * sizeof(double) );
    abort_if_NULL(LOG,line, cvec_sequence_A);
  }

  /* reset uniform RNG */
  reset_sequence(unur_get_urng(gen));

  /* generate sequence */
  for (i=0; i<sample_size; i++)
    unur_sample_vec( gen, cvec_sequence_A+(i*dim) );

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_cvec_sequence_gen_start() */

/*...........................................................................*/

int compare_cvec_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i,k;
  int ok = TRUE;
  double *x;
  int failed = 0;
  int dim;

  /* check generator object and stored sequence */
  if (gen==NULL || cvec_sequence_A==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* get dimension */
  dim = unur_get_dimension (gen);

  /* need space for sampling result */
  x = malloc( dim * sizeof(double) );
  abort_if_NULL(LOG,line, x);

  /* reset uniform RNG */
  reset_sequence_warning(unur_get_urng(gen));
  
  /* compare sequence */
  for (i=0; i<sample_size; i++) {
    unur_sample_vec( gen, x );
    for (k=0; k<dim; k++) {
      if ( !_unur_FP_equal(cvec_sequence_A[i*dim+k], x[k])) {
	ok = FALSE;
	break;
      }
    }
    if (!ok) break;
  }

  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
    fprintf(LOG,"\tx[1] = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",cvec_sequence_A[i*dim+k]);
    fprintf(LOG,"),\tx[2] = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",x[k]);
    fprintf(LOG,"),\tdiff = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",cvec_sequence_A[i*dim+k]-x[k]);
    fprintf(LOG,")\n");
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  free (x);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_cvec_sequence_gen() */

/*---------------------------------------------------------------------------*/
/* compare sequences of double matrices generated by generator */
/* only when when we can reset the uniform RNG                */

static double *matr_sequence_A = NULL;

int compare_matr_sequence_gen_start( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i;
  int dim;

  /* free sequence */
  if (matr_sequence_A) {
    free (matr_sequence_A);
    matr_sequence_A = NULL;
  }

  /* check generator object */
  if (gen==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* get dimension */
  dim = unur_get_dimension (gen);

  /* allocate memory for storing sequence */
  if (matr_sequence_A == NULL) {
    matr_sequence_A = malloc( dim * sample_size * sizeof(double) );
    abort_if_NULL(LOG,line, matr_sequence_A);
  }

  /* reset uniform RNG */
  reset_sequence(unur_get_urng(gen));
  
  /* generate sequence */
  for (i=0; i<sample_size; i++)
    unur_sample_matr( gen, matr_sequence_A+(i*dim) );

  /* there cannot be a failure */
  return UNUR_SUCCESS;

} /* end of compare_matr_sequence_gen_start() */

/*...........................................................................*/

int compare_matr_sequence_gen( FILE *LOG, int line, UNUR_GEN *gen, int sample_size )
{
  int i,k;
  int ok = TRUE;
  double *x;
  int failed = 0;
  int dim;

  /* check generator object and stored sequence */
  if (gen==NULL || matr_sequence_A==NULL) {
    /* error */
    return UNUR_FAILURE;
  }

  /* get dimension */
  dim = unur_get_dimension (gen);

  /* need space for sampling result */
  x = malloc( dim * sizeof(double) );
  abort_if_NULL(LOG,line, x);

  /* reset uniform RNG */
  reset_sequence_warning(unur_get_urng(gen));
  
  /* compare sequence */
  for (i=0; i<sample_size; i++) {
    unur_sample_matr( gen, x );
    for (k=0; k<dim; k++) {
      if ( !_unur_FP_equal(matr_sequence_A[i*dim+k], x[k])) {
	ok = FALSE;
	break;
      }
    }
    if (!ok) break;
  }

  /* print result */
  fprintf(LOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
    fprintf(LOG,"\tx[1] = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",matr_sequence_A[i*dim+k]);
    fprintf(LOG,"),\tx[2] = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",x[k]);
    fprintf(LOG,"),\tdiff = (");
    for (k=0; k<dim; k++) 
      fprintf(LOG," %g",matr_sequence_A[i*dim+k]-x[k]);
    fprintf(LOG,")\n");
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  free (x);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of compare_matr_sequence_gen() */

/*---------------------------------------------------------------------------*/
/* free memory used for comparing sequences                                  */

void compare_free_memory( void )
{
  if (double_sequence_A) free (double_sequence_A);
  if (int_sequence_A)    free (int_sequence_A);
  if (cvec_sequence_A)   free (cvec_sequence_A);
  if (matr_sequence_A)   free (matr_sequence_A);
} /* end of compare_free_memory */

/*---------------------------------------------------------------------------*/
/* print name of distribution */

void print_distr_name( FILE *LOG, const UNUR_DISTR *distr, const char *genid )
{
  int i,n_fpar;
  const double *fpar;

  /* print name of distribution */
  if (genid)
    fprintf(LOG,"%s: ",genid);
  fprintf(LOG,"%s ",unur_distr_get_name(distr));

  /* get parameters */
  n_fpar = 0;
  if ( unur_distr_is_cont(distr) )
    /* continuous distribution */
    n_fpar = unur_distr_cont_get_pdfparams( distr, &fpar );
  else if ( unur_distr_is_discr(distr) )
    /* discrete distribution */
    n_fpar = unur_distr_discr_get_pmfparams( distr, &fpar );

  /* print parameter list */
  fprintf(LOG,"(");
  if (n_fpar) {
    fprintf(LOG,"%g",fpar[0]);
    for (i=1; i<n_fpar; i++)
      fprintf(LOG,", %g",fpar[i]);
  }
  fprintf(LOG,")");

} /* end of print_distr_name() */

/*---------------------------------------------------------------------------*/
/* check p-value of statistical test and print result */

int print_pval( FILE *LOG, UNUR_GEN *gen, const UNUR_DISTR *distr, double pval, int trial, char todo )
{
  int failed = 0;
  double pval_corrected;
  int dim;
  int l;

  /* Correct p-value.   */
  /* For multivariate distributions unur_test_chi2() runs chi^2 tests on the   */
  /* multivariate distribution itself as well as on all marginal distibutions. */
  /* It returns the minimum of all p-palues. Thus we have to adjust this       */
  /* value and multiply it with (dim+1) before deciding about the significance */
  /* of the result.                                                            */
  dim = unur_distr_get_dim(distr);
  pval_corrected = (dim>1) ? pval : pval*(dim+1);


  if (pval < -1.5) {
    /* was not able to run test (CDF missing) */

    fprintf(LOG,"   not performed (missing data)\t");
    /* print distribution name */
    print_distr_name( LOG, distr, gen ? unur_get_genid(gen):"???\t" );
    fprintf(LOG,"\n");

    printf("X");

    fflush(stdout);
    fflush(LOG);

    /* test does not count as failed */
    return UNUR_SUCCESS; 
    
  }

  if (pval < 0.) {
    fprintf(LOG,"   setup failed\t\t");
  }
  else {
    fprintf(LOG,"   pval = %8.6f   ",pval);
    l = _unur_isnan(pval_corrected) 
      ? 10000
      : -(int) ((pval_corrected > 1e-6) ? (log(pval_corrected) / M_LN10) : 6.);
    switch (l) {
    case 0:
      fprintf(LOG,"      "); break;
    case 1:
      fprintf(LOG,".     "); break;
    case 2:
      fprintf(LOG,"**    "); break;
    case 3:
      fprintf(LOG,"XXX   "); break;
    case 4:
      fprintf(LOG,"XXXX  "); break;
    case 5:
      fprintf(LOG,"XXXXX "); break;
    default:
      fprintf(LOG,"######"); break;
    }
  }
  
  switch (todo) {
  case '+':
    if (pval_corrected >= PVAL_LIMIT) {
      fprintf(LOG,"\t ok");
      printf("+");
    }
    else {
      failed = 1;
      if (trial > 1) {
	fprintf(LOG,"\t Failed");
	printf("(!+)");
      }
      else {
	fprintf(LOG,"\t Try again");
	printf("?");
      }
    }
    break;
  case '0':
  case '/':
    fprintf(LOG,"\t ok (expected to fail)");
    printf("0");
    break;
  case '-':
    if (pval_corrected < PVAL_LIMIT) {
      /* in this case it is expected to fail */
      failed = 0;
      fprintf(LOG,"\t ok (expected to fail)");
      printf("-");
    }
    else {
      /* the test has failed which was not what we have expected */
      failed = 1;
      fprintf(LOG,"\t Not ok (expected to fail)");
      printf("(!-)");
    }
    break;
  default:
    fprintf(stderr, "%s:%d: invalid test symbol\n", __FILE__, __LINE__);
    exit (-1);
  }

  /* print distribution name */
  fprintf(LOG,"\t");
  print_distr_name( LOG, distr, gen ? unur_get_genid(gen):"???\t" );
  fprintf(LOG,"\n");

  fflush(stdout);
  fflush(LOG);
  return (failed ? UNUR_FAILURE : UNUR_SUCCESS);

} /* end of print_pval() */

/*---------------------------------------------------------------------------*/
/* run chi2 test */

int run_validate_chi2( FILE *LOG, int line, UNUR_GEN *gen, const UNUR_DISTR *distr, int testtype, char todo )
     /*   UNUR_SUCCESS    ... on success                                        */
     /*   UNUR_ERR_SILENT ... test failed only once                             */
     /*   UNUR_FAILURE    ... serious failure                                   */
{
#define BUFSIZE 128
  const char *distr_name;
  static char last_distr_name[BUFSIZE] = "";
  unsigned int type;
  int i;
  double pval;
  int failed = 0;

  /* get name of distribution */
  distr_name = unur_distr_get_name( distr );

  /* get type of distribution */
  type = unur_distr_get_type( distr );

  if ( strcmp(distr_name,last_distr_name) ) {
    /* different distributions */
    strncpy(last_distr_name,distr_name,BUFSIZE);
    last_distr_name[BUFSIZE-1] = '\0';
    printf(" %s",distr_name); fflush(stdout);
  }

  if (todo == '.') {
    /* nothing to do */
    printf(".");  fflush(stdout);
    return UNUR_SUCCESS;
  }

  if (todo == '0') {
    /* initialization of generator is expected to fail */
    if (gen == NULL) {
      printf("0");  fflush(stdout);
      return UNUR_SUCCESS;
    }
    else {
      /* error */
      printf("(!0)");  fflush(stdout);
      return UNUR_FAILURE;
    }
  }

  if (gen == NULL) {
    if (todo == '-') {
      print_pval(LOG,gen,distr,-0.5,10,todo);
      return UNUR_SUCCESS;
    }
    else if (todo == '/') {
      printf("/");  fflush(stdout);
      return UNUR_SUCCESS;
    }
    else {
      /* initialization failed --> cannot run test */
      print_pval(LOG,gen,distr,-0.5,10,todo);
      return UNUR_FAILURE;
    }
  }

  /* init successful */
  if ( todo == '/' ) todo = '+';

  /* run chi^2 test */
  for (i=1; i<=2; i++) {
    /* we run the test twice when it fails the first time */

    switch (type) {
    case UNUR_DISTR_DISCR:
      pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 100000, 20, CHI_TEST_VERBOSITY, LOG);
      break;
    case UNUR_DISTR_CONT:
    case UNUR_DISTR_CEMP:
      pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 0, 20, CHI_TEST_VERBOSITY, LOG);
      break;
    case UNUR_DISTR_CVEC:
      switch (testtype) {
      case 1:
	pval = unur_test_chi2_marginal( gen, CHI_TEST_INTERVALS, 0, 20, CHI_TEST_VERBOSITY, LOG);
	break;
      default:
	pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 0, 20, CHI_TEST_VERBOSITY, LOG);
	break;
      }
      break;
    default:
      fprintf(stderr, "\n%s:%d: run_validate_chi2: this should not happen\n", __FILE__, __LINE__);
      exit (EXIT_FAILURE);
    }

    if ( print_pval(LOG,gen,distr,pval,i,todo) != UNUR_SUCCESS )
      /* test failed */
      failed++;
    else
      /* test ok */
      break;
  }

  return (failed==0 ? UNUR_SUCCESS : (failed<=1 ? UNUR_ERR_SILENT : UNUR_FAILURE));

#undef BUFSIZE
} /* end of run_validate_chi2() */

/*---------------------------------------------------------------------------*/
/* run verify hat test */

#define VERIFYHAT_SAMPLESIZE 10000

int run_validate_verifyhat( FILE *LOG, int line, UNUR_GEN *gen, const UNUR_DISTR *distr, char todo )
{
#define BUFSIZE 128
  const char *distr_name;
  static char last_distr_name[BUFSIZE] = "";
  unsigned int type;
  int i;
  int failed = 0;
  int dim;
  double *x = NULL;

  /* get name of distribution */
  distr_name = unur_distr_get_name( distr );

  /* get type of distribution */
  type = unur_distr_get_type( distr );

  if (strcmp(distr_name,last_distr_name) ) {
    /* different distributions */
    strncpy(last_distr_name,distr_name,BUFSIZE);
    last_distr_name[BUFSIZE-1] = '\0';
    printf(" %s",distr_name); fflush(stdout);
  }

  if (todo == '.') {
    /* nothing to do */
    printf(".");  fflush(stdout);
    return UNUR_SUCCESS;
  }

  if (todo == '0') {
    /* initialization of generator is expected to fail */
    if (gen == NULL) {
      printf("0");  fflush(stdout);
      print_verifyhat_result(LOG,gen,distr,-1,todo);
      return UNUR_SUCCESS;
    }
    else {
      /* error */
      printf("(!0)");  fflush(stdout);
      return UNUR_FAILURE;
    }
  }

  if (gen == NULL) {
    if (todo == '-') {
      printf("-");  fflush(stdout);
      print_verifyhat_result(LOG,gen,distr,-1,todo);
      return UNUR_SUCCESS;
    }
    else {
      /* initialization failed --> cannot run test */
      printf("(!+)");  fflush(stdout);
      print_verifyhat_result(LOG,gen,distr,-1,todo);
      return UNUR_FAILURE;
    }
  }

  /* get dimension of distribution */
  dim = unur_get_dimension (gen);
  
  /* allocate working space */
  if (dim > 0) {
    x = malloc( dim * sizeof(double) );
    abort_if_NULL(LOG, line, x);
  }

  /* run verify hat test */
  for (i=0; i<VERIFYHAT_SAMPLESIZE; i++) {

    unur_errno = 0;
    switch (type) {
    case UNUR_DISTR_DISCR:
      unur_sample_discr(gen);
      break;
    case UNUR_DISTR_CONT:
    case UNUR_DISTR_CEMP:
      unur_sample_cont(gen);
      break;
    case UNUR_DISTR_CVEC:
      unur_sample_vec(gen,x);
      break;
    default:
      fprintf(stderr, "\n%s:%d: run_validate_verifyhat: this should not happen\n", __FILE__, __LINE__);
      exit (EXIT_FAILURE);
    }

    if (unur_errno) {
      /* error */
      if (unur_errno == UNUR_ERR_GEN_CONDITION)
	failed++;
      unur_errno = 0;
    }
    
  }
  
  /* free working space */
  if (x!=NULL) free(x);

  return print_verifyhat_result(LOG,gen,distr,failed,todo);
  
} /* end of run_validate_verifyhat() */

/*---------------------------------------------------------------------------*/
/* print result of verify hat test */

int print_verifyhat_result( FILE *LOG, UNUR_GEN *gen, const UNUR_DISTR *distr, int failed, char todo )
{
  int failed_test = 0;
  double failed_ratio = ((double)failed) / VERIFYHAT_SAMPLESIZE;

  if (failed >= 0) {
    fprintf(LOG,"   failures = %d (%g%%)  ",failed, 100. * failed_ratio);
    switch (todo) {
    case '+':
      if (failed > 0) {
	fprintf(LOG,"\t Failed");
	printf("(!+)");
	failed_test = 1;
      }
      else {
	fprintf(LOG,"\t ok");
	printf("+");
      }
      break;
    case '~':
      if (failed == 0) {
	fprintf(LOG,"\t ok");
	printf("+");
      }
      else if (failed_ratio <= 0.01) {
	fprintf(LOG,"\t tolerated");
	printf("(~+)");
      }
      else {
	fprintf(LOG,"\t Failed");
	printf("(!~+)");
	failed_test = 1;
      }
      break;
    case '-':
      if (failed_ratio > 0.01) {
	/* in this case it is expected to fail */
	fprintf(LOG,"\t ok (expected to fail)");
	printf("-");
      }
      else {
	/* the test has failed which was not what we have expected */
	fprintf(LOG,"\t Not ok (expected to fail)");
	printf("(!-)");
	failed_test = 1;
      }
      break;
    default:
      fprintf(stderr,"%s:%d: invalid test symbol\n", __FILE__, __LINE__);
      exit (EXIT_FAILURE);
    }
  }

  else {
    fprintf(LOG,"   setup failed\t");
    switch (todo) {
    case '0':
    case '-':
      fprintf(LOG,"\t ok (expected to fail)");
      break;
    default:
      fprintf(LOG,"\t Failed");
    }
  }

  /* print distribution name */
  fprintf(LOG,"\t");
  print_distr_name( LOG, distr, gen ? unur_get_genid(gen):"???\t" );
  fprintf(LOG,"\n");

  fflush(stdout);
  fflush(LOG);

  return (failed_test ? UNUR_ERR_SILENT : UNUR_SUCCESS);

} /* end of print_verifyhat_result() */
#undef VERIFYHAT_SAMPLESIZE

/*---------------------------------------------------------------------------*/
/* print result of timings */

void print_timing_results( FILE *LOG, int line, const UNUR_DISTR *distr,
			   double *timing_setup, double *timing_marginal, int n_results )
{
  const char *distr_name;
  static const char *last_distr_name = "";
  int i;

  /* get name of distribution */
  distr_name = unur_distr_get_name( distr );

  if (strcmp(distr_name,last_distr_name) ) {
    /* different distributions */
    last_distr_name = distr_name;
    printf(" %s=",distr_name); fflush(stdout);
  }
  else {
    printf("="); fflush(stdout);
  }
  
  /* print timings into log file */
  for (i=0; i<n_results; i++)
    if (timing_marginal[i] < 0.)
      /* no test */
      fprintf(LOG, "      /        ");
    else
      fprintf(LOG, "%5.0f /%7.2f ", timing_setup[i], timing_marginal[i]);

  /* print name of distribution into log file */
  fprintf(LOG,"\t");
  print_distr_name( LOG, distr, NULL );
  fprintf(LOG,"\n");

} /* end of print_timing_results() */

/*---------------------------------------------------------------------------*/
