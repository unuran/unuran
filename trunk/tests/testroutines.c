/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Common test routines                                                     *
 *                                                                           *
 *****************************************************************************
     $Id$ 
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "testunuran.h"

/*---------------------------------------------------------------------------*/

#define _unur_FP_equal(a,b) \
 ((a)==(b) || \
 fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b)) * 100. * DBL_EPSILON)

/*---------------------------------------------------------------------------*/
/* check for invalid NULL pointer, that should not happen in this program */

void abort_if_NULL( FILE *LOG, int line, void *ptr )
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
  exit(-1);
  
} /* abort_if_NULL() */

/*---------------------------------------------------------------------------*/
/* compare error code */

int check_errorcode( FILE *LOG, int line, unsigned errno )
{
  int failed = 0;

  fprintf(LOG,"line %4d: Error code ...\t\t",line);
  if (unur_errno != errno) {
    failed = 1;
    fprintf(LOG," Failed");
    fprintf(LOG," (observed = %#x, expected = %#x)\n",unur_errno,errno);
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_errorcode() */

/*---------------------------------------------------------------------------*/
/* check for expected NULL pointer */

int check_expected_NULL( FILE *LOG, int line, void *ptr )
{
  int failed = 0;

  fprintf(LOG,"line %4d: NULL pointer expected ...\t",line);
  if (ptr != NULL) { 
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_expected_NULL() */

/*---------------------------------------------------------------------------*/
/* check for "set failed" */

int check_expected_setfailed( FILE *LOG, int line, int ok )
{
  int failed = 0;

  fprintf(LOG,"line %4d: `failed' expected ...\t",line);
  if (ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_expected_setfailed() */

/*---------------------------------------------------------------------------*/
/* check for INFINITY */

int check_expected_INFINITY( FILE *LOG, int line, double x )
{
  int failed = 0;

  fprintf(LOG,"line %4d: INFINITY expected ...\t",line);
  if (x < UNUR_INFINITY) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_expected_INFINITY() */

/*---------------------------------------------------------------------------*/
/* check for reinit */

int check_expected_reinit( FILE *LOG, int line, int ok )
{
  int failed = 0;

  fprintf(LOG,"line %4d: reinit ...\t\t\t",line);
  if (!ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_expected_reinit() */

/*---------------------------------------------------------------------------*/
/* check for non existing reinit */

int check_expected_no_reinit( FILE *LOG, int line, int ok )
{
  int failed = 0;

  fprintf(LOG,"line %4d: no reinit ...\t\t",line);
  if (ok) {
    failed = 1;
    fprintf(LOG," Failed\n");
  }
  else
    fprintf(LOG," ok\n");

  fflush(LOG);
  return failed;

} /* end of check_expected_no_reinit() */

/*---------------------------------------------------------------------------*/
/* compare double sequences generated by generator */

static double *double_sequence_A = NULL;

int compare_double_sequence_par_start( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;

  if (double_sequence_A == NULL) {
    double_sequence_A = malloc( sample_size * sizeof(double) );
    abort_if_NULL(LOG,line, double_sequence_A);
  }

  /* generate sequence */
  prng_reset(urng);
  gen = unur_init( par ); abort_if_NULL(LOG,line,gen);
  
  for (i=0; i<sample_size; i++)
    double_sequence_A[i] = unur_sample_cont(gen);
  
  unur_free(gen); 

  /* there cannot be a failure */
  return 0;

} /* end of compare_double_sequence_par_start() */

/*...........................................................................*/

int compare_double_sequence_urng_start( FILE *LOG, int line, struct prng *urng, int sample_size )
{
  int i;
  
  if (double_sequence_A == NULL) {
    double_sequence_A = malloc( sample_size * sizeof(double) );
    abort_if_NULL(LOG, line, double_sequence_A);
  }

  /* generate sequence */
  prng_reset(urng);
  for (i=0; i<sample_size; i++)
    double_sequence_A[i] = prng_get_next(urng);

  /* there cannot be a failure */
  return 0;

} /* end of compare_double_sequence_urng() */

/*...........................................................................*/

int compare_double_sequence_par( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;
  int ok = TRUE;
  double x;
  int failed = 0;

  /* init generator */
  prng_reset(urng);
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
  }
  else
    fprintf(LOG," ok\n");
  
  fflush(LOG);
  return failed;

} /* end of compare_double_sequence_par() */

/*---------------------------------------------------------------------------*/
/* compare int sequences generated by generator */

static int *int_sequence_A = NULL;

int compare_int_sequence_par_start( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;

  if (int_sequence_A == NULL) {
    int_sequence_A = malloc( sample_size * sizeof(int) );
    abort_if_NULL(LOG, line, int_sequence_A);
  }
  
  /* generate sequence */
  prng_reset(urng);
  gen = unur_init( par ); abort_if_NULL(LOG, line,gen);

  for (i=0; i<sample_size; i++)
    int_sequence_A[i] = unur_sample_discr(gen);

  unur_free(gen); 

  /* there cannot be a failure */
  return 0;

} /* end of compare_int_sequence_par_start() */

/*...........................................................................*/

int compare_int_sequence_par( FILE *LOG, int line, struct prng *urng, UNUR_PAR *par, int sample_size )
{
  UNUR_GEN *gen;
  int i;
  int ok = TRUE;
  int failed = 0;

  /* init generator */
  prng_reset(urng);
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
  return failed;

} /* end of compare_int_sequence_par() */

/*---------------------------------------------------------------------------*/
/* print name of distribution */

void print_distr_name( FILE *LOG, UNUR_DISTR *distr, const char *genid )
{
  int i,n_fpar;
  double *fpar;

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

int print_pval( FILE *LOG, UNUR_GEN *gen, double pval, int trial, char todo )
{
  int failed = 0;
  int l;

  if (pval < -1.5) {
    /* was not able to run test (CDF missing) */

    fprintf(LOG,"   not performed (missing data)\t");
    /* print distribution name */
    print_distr_name( LOG,unur_get_distr(gen), unur_get_genid(gen) );
    fprintf(LOG,"\n");

    printf("X");

    fflush(stdout);
    fflush(LOG);

    /* test does not count as failed */
    return 0; 
    
  }

  fprintf(LOG,"   pval = %8.6f   ",pval);

  l = -(int) ((pval > 1e-6) ? (log(pval) / M_LN10) : 6.);

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
  
  switch (todo) {
  case '+':
    if (pval < PVAL_LIMIT) {
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
    else {
      fprintf(LOG,"\t ok");
      printf("+");
    }
    break;
  case '-':
    if (pval < PVAL_LIMIT) {
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
    fprintf(stderr,"invalid test symbol\n");
    exit (-1);
  }

  /* print distribution name */
  fprintf(LOG,"\t");
  print_distr_name( LOG,unur_get_distr(gen), unur_get_genid(gen) );
  fprintf(LOG,"\n");

  fflush(stdout);
  fflush(LOG);
  return failed;

} /* end of print_pval() */

/*---------------------------------------------------------------------------*/
/* run chi2 test */

int run_validate_chi2( FILE *LOG, int line, UNUR_GEN *gen, char todo )
{
  const char *distr_name;
  static const char *last_distr_name = "";
  unsigned int type;
  int i;
  double pval;
  int failed = 0;

  if (todo == '.') {
    /* nothing to do */
    printf(".");  fflush(stdout);
    return 0;
  }

  if (todo == '0') {
    /* initialization of generator is expected to fail */
    if (gen == NULL) {
      printf("0");  fflush(stdout);
      return 0;
    }
    else {
      /* error */
      printf("(!0)");  fflush(stdout);
      return 2;
    }
  }

  if (gen == NULL) {
    if (todo == '-') {
      printf("0");  fflush(stdout);
      return 0;
    }
    else if (todo == '/') {
      printf("/");  fflush(stdout);
      return 0;
    }
    else {
      /* initialization failed --> cannot run test */
      printf("(!+)");  fflush(stdout);
      return 2;
    }
  }

  /* init successful */
  if ( todo == '/' ) todo = '+';

  /* get name of distribution */
  distr_name = unur_distr_get_name( unur_get_distr(gen) );

  /* get type of distribution */
  type = unur_distr_get_type( unur_get_distr(gen) );

  if (strcmp(distr_name,last_distr_name) ) {
    /* different distributions */
    last_distr_name = distr_name;
    printf(" %s",distr_name); fflush(stdout);
  }

  /* run chi^2 test */
  for (i=1; i<=2; i++) {
    /* we run the test twice when it fails the first time */

    switch (type) {
    case UNUR_DISTR_CONT:
      pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 0, 20, CHI_TEST_VERBOSITY, LOG);
      break;
    case UNUR_DISTR_DISCR:
      pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 100000, 20, CHI_TEST_VERBOSITY, LOG);
      break;
    default:
      fprintf(stderr,"this should not happen\n");
      exit (-1);
    }

    if ( print_pval(LOG,gen,pval,i,todo) )
      /* test failed */
      failed++;
    else
      /* test ok */
      break;
  }

  return failed;

} /* end of run_validate_chi2() */

/*---------------------------------------------------------------------------*/
/* run verify hat test */

#define VERIFYHAT_SAMPLESIZE 10000

int run_validate_verifyhat( FILE *LOG, int line, UNUR_GEN *gen, char todo )
{
  const char *distr_name;
  static const char *last_distr_name = "";
  unsigned int type;
  int i;
  int failed = 0;

  if (todo == '.') {
    /* nothing to do */
    printf(".");  fflush(stdout);
    return 0;
  }

  if (todo == '0') {
    /* initialization of generator is expected to fail */
    if (gen == NULL) {
      printf("0");  fflush(stdout);
      return 0;
    }
    else {
      /* error */
      printf("(!0)");  fflush(stdout);
      return 2;
    }
  }

  if (gen == NULL) {
    if (todo == '-') {
      printf("0");  fflush(stdout);
      return 0;
    }
    else {
      /* initialization failed --> cannot run test */
      printf("(!+)");  fflush(stdout);
      return 2;
    }
  }

  /* get name of distribution */
  distr_name = unur_distr_get_name( unur_get_distr(gen) );

  /* get type of distribution */
  type = unur_distr_get_type( unur_get_distr(gen) );

  if (strcmp(distr_name,last_distr_name) ) {
    /* different distributions */
    last_distr_name = distr_name;
    printf(" %s",distr_name); fflush(stdout);
  }

  /* run verify hat test */
  for (i=0; i<VERIFYHAT_SAMPLESIZE; i++) {

    unur_errno = 0;
    switch (type) {
    case UNUR_DISTR_CONT:
      unur_sample_cont(gen);
      break;
    case UNUR_DISTR_DISCR:
      unur_sample_discr(gen);
      break;
    default:
      fprintf(stderr,"this should not happen\n");
      exit (-1);
    }

    if (unur_errno) {
      /* error */
      if (unur_errno == UNUR_ERR_GEN_CONDITION)
	failed++;
      unur_errno = 0;
    }
    
  }
  
  return print_verifyhat_result(LOG,gen,failed,todo);
  
} /* end of run_validate_verifyhat() */

/*---------------------------------------------------------------------------*/
/* print result of verify hat test */

int print_verifyhat_result( FILE *LOG, UNUR_GEN *gen, int failed, char todo )
{
  int failed_test = 0;
  double failed_ratio = ((double)failed) / VERIFYHAT_SAMPLESIZE;

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
    fprintf(stderr,"invalid test symbol\n");
    exit (-1);
  }

  /* print distribution name */
  fprintf(LOG,"\t");
  print_distr_name( LOG,unur_get_distr(gen), unur_get_genid(gen) );
  fprintf(LOG,"\n");

  fflush(stdout);
  fflush(LOG);

  return failed_test;

} /* end of print_verifyhat_result() */
#undef VERIFYHAT_SAMPLESIZE

/*---------------------------------------------------------------------------*/
/* print result of timings */

void print_timing_results( FILE *LOG, int line, UNUR_DISTR *distr, double *timing_result, int n_results )
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
    if (timing_result[i] < 0.)
      /* no test */
      fprintf(LOG, "    .  ");
    else
      fprintf(LOG, "%7.2f", timing_result[i]);

  /* print name of distribution into log file */
  fprintf(LOG,"\t");
  print_distr_name( LOG, distr, NULL );
  fprintf(LOG,"\n");

} /* end of print_timing_results() */

/*---------------------------------------------------------------------------*/
