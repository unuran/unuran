/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Common test routines                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"
#include <time.h>

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

int test_ok = TRUE;               /* all tests ok (boolean)                  */
int test_failed = 0;              /* failed tests                            */
FILE *TESTLOG = NULL;             /* test log file                           */

/*---------------------------------------------------------------------------*/

static FILE *UNURANLOG;      /* unuran log file */

/* open log files for testing */
void open_log_files( const char *file )
{
  char filename[128];
  time_t started;  

  /* open test log file */
  sprintf(filename,"%s_test.log",file);
  TESTLOG = fopen(filename,"w");
  if (TESTLOG == NULL) exit (-1);

  /* write header into log file */
  fprintf(TESTLOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");
  if (time( &started ) != -1)
    fprintf(TESTLOG,"%s",ctime(&started));
  fprintf(TESTLOG,"\n====================================================\n\n");

  /* set output stream for unuran messages */
  sprintf(filename,"%s_unuran.log",file);
  UNURANLOG = fopen(filename,"w");
  if (UNURANLOG == NULL) exit (-1);
  unur_set_stream( UNURANLOG );

} /* end of open_log_files() */

/*---------------------------------------------------------------------------*/

/* close log files */
void close_log_files( void )
{
  fprintf(TESTLOG,"\n====================================================\n\n");
  if (test_ok)
    fprintf(TESTLOG,"All tests PASSED.\n");
  else
    fprintf(TESTLOG,"Test(s) FAILED.\n");

  fclose(UNURANLOG);
  fclose(TESTLOG);
} /* end of close_log_files() */

/*---------------------------------------------------------------------------*/

/* check for invalid NULL pointer, that should not happen in this program */
void do_abort_if_NULL( int line, void *ptr )
{
  if (ptr) return; /* o.k. */

  /* There must not be a NULL pointer.
     Since we do not expect a NULL pointer something serious has
     happend. So it is better to abort the tests.
  */
  fprintf(TESTLOG,"line %4d: Unexpected NULL pointer. Panik --> Abort tests!!\t\t",line);
  printf(" Panik --> Tests aborted\n"); fflush(stdout);

  test_ok = FALSE;
  
  /* test finished */
  printf("\n");  fflush(stdout);

  /* close log files and exit */
  close_log_files();
  exit(-1);
} /* abort_if_NULL() */

/*---------------------------------------------------------------------------*/

/* compare error code */
void do_check_errorcode( int line, unsigned errno )
{
  fprintf(TESTLOG,"line %4d: Error code ...\t\t",line);
  if (unur_errno != errno) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
} /* end of do_check_errorcode() */

/*---------------------------------------------------------------------------*/

/* check for expected NULL pointer */
void do_check_expected_NULL( int line, void *ptr )
{
  fprintf(TESTLOG,"line %4d: NULL pointer expected ...\t",line);
  if (ptr != NULL) { 
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
} /* end of do_check_expected_NULL() */

/*---------------------------------------------------------------------------*/

/* check for "set failed" */
void do_check_expected_setfailed( int line, int ok )
{
  fprintf(TESTLOG,"line %4d: `failed' expected ...\t",line);
  if (ok) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
} /* end of do_check_expected_setfailed() */

/*---------------------------------------------------------------------------*/

#define double_not_equal(a,b) \
 ( ((a) > 0. && ((a) > (1.+DBL_EPSILON) * (b) || \
                 (a) < (1.-DBL_EPSILON) * (b)) ) || \
   ((a) < 0. && ((a) < (1.+DBL_EPSILON) * (b) || \
                 (a) > (1.-DBL_EPSILON) * (b)) ) )

/* compare two sequences of type double */
void do_compare_double_sequences( int line, double *a, double *b, int n )
{
  int i;
  int ok = TRUE;

  for (i=0; i<n; i++)
    if (double_not_equal(a[i],b[i])) {
      ok = FALSE;
      break;
    }
  fprintf(TESTLOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
  
} /* end of do_compare_double_sequences() */

#undef double_not_equal(a,b)

/*---------------------------------------------------------------------------*/

/* compare two sequences of type int */
void do_compare_int_sequences( int line, int *a, int *b, int n )
{
  int i;
  int ok = TRUE;

  for (i=0; i<n; i++)
    if (a[i] != b[i]) {
      ok = FALSE;
      break;
    }
  fprintf(TESTLOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
  
} /* end of do_compare_int_sequences() */

/*---------------------------------------------------------------------------*/

/* check p-value of statistical test */
void do_check_pval( int line, UNUR_GEN *gen, double pval, int trial )
{
  fprintf(TESTLOG,"line %4d: ",line);
  print_pval(pval,trial);

  /* print distribution name */
  fprintf(TESTLOG,"\t");
  print_distr_name( unur_get_distr(gen), unur_get_genid(gen) );
  fprintf(TESTLOG,"\n");

} /* end of do_check_pval() */
  
/*---------------------------------------------------------------------------*/

/* print name of distribution */
void print_distr_name( UNUR_DISTR *distr, const char *genid )
{
  int i,n_fpar;
  double *fpar;

  fprintf(TESTLOG,"%s: %s ",genid,unur_distr_get_name(distr));

  if ( unur_distr_is_cont(distr) ) {
    n_fpar = unur_distr_cont_get_pdfparams( distr, &fpar );
    fprintf(TESTLOG,"(");
    if (n_fpar) {
      fprintf(TESTLOG,"(%g",fpar[0]);
      for (i=1; i<n_fpar; i++)
	fprintf(TESTLOG,", %g",fpar[i]);
    }
    fprintf(TESTLOG,")");
  }

} /* end of print_distr_name() */

/*---------------------------------------------------------------------------*/

/* print p-value of statistical test */
void print_pval( double pval, int trial )
{
  int l;

  if (pval < 0.) {
    /* test has not been executed */
    fprintf(TESTLOG,"pval = (0)\t      Not performed");
    printf(".");
    return;
  }

  fprintf(TESTLOG,"pval = %8.6f   ",pval);

  l = -(int) ((pval > 1e-6) ? (log(pval) / M_LN10) : 6.);

  switch (l) {
  case 0:
    fprintf(TESTLOG,"      "); break;
  case 1:
    fprintf(TESTLOG,".     "); break;
  case 2:
    fprintf(TESTLOG,"**    "); break;
  case 3:
    fprintf(TESTLOG,"XXX   "); break;
  case 4:
    fprintf(TESTLOG,"XXXX  "); break;
  case 5:
    fprintf(TESTLOG,"XXXXX "); break;
  default:
    fprintf(TESTLOG,"######"); break;
  }
    
  if (pval < PVAL_LIMIT) {
    ++test_failed;
    if (trial > 1) {
      fprintf(TESTLOG,"\t Failed");
      printf("X");
    }
    else {
      fprintf(TESTLOG,"\t Try again");
      printf("?");
    }
  }
  else {
    fprintf(TESTLOG,"\t ok");
    printf("+");
  }
  fflush(stdout);

} /* end of print_pval() */
  
/*---------------------------------------------------------------------------*/

/* run level 2 test on collected p-values */
void run_level2( int line, double *pvals, int n_pvals )
{
  int i;
  int *classes;
  int n_classes;
  double pval2;

  if (pvals==NULL)
    /* nothing to do */
    return;

  /* number classes */
  n_classes = (int) (sqrt(n_pvals)+0.5);
  if (n_pvals/n_classes < 6)
    /* classes would have too few entries */
    n_classes = n_pvals / 6;

  /* allocate memory for classes */
  classes = calloc( n_classes+1, sizeof(int) );

  /* count bins */
  for (i=0; i<n_pvals; i++)
    ++(classes[ (int)(pvals[i] * n_classes) ]);

  /* run test */
  pval2 = _unur_test_chi2test( NULL, classes, n_classes, 5, 0 );

  /* print result */
  printf(" Level-2-test(%d)",n_pvals);
  fprintf(TESTLOG,"line %4d: ",__LINE__);
  print_pval(pval2,100);
  fprintf(TESTLOG,"\tLevel 2 Test (n=%d) (not powerful)\n",n_pvals);

  /* clear */
  free(classes);

} /* end of run_level2() */

/*---------------------------------------------------------------------------*/

