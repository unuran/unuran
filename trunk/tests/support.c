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

#define double_equal(a,b) \
 ( ((a) > 0. && ((a) > (1.+DBL_EPSILON) * (b) || \
                 (a) < (1.-DBL_EPSILON) * (b)) ) || \
   ((a) < 0. && ((a) < (1.+DBL_EPSILON) * (b) || \
                 (a) > (1.-DBL_EPSILON) * (b)) ) )

/* compare two sequences */
void do_compare_sequences( int line, double *a, double *b, int n )
{
  int i;
  int ok = TRUE;

  for (i=0; i<n; i++)
    if (double_equal(a[i],b[i])) {
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
  
} /* end of do_compare_sequences() */

#undef double_equal(a,b)

/*---------------------------------------------------------------------------*/

/* check p-value of statistical test */
void do_check_pval( int line, UNUR_GEN *gen, double pval, int trial )
{
  int i,n_fpar;
  double *fpar;
  UNUR_DISTR *distr;

  fprintf(TESTLOG,"line %4d: ",line);
  print_pval(pval,trial);

  /* print distribution name */
  distr = unur_get_distr(gen);
  n_fpar = unur_distr_cont_get_pdfparams( distr, &fpar );
  fprintf(TESTLOG,"\t%s: %s (",unur_get_genid(gen),unur_distr_get_name(distr));
  if (n_fpar) {
    fprintf(TESTLOG,"%g",fpar[0]);
    for (i=1; i<n_fpar; i++)
      fprintf(TESTLOG,", %g",fpar[i]);
  }
  fprintf(TESTLOG,")\n");

} /* end of do_check_pval() */
  
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
  fprintf(TESTLOG,"\tLevel 2 Test (n=%d)\n",n_pvals);

  /* clear */
  free(classes);

} /* end of run_level2() */

/*---------------------------------------------------------------------------*/

