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
  int i,l,n_fpar;
  double *fpar;
  UNUR_DISTR *distr;

  l = -(int) ((pval > 1e-6) ? (log(pval) / M_LN10) : 6.);

  fprintf(TESTLOG,"line %4d: pval = %8.6f   ",line,pval);

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
