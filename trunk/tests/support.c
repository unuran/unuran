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

UNUR_DISTR **list_of_distr;         /* pointer to list of distributions      */
int n_distr;                        /* number of distributions               */

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

/* compare two sequences */
void do_compare_sequences( int line, double *a, double *b, int n )
{
  int i;
  int ok = TRUE;

  for (i=0; i<n; i++) {
    if ( (a[i] > 0. && a[i] > (1.+DBL_EPSILON) * b[i]) || 
	 (a[i] > 0. && a[i] < (1.-DBL_EPSILON) * b[i]) || 
	 (a[i] < 0. && a[i] < (1.+DBL_EPSILON) * b[i]) || 
	 (a[i] < 0. && a[i] > (1.-DBL_EPSILON) * b[i]) ) {
      ok = FALSE;
      break;
    }
  }
  fprintf(TESTLOG,"line %4d: random seqences ...\t\t",line);
  if (!ok) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");
  
} /* end of do_compare_sequences() */

/*---------------------------------------------------------------------------*/

/* check p-value of statistical test */
void check_pval( UNUR_GEN *gen, double pval )
{
  int i,l,n_fpar;
  double *fpar;
  UNUR_DISTR *distr;

  l = -(int) ((pval > 1e-6) ? (log(pval) / M_LN10) : 6.);

  fprintf(TESTLOG,"%s: pval = %8.6f   ",unur_get_genid(gen),pval);

  switch (l) {
  case 0:
    fprintf(TESTLOG,"      "); break;
  case 1:
    fprintf(TESTLOG,".     "); break;
  case 2:
    fprintf(TESTLOG,"**    "); break;
  case 3:
    fprintf(TESTLOG,"***   "); break;
  case 4:
    fprintf(TESTLOG,"****  "); break;
  case 5:
    fprintf(TESTLOG,"***** "); break;
  default:
    fprintf(TESTLOG,"######"); break;
  }
    
  if (pval < 1e-3) {
    ++test_failed;
    fprintf(TESTLOG,"\t Failed");
  }
  else
    fprintf(TESTLOG,"\t ok");

  /* print distribution name */
  distr = unur_get_distr(gen);
  n_fpar = unur_distr_cont_get_pdfparams( distr, &fpar );
  fprintf(TESTLOG,"\t%s (",unur_distr_get_name(distr));
  if (n_fpar) {
    fprintf(TESTLOG,"%g",fpar[0]);
    for (i=1; i<n_fpar; i++)
      fprintf(TESTLOG,", %g",fpar[i]);
  }
  fprintf(TESTLOG,")\n");

} /* end of do_check_pval() */
  

/*---------------------------------------------------------------------------*/

/* make list of distributions                                                */
void make_list_of_distributions( void )
{
  double fpar[10];

  /* reset counter for distributions */
  n_distr = 0;

  /* allocate memory for list */
  list_of_distr = malloc( N_DISTRIBUTIONS * sizeof(UNUR_DISTR *));

#ifdef D_GAMMA
  /** Gamma distributions **/
  fpar[0] = 1.;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,1);

  fpar[0] = 2.;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,1);

  fpar[0] = 3.;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,1);

  fpar[0] = 10.;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,1);

  fpar[0] = 1000.;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,1);

  fpar[0] = 5.;
  fpar[1] = 1.e+10;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,2);

  fpar[0] = 5.;
  fpar[1] = 1.e-10;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,2);

  fpar[0] = 5.;
  fpar[1] = 10.;
  fpar[2] = 1.e+10;
  list_of_distr[n_distr++] = unur_distr_gamma(fpar,3);
#endif

#ifdef D_NORMAL
  /** Normal distributions **/
  list_of_distr[n_distr++] = unur_distr_normal(NULL,0);

  fpar[0] = 10.;
  fpar[1] = 1.e-10;
  list_of_distr[n_distr++] = unur_distr_normal(fpar,2);

  fpar[0] = 0.;
  fpar[1] = 1.e+10;
  list_of_distr[n_distr++] = unur_distr_normal(fpar,2);
#endif

  /* check N_DISTRIBUTIONS (compile time constant !!) */
  if (n_distr > N_DISTRIBUTIONS) {
    /* this is very bad */
    fprintf(stderr,"N_DISTRIBUTIONS must be raised to %d\n",n_distr);
    abort();
  }

} /* end of make_list_of_distributions() */

/*---------------------------------------------------------------------------*/

