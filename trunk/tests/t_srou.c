/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests                                                                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/
#ifdef T_SROU
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

int test_ok = TRUE;               /* all tests ok (boolean)                  */
int test_failed = 0;              /* failed tests                            */
FILE *TESTLOG = NULL;             /* test log file                           */

static FILE *UNURANLOG = NULL;    /* unuran log file                         */

static struct prng *urng;         /* uniform random number generator         */

#define FAILED_LIMIT 1            /* maximal number of single failed 
				     statistical tests (if there are more 
				     the generator has failed the test.)     */

double *list_pvals = NULL;        /* list of collected p-values              */
int size_pvals = 0;               /* size of list                            */
int n_pvals = 0;                  /* number of collected p-values            */


/*---------------------------------------------------------------------------*/

void test_srou_new( void );
void test_srou_set( void );
void test_srou_chg( void );
void test_srou_init( void );
void test_srou_reinit( void );
void test_srou_sample( void );

void test_srou_validate( void );
void test_srou_validate_chi2( UNUR_DISTR *distr );
void run_chi2( UNUR_PAR *par, UNUR_DISTR *distr, int line );

double pdf( double x, UNUR_DISTR *distr );

/*---------------------------------------------------------------------------*/

int main()
{ 
  char filename[128];

  /* open test log file */
  sprintf(filename,"%s_test.log",__FILE__);
  TESTLOG = fopen(filename,"w");
  if (TESTLOG == NULL) exit (-1);

  /* set output stream for unuran messages */
  sprintf(filename,"%s_unuran.log",__FILE__);
  UNURANLOG = fopen(filename,"w");
  if (UNURANLOG == NULL) exit (-1);
  unur_set_stream( UNURANLOG );

  /* make a list of distributions */
  make_list_of_distributions();

  /* we use Mersenne Twister as uniform random number generator */
  urng = prng_new("mt19937(5678)");
  unur_set_default_urng(urng);

  /* start test */
  printf("%s: ",__FILE__);

  /* run tests */
  test_srou_new();
  test_srou_set();
  test_srou_chg();
  test_srou_init();
  test_srou_sample();
  test_srou_reinit();
  test_srou_validate();

  /* test finished */
  printf("\n");  fflush(stdout);

  /* close log files and exit */
  fclose(TESTLOG);
  fclose(UNURANLOG);
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

void test_srou_new( void )
{
  UNUR_DISTR *distr;

  /* start test */
  printf("[new "); fflush(stdout);
  fprintf(TESTLOG,"\n[new]\n");

  test_failed = 0;

  /* check error handling */

  /* error: NULL pointer */
  distr = NULL;   
  check_expected_NULL( unur_srou_new(distr) );
  check_errorcode( UNUR_ERR_NULL );
     
  /* error: invalid distribution type */
  distr = unur_distr_discr_new();
  check_expected_NULL( unur_srou_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_INVALID );
  unur_distr_free(distr);

  /* error: data missing in distribution object */
  distr = unur_distr_cont_new();               /* pdf, mode, pdfarea */
  check_expected_NULL( unur_srou_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_REQUIRED );
  unur_distr_cont_set_pdf(distr,pdf);          /* mode, pdfarea */
  check_expected_NULL( unur_srou_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_REQUIRED );
  unur_distr_cont_set_mode(distr,1.);          /* pdfarea */
  check_expected_NULL( unur_srou_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_REQUIRED );
  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_srou_new() */

/*---------------------------------------------------------------------------*/

void test_srou_set( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  double fpar[2] = {0.,1.};

  /* start test */
  printf("[set "); fflush(stdout);
  fprintf(TESTLOG,"\n[set]\n");

  test_failed = 0;

  /* check error handling */

  /* error: NULL pointer */
  par = NULL;   
  check_expected_setfailed( unur_srou_set_cdfatmode(par,0.) );
  check_errorcode( UNUR_ERR_NULL );
  check_expected_setfailed( unur_srou_set_pdfatmode(par,1.) );
  check_errorcode( UNUR_ERR_NULL );
  check_expected_setfailed( unur_srou_set_verify(par,1) );
  check_errorcode( UNUR_ERR_NULL );
  check_expected_setfailed( unur_srou_set_usesqueeze(par,1) );
  check_errorcode( UNUR_ERR_NULL );
  check_expected_setfailed( unur_srou_set_usemirror(par,1) );
  check_errorcode( UNUR_ERR_NULL );

  /* error: invalid parameter object */
  distr = unur_distr_normal(fpar,2);
  par = unur_arou_new(distr);
  check_expected_setfailed( unur_srou_set_cdfatmode(par,0.) );
  check_errorcode( UNUR_ERR_PAR_INVALID );
  check_expected_setfailed( unur_srou_set_pdfatmode(par,1.) );
  check_errorcode( UNUR_ERR_PAR_INVALID );
  check_expected_setfailed( unur_srou_set_verify(par,1) );
  check_errorcode( UNUR_ERR_PAR_INVALID );
  check_expected_setfailed( unur_srou_set_usesqueeze(par,1) );
  check_errorcode( UNUR_ERR_PAR_INVALID );
  check_expected_setfailed( unur_srou_set_usemirror(par,1) );
  check_errorcode( UNUR_ERR_PAR_INVALID );
  free(par); 
  unur_distr_free(distr);

  /* error: invalid parameters */
  distr = unur_distr_normal(fpar,2);
  par = unur_srou_new(distr);
  check_expected_setfailed( unur_srou_set_cdfatmode(par,-1.) );
  check_errorcode( UNUR_ERR_PAR_SET );
  check_expected_setfailed( unur_srou_set_pdfatmode(par,-1.) );
  check_errorcode( UNUR_ERR_PAR_SET );
  free(par); 
  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_srou_set() */

/*---------------------------------------------------------------------------*/

void test_srou_chg( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  double fpar[2] = {0.,1.};

  /* start test */
  printf("[chg ");  fflush(stdout);
  fprintf(TESTLOG,"\n[chg]\n");

  test_failed = 0;

  /* check error handling */

  /* error: invalid generator object */
  distr = unur_distr_normal(fpar,2);
  par = unur_arou_new(distr);
  unur_set_debug(par,0);
  gen = unur_init( par );
  check_expected_setfailed( unur_srou_chg_pdfparams(gen,fpar,2) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  check_expected_setfailed( unur_srou_chg_domain(gen,0.,1.) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  check_expected_setfailed( unur_srou_chg_mode(gen,0.) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  check_expected_setfailed( unur_srou_chg_cdfatmode(gen,1.) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  check_expected_setfailed( unur_srou_chg_pdfatmode(gen,1.) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  check_expected_setfailed( unur_srou_chg_pdfarea(gen,1.) );
  check_errorcode( UNUR_ERR_GEN_INVALID );
  unur_free(gen); 
  unur_distr_free(distr);

  /* error: invalid parameters */
  distr = unur_distr_normal(fpar,2);
  par = unur_srou_new(distr);
  unur_set_debug(par,0);
  gen = unur_init( par );
  check_expected_setfailed( unur_srou_chg_pdfparams(gen,NULL,1) );
  check_errorcode( UNUR_ERR_NULL );
  check_expected_setfailed( unur_srou_chg_pdfparams(gen,fpar,UNUR_DISTR_MAXPARAMS+10) );
  check_errorcode( UNUR_ERR_DISTR_NPARAMS );
  check_expected_setfailed( unur_srou_chg_cdfatmode(gen,-1.) );
  check_errorcode( UNUR_ERR_PAR_SET );
  check_expected_setfailed( unur_srou_chg_pdfatmode(gen,-1.) );
  check_errorcode( UNUR_ERR_PAR_SET );
  check_expected_setfailed( unur_srou_chg_domain(gen,1.,-1.) );
  check_errorcode( UNUR_ERR_DISTR_SET );
  check_expected_setfailed( unur_srou_chg_pdfarea(gen,-1.) );
  check_errorcode( UNUR_ERR_DISTR_SET );
  unur_free(gen); 
  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_srou_chg() */

/*---------------------------------------------------------------------------*/

void test_srou_init( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;

  /* start test */
  printf("[init ");  fflush(stdout);
  fprintf(TESTLOG,"\n[init]\n");

  test_failed = 0;

  /* check error handling */
  distr = unur_distr_cont_new();
  unur_distr_cont_set_pdf(distr,pdf);
  unur_distr_cont_set_mode(distr,0.);
  unur_distr_cont_set_pdfarea(distr,1.);
  par = unur_srou_new(distr);
  check_expected_NULL( unur_init(par) );
  check_errorcode( UNUR_ERR_GEN_DATA );
  unur_distr_free(distr);

  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_srou_init() */

/*---------------------------------------------------------------------------*/

void test_srou_reinit( void )
{
  /* nothing to do */

  printf("[reinit ");  fflush(stdout);
  fprintf(TESTLOG,"\n[reinit]\n");

  test_failed = 0;

  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_srou_reinit() */

/*---------------------------------------------------------------------------*/

void test_srou_sample( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  int i;
  double fpar[2] = {0.,1.};

#define N_SAMPLE 500
  double sa[N_SAMPLE], sb[N_SAMPLE];

  /* start test */
  printf("[sample ");  fflush(stdout);
  fprintf(TESTLOG,"\n[sample]\n");

  test_failed = 0;

  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* we use the normal distributions for the next tests */
  distr = unur_distr_normal(fpar,2);

  /* the following generator should generate the same sequence */

  /* basic algorithm */
  prng_reset(urng);
  par = unur_srou_new(distr);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sa[i] = unur_sample_cont(gen);
  unur_free(gen); 

  /* basic algorithm - verifying mode */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = unur_sample_cont(gen);
  unur_free(gen); 
  compare_sequences(sa,sb,N_SAMPLE);

  /* the following generator should generate the same sequence */

  /* use cdf at mode */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,0.5);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sa[i] = unur_sample_cont(gen);
  unur_free(gen); 

  /* use cdf at mode and squeeze */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,0.5);
  unur_srou_set_usesqueeze(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = unur_sample_cont(gen);
  unur_free(gen); 
  compare_sequences(sa,sb,N_SAMPLE);

  /* use cdf at mode - verifying mode */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,0.5);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = unur_sample_cont(gen);
  unur_free(gen); 
  compare_sequences(sa,sb,N_SAMPLE);

  /* use cdf at mode and squeeze - verifying mode */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,0.5);
  unur_srou_set_usesqueeze(par,1);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = unur_sample_cont(gen);
  unur_free(gen); 
  compare_sequences(sa,sb,N_SAMPLE);

  /* the following generator should generate the same sequence */

  /* use mirror principle */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_usemirror(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sa[i] = unur_sample_cont(gen);
  unur_free(gen); 

  /* use mirror principle - verifying mode */
  prng_reset(urng);
  par = unur_srou_new(distr);
  unur_srou_set_usemirror(par,1);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = unur_sample_cont(gen);
  unur_free(gen); 
  compare_sequences(sa,sb,N_SAMPLE);


  /* the following distribution violates the condition of the method:
     pdf at mode is too small --> hat < pdf near mode */
  unur_errno = 0;
  par = unur_srou_new(distr);
  unur_srou_set_pdfatmode(par,0.1);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<20; i++)
    unur_sample_cont(gen);
  unur_free(gen); 
  check_errorcode( UNUR_ERR_GEN_CONDITION );
  
  unur_errno = 0;
  par = unur_srou_new(distr);
  unur_srou_set_pdfatmode(par,0.1);
  unur_srou_set_cdfatmode(par,0.5);
  unur_srou_set_usesqueeze(par,1);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  for (i=0; i<20; i++)
    unur_sample_cont(gen);
  unur_free(gen); 
  check_errorcode( UNUR_ERR_GEN_CONDITION );
  
  /* the following distribution violates the condition of the method:
     pdf at mode is too large --> squeeze > pdf near mode */
  unur_errno = 0;
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,0.5);
  unur_srou_set_pdfatmode(par,10.);
  unur_srou_set_usesqueeze(par,1);
  unur_srou_set_verify(par,1);
  gen = unur_init( par );
  unur_srou_chg_pdfarea(gen,10.);
  unur_reinit(gen);
  for (i=0; i<20; i++)
    unur_sample_cont(gen);
  unur_free(gen); 
  check_errorcode( UNUR_ERR_GEN_CONDITION );

  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

#undef N_SAMPLE
} /* end of test_srou_sample() */

/*---------------------------------------------------------------------------*/

void test_srou_validate(void)
{
  int n;
  UNUR_DISTR *distr;
  const char *distr_name;
  const char *last_distr_name = "";
  
  /* start test */
  printf("[validate: ");  fflush(stdout);
  fprintf(TESTLOG,"\n[validate]\n");

  test_failed = 0;

  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* run chi^2 tests on test distributions */
  for (n=0; n<n_distr; n++) {

    /* test type of distribution */
    if( (list_of_distr[n].type != T_TYPE_TDR) ||
	(list_of_distr[n].c_max < -0.5) )
      /* we cannot use this method for this distribution */
      continue;

    /* get pointer to distribution */
    distr = list_of_distr[n].distr;

    /* get name of distribution */
    distr_name = unur_distr_get_name(distr);
    if (strcmp(distr_name,last_distr_name) ) {
      /* different distributions */
      last_distr_name = distr_name;
      printf(" %s",distr_name); fflush(stdout);
    }

    /* run tests */
    test_srou_validate_chi2(distr);
  }

  /* run level 2 test on collected p-values */
  run_level2(__LINE__, list_pvals, n_pvals);

  /* test finished */
  test_ok &= (test_failed > FAILED_LIMIT) ? 0 : 1;
  (test_failed>FAILED_LIMIT) ? printf(" --> failed] ") : printf(" --> ok] ");

} /* end of test_srou_validate() */

/*...........................................................................*/

void test_srou_validate_chi2( UNUR_DISTR *distr )
{
  UNUR_PAR *par;
  double cdfatmode;

  /* basic algorithm */
  par = unur_srou_new(distr);
  run_chi2(par,distr,__LINE__); 

  /* use mirror principle */
  par = unur_srou_new(distr);
  unur_srou_set_usemirror(par,1);
  run_chi2(par,distr,__LINE__); 

  /* get cdf at mode */
  cdfatmode = unur_distr_cont_cdf( unur_distr_cont_get_mode(distr), distr );
  if (cdfatmode < 0.)
    /* cdf not available */
    return;

  /* use cdf at mode and squeeze */
  par = unur_srou_new(distr);
  unur_srou_set_cdfatmode(par,cdfatmode);
  unur_srou_set_usesqueeze(par,1);
  run_chi2(par,distr,__LINE__); 

} /* end of test_srou_validate_chi2() */  

/*...........................................................................*/

void run_chi2( UNUR_PAR *par, UNUR_DISTR *distr, int line )
{
  UNUR_GEN *gen;
  double pval;
  int i;

  gen = unur_init(par);
  if (gen==NULL) {
    /* this must not happen */
    ++test_failed;
    fprintf(TESTLOG,"line %4d: pval =     Initialization failed\t\t",line);
    print_distr_name( distr,"");
    fprintf(TESTLOG,"\n");
    printf("0");
    return;
  }

  /* run chi^2 test */
  for (i=1; i<=2; i++) {
    /* we run the test twice when it fails the first time */
    pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 0, 20, 0);
    do_check_pval(line,gen,pval,i);

    /* store p-value */
    if (n_pvals >= size_pvals) {
      /* no space left in list: enlarge with 100 entries */
      size_pvals += 100;
      list_pvals = realloc( list_pvals, size_pvals * sizeof(double) );
    }
    /* append to list */
    list_pvals[n_pvals++] = pval;

    if (pval >= PVAL_LIMIT || pval < 0.) 
      /* test succeeded or not performed */
      break; 
  }

  unur_free(gen);

} /* end of run_chi2() */

/*---------------------------------------------------------------------------*/

/* pdf that does not work */
double pdf( double x, UNUR_DISTR *distr )
{ 
  return ((x==0.) ? 0. : exp(-x*x/2.));
} /* end of pdf */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_SROU */
/*---------------------------------------------------------------------------*/
