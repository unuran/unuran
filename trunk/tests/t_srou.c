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
int test_ok_local = TRUE;         /* running test ok (boolean)               */
FILE *TESTLOG = NULL;             /* test log file                           */

static FILE *UNURANLOG = NULL;    /* unuran log file                         */

/*---------------------------------------------------------------------------*/

void test_srou_new( void );
void test_srou_set( void );
void test_srou_chg( void );
void test_srou_init( void );
void test_srou_reinit( void );
void test_srou_sample( void );

void test_srou_validate( void );
void test_srou_validate_chi2( UNUR_DISTR *distr );

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
  printf("\n");

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
  printf("[new ");
  fprintf(TESTLOG,"\n[new]\n");

  test_ok_local = TRUE;

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
  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");
} /* end of test_srou_new() */

/*---------------------------------------------------------------------------*/

void test_srou_set( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  double fpar[2] = {0.,1.};

  /* start test */
  printf("[set ");
  fprintf(TESTLOG,"\n[set]\n");

  test_ok_local = TRUE;

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
  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

} /* end of test_srou_set() */

/*---------------------------------------------------------------------------*/

void test_srou_chg( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  double fpar[2] = {0.,1.};

  /* start test */
  printf("[chg ");
  fprintf(TESTLOG,"\n[chg]\n");

  test_ok_local = TRUE;

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
  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

} /* end of test_srou_chg() */

/*---------------------------------------------------------------------------*/

void test_srou_init( void )
{
  /* nothing to do */

  printf("[init ");
  fprintf(TESTLOG,"\n[init]\n");

  test_ok_local = TRUE;

  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

} /* end of test_srou_init() */

/*---------------------------------------------------------------------------*/

void test_srou_reinit( void )
{
  /* nothing to do */

  printf("[reinit ");
  fprintf(TESTLOG,"\n[reinit]\n");

  test_ok_local = TRUE;

  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

} /* end of test_srou_reinit() */

/*---------------------------------------------------------------------------*/

void test_srou_sample( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  struct prng *urng;

  int i;
  double fpar[2] = {0.,1.};

#define N_SAMPLE 500
  double sa[N_SAMPLE], sb[N_SAMPLE];

  /* start test */
  printf("[sample ");
  fprintf(TESTLOG,"\n[sample]\n");

  test_ok_local = TRUE;

  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* we need a uniform random number generator */
  urng = prng_new("mt19937(5678)");
  unur_set_default_urng(urng);

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
  prng_free(urng);

  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

#undef N_SAMPLE
} /* end of test_srou_sample() */

/*---------------------------------------------------------------------------*/

void test_srou_validate(void)
{
  struct prng *urng;
  int n;

  /* start test */
  printf("[validate ");
  fprintf(TESTLOG,"\n[validate]\n");

  test_ok_local = TRUE;

  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* we need a uniform random number generator */
  urng = prng_new("mt19937(5678)");
  unur_set_default_urng(urng);

  /* run chi^2 tests on test distributions */
  for (n=0; n<n_distr; n++)
    test_srou_validate_chi2(list_of_distr[n]);
    
  /* test finished */
  prng_free(urng);

  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");

} /* end of test_srou_validate() */

void test_srou_validate_chi2( UNUR_DISTR *distr )
{
  UNUR_PAR *par;
  UNUR_GEN *gen;
  double pval;

#define TEST_CHI2 \
 gen=unur_init(par);\
 pval=unur_test_chi2(gen,CHI_TEST_INTERVALS,0,20,0);\
 check_pval(gen,pval);\
 unur_free(gen);\
 unur_distr_free(distr)

  par = unur_srou_new(distr);
  TEST_CHI2; 

#undef TEST_CHI2
} /* end of test_srou_validate_chi2() */  

/*---------------------------------------------------------------------------*/

double pdf( double x, UNUR_DISTR *distr )
{ 
  return exp(-x*x/2.);
} /* end of pdf */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_SROU */
/*---------------------------------------------------------------------------*/
