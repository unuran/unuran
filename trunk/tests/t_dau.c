/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests for DAU                                                            *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/
#ifdef T_DAU
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static struct prng *urng;         /* uniform random number generator         */

#define FAILED_LIMIT 1            /* maximal number of single failed 
				     statistical tests (if there are more 
				     the generator has failed the test.)     */

static double *list_pvals = NULL; /* list of collected p-values              */
static int size_pvals = 0;        /* size of list                            */
static int n_pvals = 0;           /* number of collected p-values            */

/*---------------------------------------------------------------------------*/

void test_dau_new( void );
void test_dau_set( void );
void test_dau_chg( void );
void test_dau_init( void );
void test_dau_reinit( void );
void test_dau_sample( void );

void test_dau_validate( void );
void test_dau_validate_chi2( UNUR_DISTR *distr );
void run_chi2( UNUR_PAR *par, UNUR_DISTR *distr, int line );

UNUR_DISTR *get_distr_with_pv( void );
UNUR_DISTR *get_distr_with_invalid_pv( void );

/*---------------------------------------------------------------------------*/

int main()
{ 
  /* open log files for testing */
  open_log_files(__FILE__);

  /* make a list of distributions */
  make_list_of_distributions();

  /* we use Mersenne Twister as uniform random number generator */
  urng = prng_new("mt19937(5678)");
  unur_set_default_urng(urng);

  /* set default debugging flag */
  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* start test */
  printf("%s: ",__FILE__);

  /* run tests */
  test_dau_new();
  test_dau_set();
  test_dau_chg();
  test_dau_init();
  test_dau_reinit();
  test_dau_sample();
  test_dau_validate();

  /* test finished */
  printf("\n");  fflush(stdout);

  /* close log files and exit */
  close_log_files();
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

void test_dau_new( void )
{
  UNUR_DISTR *distr;
  double fpar[1] = {0.5};

  /* start test */
  printf("[new "); fflush(stdout);
  fprintf(TESTLOG,"\n[new]\n");

  test_failed = 0;

  /* check error handling */

  /* error: NULL pointer */
  distr = NULL;   
  check_expected_NULL( unur_dau_new(distr) );
  check_errorcode( UNUR_ERR_NULL );
     
  /* error: invalid distribution type */
  distr = unur_distr_cont_new();

  check_expected_NULL( unur_dau_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_INVALID );

  unur_distr_free(distr);

  /* error: data missing in distribution object */
  distr = unur_distr_geometric(fpar,1);   /* no probability vector */
  check_expected_NULL( unur_dau_new(distr) );
  check_errorcode( UNUR_ERR_DISTR_REQUIRED );

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_new() */

/*---------------------------------------------------------------------------*/

void test_dau_set( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;

  /* start test */
  printf("[set "); fflush(stdout);
  fprintf(TESTLOG,"\n[set]\n");

  test_failed = 0;

  /* check error handling */

  /* error: NULL pointer */
  par = NULL;   

  check_expected_setfailed( unur_dau_set_urnfactor(par,1.) );
  check_errorcode( UNUR_ERR_NULL );

  /* error: invalid parameter object */
  distr = get_distr_with_pv();
  par = unur_dgt_new(distr);

  check_expected_setfailed( unur_dau_set_urnfactor(par,1.) );
  check_errorcode( UNUR_ERR_PAR_INVALID );

  free(par); 
  unur_distr_free(distr);

  /* error: invalid parameters */
  distr = get_distr_with_pv();
  par = unur_dau_new(distr);

  check_expected_setfailed( unur_dau_set_urnfactor(par,0.8) );
  check_errorcode( UNUR_ERR_PAR_SET );

  free(par); 
  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_set() */

/*---------------------------------------------------------------------------*/

void test_dau_chg( void )
{
  /* start test */
  printf("[chg ");  fflush(stdout);
  fprintf(TESTLOG,"\n[chg]\n");

  test_failed = 0;

  /* nothing to do */

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_chg() */

/*---------------------------------------------------------------------------*/

void test_dau_init( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;

  /* start test */
  printf("[init ");  fflush(stdout);
  fprintf(TESTLOG,"\n[init]\n");

  test_failed = 0;

  /* check error handling */
  distr = get_distr_with_invalid_pv();
  par = unur_dau_new(distr);

  check_expected_NULL( unur_init(par) );
  check_errorcode( UNUR_ERR_GEN_DATA );

  unur_distr_free(distr);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_init() */

/*---------------------------------------------------------------------------*/

void test_dau_reinit( void )
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  /* start test */
  printf("[reinit ");  fflush(stdout);
  fprintf(TESTLOG,"\n[reinit]\n");

  test_failed = 0;

  /* check whether succesfully executed */
  distr = get_distr_with_pv();
  par = unur_dau_new(distr);
  gen = unur_init(par);

  fprintf(TESTLOG,"line %4d: no reinit ...\t\t",__LINE__);  
  if (unur_reinit(gen)) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_reinit() */

/*---------------------------------------------------------------------------*/

void test_dau_sample( void )
{

  /* start test */
  printf("[sample ");  fflush(stdout);
  fprintf(TESTLOG,"\n[sample]\n");

  test_failed = 0;

  /* nothing to do */

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_dau_sample() */

/*---------------------------------------------------------------------------*/

void test_dau_validate(void)
{
  int n;
  UNUR_DISTR *distr;
  const char *distr_name;
  const char *last_distr_name = "";
  
  /* start test */
  printf("[validate: ");  fflush(stdout);
  fprintf(TESTLOG,"\n[validate]\n");

  test_failed = 0;

  /* reset default debugging flag */
  unur_set_default_debug(1);

  /* run chi^2 tests on test distributions */
  for (n=0; n<n_distr; n++) {

    /* test type of distribution */
    if( !(list_of_distr[n].type & T_TYPE_PV) )
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
    test_dau_validate_chi2(distr);
  }

  /* run level 2 test on collected p-values */
  run_level2(__LINE__, list_pvals, n_pvals);

  /* test finished */
  test_ok &= (test_failed > FAILED_LIMIT) ? 0 : 1;
  (test_failed>FAILED_LIMIT) ? printf(" --> failed] ") : printf(" --> ok] ");

} /* end of test_dau_validate() */

/*...........................................................................*/

void test_dau_validate_chi2( UNUR_DISTR *distr )
{
  UNUR_PAR *par;

  /* default algorithm */
  par = unur_dau_new(distr);
  run_chi2(par,distr,__LINE__); 

  /* larger urn table size */
  par = unur_dau_new(distr);
  unur_dau_set_urnfactor(par,10);
  run_chi2(par,distr,__LINE__); 

} /* end of test_dau_validate_chi2() */  

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
    pval = unur_test_chi2( gen, CHI_TEST_INTERVALS, 100000, 20, 0);
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

UNUR_DISTR *get_distr_with_pv( void )
{
#define PV_SIZE 200
  double prob[PV_SIZE];
  UNUR_DISTR *distr;
  int i;

  for (i=0; i<PV_SIZE; i++)
    prob[i] = prng_get_next(urng);

  distr = unur_distr_discr_new();
  unur_distr_discr_set_prob(distr,prob,PV_SIZE);
  abort_if_NULL(distr);

  return distr;

#undef PV_SIZE
} /* end of get_distr_with_pv() */

/*---------------------------------------------------------------------------*/

UNUR_DISTR *get_distr_with_invalid_pv( void )
{
#define PV_SIZE 10
  double prob[PV_SIZE];
  UNUR_DISTR *distr;
  int i;

  prob[0] = -1.;    /* invalid ! */

  for (i=1; i<PV_SIZE; i++)
    prob[i] = prng_get_next(urng);

  distr = unur_distr_discr_new();
  unur_distr_discr_set_prob(distr,prob,PV_SIZE);
  abort_if_NULL(distr);

  return distr;
  
#undef PV_SIZE
} /* end of get_distr_with_invalid_pv() */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_DAU */
/*---------------------------------------------------------------------------*/
