/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests for UNIF                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/
#ifdef T_UNIF
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static struct prng *urng;         /* uniform random number generator         */

/*---------------------------------------------------------------------------*/

void test_unif_new( void );
void test_unif_set( void );
void test_unif_chg( void );
void test_unif_init( void );
void test_unif_reinit( void );
void test_unif_sample( void );
void test_unif_validate( void );

/*---------------------------------------------------------------------------*/

int main()
{ 
  /* open log files for testing */
  open_log_files(__FILE__);

  /* make a list of distributions */
  make_list_of_distributions();

  /* we use Mersenne Twister as uniform random number generator */
  urng = prng_new("mt19937(6202)");
  unur_set_default_urng(urng);

  /* set default debugging flag */
  unur_set_default_debug(UNUR_DEBUG_ALL);

  /* start test */
  printf("%s: ",__FILE__);

  /* run tests */
  test_unif_new();
  test_unif_set();
  test_unif_chg();
  test_unif_init();
  test_unif_sample();
  test_unif_reinit();
  test_unif_validate();

  /* test finished */
  printf("\n");  fflush(stdout);

  /* close log files and exit */
  close_log_files();
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

void test_unif_new( void )
{
  /* start test */
  printf("[new "); fflush(stdout);
  fprintf(TESTLOG,"\n[new]\n");

  test_failed = 0;

  /* nothing to do */

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_unif_new() */

/*---------------------------------------------------------------------------*/

void test_unif_set( void )
{
  /* start test */
  printf("[set "); fflush(stdout);
  fprintf(TESTLOG,"\n[set]\n");

  test_failed = 0;

  /* nothing to do */

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_unif_set() */

/*---------------------------------------------------------------------------*/

void test_unif_chg( void )
{
  /* start test */
  printf("[chg ");  fflush(stdout);
  fprintf(TESTLOG,"\n[chg]\n");

  test_failed = 0;

  /* nothing to do */

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_unif_chg() */

/*---------------------------------------------------------------------------*/

void test_unif_init( void )
{
  UNUR_PAR *par;

  /* start test */
  printf("[init ");  fflush(stdout);
  fprintf(TESTLOG,"\n[init]\n");

  test_failed = 0;

  /* check whether succesfully executed */
  par = unur_unif_new();

  fprintf(TESTLOG,"line %4d: init ...\t\t\t",__LINE__);  
  if (!unur_init(par)) {
    ++test_failed;
    fprintf(TESTLOG," Failed\n");
  }
  else
    fprintf(TESTLOG," ok\n");

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

} /* end of test_unif_init() */

/*---------------------------------------------------------------------------*/

void test_unif_reinit( void )
{
  UNUR_PAR *par;
  UNUR_GEN *gen;

  /* start test */
  printf("[reinit ");  fflush(stdout);
  fprintf(TESTLOG,"\n[reinit]\n");

  test_failed = 0;

  /* check whether succesfully executed */
  par = unur_unif_new();
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

} /* end of test_unif_reinit() */

/*---------------------------------------------------------------------------*/

void test_unif_sample( void )
{
  UNUR_PAR *par;
  UNUR_GEN *gen;

  int i;

#define N_SAMPLE 500
  double sa[N_SAMPLE], sb[N_SAMPLE];

  /* start test */
  printf("[sample ");  fflush(stdout);
  fprintf(TESTLOG,"\n[sample]\n");

  test_failed = 0;

  /* the following generators should generate the same sequence */

  /* UNIF */
  prng_reset(urng);
  par = unur_unif_new();
  gen = unur_init( par ); abort_if_NULL(gen);

  for (i=0; i<N_SAMPLE; i++)
    sa[i] = unur_sample_cont(gen);
  unur_free(gen); 

  /* the urng */
  prng_reset(urng);
  for (i=0; i<N_SAMPLE; i++)
    sb[i] = prng_get_next(urng);

  compare_double_sequences(sa,sb,N_SAMPLE);

  /* test finished */
  test_ok &= (test_failed) ? 0 : 1;
  (test_failed) ? printf("--> failed] ") : printf("--> ok] ");

#undef N_SAMPLE
} /* end of test_unif_sample() */

/*---------------------------------------------------------------------------*/

void test_unif_validate(void)
{
  /* start test */
  printf("[validate: ");  fflush(stdout);
  fprintf(TESTLOG,"\n[validate]\n");

  test_failed = 0;

  /* nothing to do:
     we do not test the underlying uniform random number generator.
  */

  /* test finished */
  (test_failed) ? printf(" --> failed] ") : printf(" --> ok] ");

} /* end of test_unif_validate() */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_UNIF */
/*---------------------------------------------------------------------------*/
