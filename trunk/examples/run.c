/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Examples                                                                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <malloc.h>

#include <unuran.h>
#include <unuran_tests.h>

#define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  int i;
  double vec[3];

  double fpm[10];

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  distr = unur_distr_normal(NULL,0);
  unur_distr_cont_set_domain(distr,1,20);

  par = unur_tdr_new(distr);
  unur_tdr_set_variant_ps(par);

  gen = unur_init(par);

  // unur_tdr_chg_truncated(gen,3,3.1);

  unur_test_chi2( gen, 100, 0, 20, 1 );


  // unur_run_tests( par, RUN_TESTS) ;

  unur_distr_free(distr);

  exit (0);
}

/*---------------------------------------------------------------------------*/







