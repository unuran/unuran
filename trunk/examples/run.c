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

int main()
{

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      
  double fpar[10];

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  fpar[0] = 1.;
  fpar[1] = 2.;
  distr = unur_distr_beta(fpar,2);

  //  distr = unur_distr_normal(NULL,0);

  par = unur_cstd_new(distr);
  unur_run_tests(par, RUN_TESTS);

  par = unur_cstd_new(distr);
  unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION);
  unur_run_tests(par, RUN_TESTS);


  par = unur_cstd_new(distr);
  unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION);
  gen = unur_init(par);
  unur_cstd_chg_truncated(gen,0.5,0.55);
  unur_test_chi2( gen, 100, 0, 20, 1 );
  unur_free(gen);

  unur_distr_free(distr);

  exit (0);
}

/*---------------------------------------------------------------------------*/







