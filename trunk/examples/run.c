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

  UNUR_DISTR *distr, *os;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      
  double fpar[10];

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

#if 0
  fpar[0] = 3.;
  fpar[1] = 6.;
  distr = unur_distr_beta(fpar,2);
  unur_distr_cont_set_domain(distr,-1.,3.);
  par = unur_tdr_new(distr);
  unur_tdr_set_variant_ps(par);
  unur_run_tests( par, RUN_TESTS);
#endif

#if 1
  fpar[0] = 0.;
  fpar[1] = 100.;
  distr = unur_distr_normal(fpar,2);
  os = unur_distr_corder_new( distr, 10, 8);

  par = unur_tdr_new(os);
  // unur_tdr_set_variant_ps(par);
  unur_tdr_set_max_intervals(par,100);
  unur_tdr_set_cpoints(par,47,NULL);
  //  unur_tdr_set_c(par,0.);
  unur_run_tests( par, RUN_TESTS);
#endif

  exit (0);
}

/*---------------------------------------------------------------------------*/







