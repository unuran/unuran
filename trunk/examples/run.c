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

  fpar[0] = 2.;
  distr = unur_distr_gamma(fpar,1);
  par = unur_arou_new(distr);
  unur_arou_set_cpoints( par, 8, NULL );
  unur_arou_set_max_sqhratio(par,0.5);

  unur_run_tests(par,RUN_TESTS);

  exit (0);
}

/*---------------------------------------------------------------------------*/







