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
#include <stdlib.h>

#include <unuran.h>
#include <unuran_tests.h>

#define RUN_TESTS       (~0x0u)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;
  double x;

  double fpm[]={0.01,2.};

  unur_set_default_debug(~0u);

  distr = unur_distr_beta(fpm,2);

  par = unur_cstd_new(distr);
  gen = unur_init(par);

  for (i=0;i<1000;i++)
    unur_sample_cont(gen);

  unur_test_chi2(gen,100,0,20,3,stdout);

  unur_distr_free(distr);
  unur_free(gen);

/*    unur_run_tests(par,UNUR_TEST_CHI2); */

  return 0;

}

/*---------------------------------------------------------------------------*/













