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

/* #define RUN_TESTS       (~0x0u) */
#define RUN_TESTS       UNUR_TEST_SAMPLE

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  double rk[] = { 1., 0.5, 0.1,  0.5, 1., 0.3,  0.1, 0.3, 1. };

  double erc[9];

  unur_set_default_debug(~0U);

  distr =unur_distr_copula(3,rk);
  par = unur_norta_new(distr);
  gen = unur_init(par);

/*   unur_run_tests( par, RUN_TESTS); */

  unur_test_cvec_rankcorr( erc, gen, 100000, 2, stdout );

  unur_free(gen);
  unur_distr_free(distr);
  
  return 0;
}

/*---------------------------------------------------------------------------*/

