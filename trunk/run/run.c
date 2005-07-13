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
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  double fpm[10];

  unur_set_default_debug(~0U);


  fpm[0] = 0.2;
  distr = unur_distr_gamma(fpm,1);
  par = unur_ninv_new(distr);
  unur_ninv_set_usenewton(par);
  unur_ninv_set_x_resolution(par,1.e-14);
  unur_ninv_set_max_iter (par, 100);

/*   unur_run_tests( par, RUN_TESTS); */

  gen = unur_init(par);
  unur_test_chi2(gen,100,0,0,2,stdout);
  unur_distr_free(distr);
  
  return 0;
}

/*---------------------------------------------------------------------------*/

