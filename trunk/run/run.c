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
  UNUR_URNG *qrng, *srng, *pg;

  double fpm[10];

  unur_set_default_debug(~0U);

  qrng = unur_urng_gslqrng_new(gsl_qrng_sobol,3);
  srng = unur_get_default_urng();
  pg = unur_urng_randomshift_new(qrng,srng,3);

  distr = unur_distr_normal(NULL,0);
  par = unur_hinv_new(distr);

  unur_set_urng(par,pg);

  
/*   gen = unur_init(par); */

  unur_run_tests( par, RUN_TESTS);

  unur_distr_free(distr);
  unur_urng_free(pg);
  unur_urng_free(qrng);
  unur_urng_free(srng);

  return 0;
}

/*---------------------------------------------------------------------------*/

