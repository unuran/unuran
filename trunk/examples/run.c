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

  fpar[0] = 1.;
  fpar[1] = 1.e-5;
  distr = unur_distr_normal(fpar,2);

  par = unur_ninv_new(distr);
  //  unur_ninv_set_usenewton(par);
  unur_ninv_set_table(par, 100);
  gen = unur_test_timing(par,5,fpar+8,fpar+9,1);
  unur_free(gen);

  exit (0);
}

/*---------------------------------------------------------------------------*/







