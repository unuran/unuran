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
#define dim (3)
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  double timing_setup, timing_sample;

  unur_set_default_debug(~0u);

  par = unur_unif_new(NULL);
  gen = unur_test_timing( par, 6, &timing_setup, &timing_sample, TRUE, stdout);

  unur_free(gen);
  
    return 0;
#undef dim
}

/*---------------------------------------------------------------------------*/

