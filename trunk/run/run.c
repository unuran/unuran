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

  unur_set_default_debug(~0u);

  distr =unur_distr_normal(NULL,0);
  par = unur_arou_new(distr);

  unur_run_tests( par, ~0U);
  
  return 0;
}

/*---------------------------------------------------------------------------*/

