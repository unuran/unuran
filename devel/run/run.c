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
  int i;

  unur_set_default_debug(~0U);

  distr = unur_distr_normal(NULL,0);
  par = unur_ninv_new(distr);
  unur_ninv_set_u_resolution(par,1.e-10);
  unur_ninv_set_x_resolution(par,-1);
  gen = unur_init(par);
  /* unur_run_tests(par,UNUR_TEST_N_PDF|UNUR_TEST_TIME,stdout); */

  /* gen = unur_str2gen("normal & method=ninv; u_resolution=1e-10"); */
  /* for (i=0; i<1e1; i++)  */
  /*   unur_sample_cont(gen); */

  return 0;
}

/*---------------------------------------------------------------------------*/

