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

  double fpar[] = {1.1,3.2};

  unur_set_default_debug(~0U);

  distr = unur_distr_F(fpar,2);
  unur_distr_cont_set_center(distr,1);
  par = unur_pinv_new(distr);
  unur_run_tests(par,UNUR_TEST_N_PDF|UNUR_TEST_TIME,stdout);

  /* gen = unur_str2gen("normal & method=ninv; u_resolution=1e-10"); */
  /* for (i=0; i<1e1; i++)  */
  /*   unur_sample_cont(gen); */

  return 0;
}

/*---------------------------------------------------------------------------*/

