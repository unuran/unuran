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

#define RUN_TESTS       (~0x0u)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  unur_set_default_debug(~0u);

  distr = unur_distr_normal(NULL,0);
/*    unur_distr_cont_set_domain(distr,0,1000); */

  par = unur_sinv_new(distr);
  unur_sinv_set_u_resolution(par,1.e-8);
  unur_sinv_set_order(par,1);

  gen = unur_init(par);

  for (i=0;i<10;i++)
    printf("%g\n",unur_sample_cont(gen));

/*    unur_sinv_chg_truncated(gen,1,1.001); */

/*    for (i=0;i<10;i++) */
/*      printf("%g\n",unur_sample_cont(gen)); */


  unur_distr_free(distr);
  unur_free(gen);

/*    unur_run_tests(par,RUN_TESTS); */

  return 0;

}

/*---------------------------------------------------------------------------*/













