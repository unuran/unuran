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
  double stp[]={-2,-1,0,1,2,3};

  unur_set_default_debug(~0u);

  distr = unur_distr_normal(NULL,0);
/*    unur_distr_cont_set_domain(distr,0,1000); */

  par = unur_hinv_new(distr);
  unur_hinv_set_u_resolution(par,1.e-8);
  unur_hinv_set_order(par,1);
  unur_hinv_set_cpoints(par,stp,6);
  unur_hinv_set_max_intervals(par,1000);

  gen = unur_init(par);

  for (i=0;i<10;i++)
    printf("%g\n",unur_sample_cont(gen));

/*    unur_hinv_chg_truncated(gen,1,1.001); */

/*    for (i=0;i<10;i++) */
/*      printf("%g\n",unur_sample_cont(gen)); */

  printf("N = %d\n",unur_hinv_get_n_intervals(gen));

  unur_distr_free(distr);
  unur_free(gen);

/*    unur_run_tests(par,RUN_TESTS); */

  return 0;

}

/*---------------------------------------------------------------------------*/













