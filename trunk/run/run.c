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
  struct unur_slist *mlist;

  double time,time_0;

  double unur_urng_mstd(void);

  double fpm[] = {1,5};

  unur_set_default_debug(0u);
/*    unur_set_stream(stdout); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ia(par); */

  distr = unur_str2distr("normal");
  par = _unur_str2par(distr,"method=tdr; variant_ia",&mlist);

/*    par = unur_unif_new(NULL); */

/*    unur_set_urng(par,unur_urng_mstd); */

/*    gen = unur_init(par); */
/*    unur_srou_chg_domain(gen,0.9,0.95); */
/*    unur_srou_upd_pdfarea(gen); */
/*    unur_srou_upd_mode(gen); */
/*    unur_srou_reinit(gen);   */

/*    unur_test_chi2(gen, 10, 100, 5, 1, stdout); */

  time_0 = unur_test_timing_total(par, 1, 0.1 );
  printf("time_0 = %g\n\n",time_0);
  
  time = unur_test_timing_total(par, 10, 0.1 );
  printf("time_1 = %g (%g)\n\n",time,(time-time_0)/9);

  time = unur_test_timing_total(par, 100, 0.1 );
  printf("time_2 = %g (%g)\n\n",time,(time-time_0)/99);

  time = unur_test_timing_total(par, 1000, 0.1 );
  printf("time_3 = %g (%g)\n\n",time,(time-time_0)/999);
  
  time = unur_test_timing_total(par, 10000, 0.1 );
  printf("time_4 = %g (%g)\n\n",time,(time-time_0)/9999);

  time = unur_test_timing_total(par, 100000, 0.1 );
  printf("time_5 = %g (%g)\n\n",time,(time-time_0)/99999);

  time = unur_test_timing_total(par, 1000000, 0.1 );
  printf("time_6 = %g (%g)\n\n",time,(time-time_0)/999999);

  time = unur_test_timing_total(par, 10000000, 0.1 );
  printf("time_7 = %g (%g)\n\n",time,(time-time_0)/9999999);

  time = unur_test_timing_total(par, 100000000, 0.1 );
  printf("time_8 = %g (%g)\n\n",time,(time-time_0)/99999999);

  unur_run_tests(par,RUN_TESTS);

  _unur_slist_free(mlist);

  return 0;

}

/*---------------------------------------------------------------------------*/













