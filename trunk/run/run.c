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

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  double unur_urng_mstd(void);

  double fpm[] = {1,5};

  unur_set_default_debug(~0u);
/*    unur_set_stream(stdout); */

  distr = unur_distr_normal(NULL,0);

  par = unur_tdr_new(distr);

/*    unur_set_urng(par,unur_urng_mstd); */

/*    gen = unur_init(par); */
/*    unur_srou_chg_domain(gen,0.9,0.95); */
/*    unur_srou_upd_pdfarea(gen); */
/*    unur_srou_upd_mode(gen); */
/*    unur_srou_reinit(gen);   */

/*    unur_test_chi2(gen, 10, 100, 5, 1, stdout); */

  
  unur_run_tests(par,RUN_TESTS);

  return 0;

}

/*---------------------------------------------------------------------------*/













