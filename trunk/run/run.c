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

  double fpm[] = {1,5};

  unur_set_default_debug(~0u);
  unur_set_stream(stdout);

  distr = unur_distr_beta(fpm,2);
/*    unur_distr_cont_set_domain(distr,0.1,1.); */
/*    unur_distr_cont_upd_mode(distr); */
/*    unur_distr_cont_upd_pdfarea(distr); */

  par = unur_srou_new(distr);
  unur_srou_set_r(par,10.);
  unur_srou_set_verify(par,1);



/*    par = unur_srou_new(distr); */
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













