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

  double fpm[10];

  unur_set_default_debug(~0u);
  unur_set_stream(stdout);

  distr = unur_distr_normal(NULL,0);
  unur_distr_cont_set_domain(distr,0.1,1.);
  unur_distr_cont_upd_mode(distr);
  unur_distr_cont_upd_pdfarea(distr);

  {
    double cdfatmode = unur_distr_cont_eval_cdf( unur_distr_cont_get_mode(distr), distr );
    par = unur_srou_new(distr);
    unur_srou_set_cdfatmode(par,cdfatmode);
/*      unur_srou_set_usesqueeze(par,1);  */
  }

/*    par = unur_srou_new(distr); */
/*    gen = unur_init(par); */
  
  unur_run_tests(par,RUN_TESTS);

  return 0;

}

/*---------------------------------------------------------------------------*/













