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

  double fpm[] = { 2., 3. };

  unur_set_default_debug(~0u);
  unur_set_stream(stdout);

#if 1
  fpm[0] = 0.;
  fpm[1] = 1e-5;
  //distr = unur_distr_normal(NULL,0);
  distr = unur_distr_normal(fpm,2);
  par = unur_tdr_new( distr );
  unur_tdr_set_cpoints(par,10,NULL);
  unur_tdr_set_max_intervals(par,100);
  //unur_tdr_set_usedars(par,FALSE);
  //unur_tdr_set_darsfactor(par,0.);
  //unur_tdr_set_variant_gw(par);
  //unur_tdr_set_c(par,0.);
  gen = unur_init(par);
  unur_distr_free(distr);
#else
  distr = unur_distr_beta(fpm,2);
  par = unur_tdr_new( distr );
  unur_tdr_set_cpoints(par,10,NULL);
  //  unur_tdr_set_darsfactor(par,0.);
  unur_tdr_set_variant_gw(par);
  unur_tdr_set_c(par,0.);
  //  unur_tdr_set_usedars(par,1);
  gen = unur_init(par);
  unur_distr_free(distr);
#endif

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new( distr ); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,1u); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,3u); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,UNUR_STDGEN_DEFAULT); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_exponential(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_powerexponential(fpm,1); */
/*    par = unur_tdr_new( distr ); */
/*    unur_tdr_set_variant_gw(par); */
/*    gen = unur_init( par ); */
/*    unur_distr_free(distr); */
/*    unur_free(gen); */

  exit (0);
}

/*---------------------------------------------------------------------------*/














