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

double mypdf(double x,UNUR_DISTR *distr)
{
  // if (x<0.1 || x > 1.) return 0.;
  return 1./(1.+x*x);
}

int main()
{

  UNUR_DISTR *distr1, *distr2;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      
  double fpar[10];

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  fpar[0] = 0.;
  distr1 = unur_distr_cauchy(fpar,1);
  //  unur_distr_cont_set_domain(distr1,0.1,1.);
  unur_distr_cont_upd_mode(distr1);
  unur_distr_cont_upd_pdfarea(distr1);

  par = unur_srou_new(distr1);
  unur_srou_set_usemirror(par,1);
  unur_srou_set_verify(par,1);
  //  unur_run_tests(par, RUN_TESTS);


  distr2 = unur_distr_cont_new();
  unur_distr_cont_set_pdf(distr2,mypdf);
  unur_distr_cont_set_cdf(distr2, unur_distr_cont_get_cdf(distr1));
#if 0
  unur_distr_cont_set_mode(distr2,0.);
  unur_distr_cont_set_pdfarea(distr2, 3.141592655);
#else
  unur_distr_cont_set_domain(distr2,0.1,1.);
  unur_distr_cont_set_mode(distr2,0.1);
  unur_distr_cont_set_pdfarea(distr2, 0.6857295109);
#endif

  par = unur_srou_new(distr2);
  unur_srou_set_usemirror(par,1);
  unur_srou_set_verify(par,1);
  unur_run_tests(par, RUN_TESTS);


#if 0
  par = unur_cstd_new(distr);
  unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION);
  unur_run_tests(par, RUN_TESTS);


  par = unur_cstd_new(distr);
  unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION);
  gen = unur_init(par);
  unur_cstd_chg_truncated(gen,0.5,0.55);
  unur_test_chi2( gen, 100, 0, 20, 1 );
  unur_free(gen);
#endif

  exit (0);
}

/*---------------------------------------------------------------------------*/







