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

int main()
{

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      
  double fpar[10];

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  fpar[0] = 1.;
  fpar[1] = 2.;
  distr = unur_distr_beta(fpar,2);
  unur_distr_cont_set_domain(distr,0.5,1.);
  unur_distr_cont_upd_mode(distr);
  unur_distr_cont_upd_pdfarea(distr);


  par = unur_srou_new(distr);
  gen = unur_init(par);

  unur_test_printsample( gen, 3, 10 );
  unur_test_count_urn( gen, 100000 );
  unur_test_chi2( gen, 100, 0, 20, 1 );

  unur_srou_chg_domain(gen,0.9,0.91);
  unur_srou_upd_mode(gen);
  unur_srou_upd_pdfarea(gen);
  unur_srou_reinit(gen);	

  unur_test_printsample( gen, 3, 10 );
  unur_test_count_urn( gen, 100000 );
  unur_test_chi2( gen, 100, 0, 20, 1 );

  unur_free(gen);
  unur_distr_free(distr);

  exit (0);
}

/*---------------------------------------------------------------------------*/







