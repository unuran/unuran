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

/* #define RUN_TESTS       (~0x0u) */
#define RUN_TESTS       UNUR_TEST_SAMPLE

/*---------------------------------------------------------------------------*/

int main()
{
#define dim (3)
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  double umin[]={-1,-1,-1};
  double umax[]={1,1,1};

  /*
  int i;
  double mean[dim], covar[dim*dim];
  UNUR_DISTR *covar_distr = unur_distr_correlation(dim);
  UNUR_GEN *covar_gen = unur_init(unur_mcorr_new(covar_distr));
  UNUR_GEN *mean_gen = unur_str2gen("normal(5,1)");
  unur_sample_matr(covar_gen,covar);
  for (i=0; i<dim; i++) 
    mean[i] = unur_sample_cont(mean_gen);
  distr = unur_distr_multinormal(dim,mean,covar);
  unur_distr_free(covar_distr);
  unur_free(covar_gen);
  unur_free(mean_gen); 
  */

  unur_set_default_debug(~0u);

/*   distr = unur_distr_multinormal(dim,mean,covar); */
  distr = unur_distr_multinormal(dim,NULL,NULL);
  par = unur_vnrou_new(distr);
  unur_vnrou_set_rect_u(par,umin,umax);
  unur_vnrou_set_rect_v(par,1);
  gen = unur_init(par);

  unur_test_printsample( gen, 3, 3, stdout);
  unur_test_count_urn( gen, 1000, 1, stdout );

  unur_test_chi2( gen, 10, 1000, 10, 2, stdout );


/*   unur_sample_matr(gen,x); */

/*   unur_run_tests(par,UNUR_TEST_SAMPLE & UNUR_TEST_CHI2); */
/*   unur_run_tests(par,UNUR_TEST_SAMPLE & UNUR_TEST_CHI2); */
/*   unur_run_tests(par,RUN_TESTS); */

  unur_free(gen);
  unur_distr_free(distr);
  
    return 0;
#undef dim
}

/*---------------------------------------------------------------------------*/

