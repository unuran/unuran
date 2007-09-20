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
#include <testdistributions.h>
#include <mvtdr.h>

#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

#define unur_distr_multinormal  unur_distr_multinormal_w_marginals


/*---------------------------------------------------------------------------*/

/* double pdf( double x, const UNUR_DISTR *distr ) */
/* { */
/* /\*   return (1.e32 * pow(1-x*x,-50) * exp((78.043 * x - 107.415)/(1 - x*x))); *\/ */
/*   return exp(log(1.e32)+  50*log(1-x*x) + ((78.043 * x - 107.415)/(1 - x*x))); */
/* } */
 

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  double ll[] = { -1e-4, -1e-4};
  double ur[] = { 1, 1};

/*   double ll[] = {-UNUR_INFINITY,-UNUR_INFINITY}; */
/*   double ur[] = { UNUR_INFINITY, UNUR_INFINITY}; */
  
  double mean[] = {0.,0.};

  unur_set_default_debug(~0U);

/*   distr = unur_distr_multinormal(2,mean,NULL); */
  distr = unur_distr_multinormal(2,NULL,NULL);
  unur_distr_cvec_set_domain_rect(distr,ll,ur);
  unur_distr_cvec_set_mode(distr,mean);

  par = unur_mvtdr_new(distr);
  unur_mvtdr_set_stepsmin(par, 0 );
  unur_run_tests(par,RUN_TESTS,stdout);

/*   gen = unur_init(par); */
/*   unur_test_printsample (gen, 100, 1, stdout); */
/*   unur_test_chi2( gen, 100, 0, 20, 1, stdout); */


  unur_distr_free(distr);

  return 0;
}

/*---------------------------------------------------------------------------*/

