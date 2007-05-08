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

#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

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
/*   double fpar[2] = {0.002, 5.}; */
/*   double fpar[4] = {3.,0.5, -1., 0.}; */

  int i;

  unur_set_default_debug(~0U);

  distr = unur_distr_cont_new();
  unur_distr_cont_set_cdfstr(distr, "(x<=2)*(0.5-x+0.5*x*x)+(x>2)*(-3.5+3*x-0.5*x*x)");
  unur_distr_cont_set_domain(distr,1.,3.);

  par = unur_hinv_new(distr);
  unur_hinv_set_order(par,3);
  unur_hinv_set_u_resolution(par,1.e-10);
  unur_set_debug(par,1u);

  gen = unur_init(par);

  unur_test_printsample(gen, 10, 10, stdout);

/*   { */
/*     double max_error, MAE; */
/*     unur_hinv_estimate_error(gen, 10000, &max_error, &MAE); */
/*     printf("max error = %g\n",max_error); */
/*     printf("MAE       = %g\n",MAE); */
/*   } */

  /*   unur_run_tests(par,RUN_TESTS); */
  
  unur_free(gen);
  unur_distr_free(distr);

  return 0;
}

/*---------------------------------------------------------------------------*/

