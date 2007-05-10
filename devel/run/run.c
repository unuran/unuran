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
/*   double fpar[2] = {7,0.1}; */
/*   double fpar[4] = {3.,0.5, -1., 0.}; */

  int i;

  unur_set_default_debug(~0U);

  gen = unur_str2gen("distr=cont;cdf='(x<=3)*(0.0555555555555556 +( -0.111111111111111 )*x+( 0.0555555555555556 )*x*x)+(x> 3 )*( -0.587301587301587 +( 0.317460317460317 )*x+( -0.0158730158730159 )*x*x)'; domain=( 2 , 11 ) & method=hinv");


  unur_test_printsample(gen, 10, 10, stdout);

  /*   unur_run_tests(par,RUN_TESTS); */
  
  unur_free(gen);
  unur_distr_free(distr);

  return 0;
}

/*---------------------------------------------------------------------------*/

