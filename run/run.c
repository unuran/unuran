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

#include <experimental/itdr.h>

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
  double fpar[3] = { 0.5, 1., 0. };
 
  unur_set_default_debug(~0U);

  /* standard normal */
  distr = unur_distr_gamma(fpar,3);
/*   distr = unur_distr_beta(fpar,2); */
/*   unur_distr_cont_set_domain(distr,fpar[2],10); */
  par = unur_itdr_new(distr);
/*   unur_itdr_set_bx(par,fpar[0]+fpar[2]); */
/*   unur_itdr_set_cp(par,-0.5); */
/*   unur_itdr_set_ct(par,-0.5); */

/*   gen = unur_init(par); */

  unur_run_tests(par,~0u);
/*   unur_run_tests(par,UNUR_TEST_SAMPLE); */
  
/*   unur_free(gen); */
  unur_distr_free(distr);

  return 0;
}

/*---------------------------------------------------------------------------*/

