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

#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

/*---------------------------------------------------------------------------*/

double pdf( double x, const UNUR_DISTR *distr )
{
/*   return (1.e32 * pow(1-x*x,-50) * exp((78.043 * x - 107.415)/(1 - x*x))); */
  return exp(log(1.e32)+  50*log(1-x*x) + ((78.043 * x - 107.415)/(1 - x*x)));
}


/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  unur_set_default_debug(~0U);

  distr = unur_distr_cont_new();
  unur_distr_cont_set_mode(distr,0.710413);
  unur_distr_cont_set_pdf(distr, pdf);

  par = unur_nrou_new(distr);

  gen = unur_init(par);

/*   unur_run_tests(par,~0u); */
/*   unur_run_tests(par,UNUR_TEST_SAMPLE); */
  
  unur_free(gen);
/*   unur_distr_free(normal); */

  return 0;
}

/*---------------------------------------------------------------------------*/

