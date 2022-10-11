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

double mylogpdf( double x, const UNUR_DISTR *distr )
{
  /* return fmin (1.-x*x/2, 0.); */
  /* return fmin (1.-x*x, 0.); */
  return fmin (-x*x, 0.);
}

double mypdf2( double x, const UNUR_DISTR *distr )
{
  return (x <= 1e-5)  ? exp(-x*x) : 3. * exp(-x*x);
}

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  unur_set_default_debug(~0U);

  /* Get empty distribution object for a continuous distribution */
  distr = unur_distr_cont_new();

  /* Fill the distribution object -- the provided information    */
  /* must fulfill the requirements of the method choosen below.  */
  /* unur_distr_cont_set_logpdf(distr,  mylogpdf);     /\* PDF             *\/ */
  unur_distr_cont_set_pdf(distr,  mypdf2);     /* PDF             */
  unur_distr_cont_set_domain(distr,-1.5,2.);

  
  par = unur_pinv_new(distr);
  /* unur_pinv_set_usecdf(par); */
  unur_pinv_set_order(par,5);

  gen = unur_init(par);

  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  for (i=1; i<10; i++)
    printf("%g\n",unur_sample_cont(gen));

  return 0;
}

/*---------------------------------------------------------------------------*/

