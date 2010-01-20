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
  double bs = 0.001;
  double bn = 0.001/3.;
  double s  = 12.;

  return ((bs-1.)*log(x)+log(1.-x)*(bn-1.)+s*(1.-x));
}

double mydlogpdf( double x, const UNUR_DISTR *distr )
{
  double bs = 0.001;
  double bn = 0.001/3.;
  double s  = 12.;

  return ((bs-1.)/x - (bn-1.)/(1.-x)-s);
}

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  unur_set_default_debug(~0U);

  distr = unur_distr_cont_new();
  unur_distr_cont_set_logpdf( distr, mylogpdf );
  unur_distr_cont_set_dlogpdf( distr, mydlogpdf );
  unur_distr_cont_set_domain( distr, 0., 0.5);
  unur_distr_cont_set_mode( distr, 0.);

  par = unur_itdr_new(distr);
  unur_itdr_set_cp(par,-0.999);

  gen = unur_init(par);

  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  for (i=1; i<100; i++)
    printf("%g\n",unur_sample_cont(gen));

  return 0;
}

/*---------------------------------------------------------------------------*/

