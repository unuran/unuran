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

/* pdf with piecewise linear function as transformed density with T = -1/sqrt */
double pdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
  double y;
  x -= 5.;
  y = 1./(fabs(x)+1.);
  return y*y;
}
double dpdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
  double y;
  x -= 5.;
  y = 1./(fabs(x)+1.);
  y = 2.*y*y*y;
  return ((x<0.) ? y : - y);
}
double cdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
  x -= 5.;
  if (x<=0.)
    return 0.5/(1.-x);
  else
    return (1.-0.5/(1.+x));
}

/*---------------------------------------------------------------------------*/

int main()
{
  int i;
  double vec[3];

  double fpm[10];

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  distr = unur_distr_cont_new();
  unur_distr_cont_set_pdf(distr,pdf_sqrtlin);
  unur_distr_cont_set_dpdf(distr,dpdf_sqrtlin);
  unur_distr_cont_set_cdf(distr,cdf_sqrtlin);
  unur_distr_set_name(distr,"sqrtlin");

  unur_errno = 0;
  par = unur_arou_new(distr);

  unur_run_tests( par, RUN_TESTS) ;

  unur_distr_free(distr);

  exit (0);
}

/*---------------------------------------------------------------------------*/







