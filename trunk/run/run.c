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

/*---------------------------------------------------------------------------*/

double pmf(int k, const UNUR_DISTR *distr)
{
  double fpm[]={1,2,3,4,5,6,5,3,6,7,5,4,3,4,6,7,5,3,4,3,8,6,5,4,7,7,2,2,1,4,7,4,3,7,3,2};

  if (k<0 || k>35) return 0.;
  else return fpm[k];
}

/*---------------------------------------------------------------------------*/

double cdf(int k, const UNUR_DISTR *distr)
{
  double fpm[]={1,2,3,4,5,6,5,3,6,7,5,4,3,4,6,7,5,3,4,3,8,6,5,4,7,7,2,2,1,4,7,4,3,7,3,2};
  double sum;
  int j;

  if (k<0) return 0.;
  if (k>=35) return 1.;

  for (sum=0.,j=0; j<=k; j++)
    sum += fpm[j];

  return sum/158.; 
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;
/*    double x; */

  double fpm[]={1,2,3,4,5,6,5,3,6,7,5,4,3,4,6,7,5,3,4,3,8,6,5,4,7,7,2,2,1,4,7,4,3,7,3,2};

  unur_set_default_debug(~0u);

  distr = unur_distr_discr_new();
  unur_distr_discr_set_domain(distr,0,35);
/*    unur_distr_discr_set_pv(distr,fpm,36); */
/*    unur_distr_discr_set_pmf(distr,pmf); */
  unur_distr_discr_set_cdf(distr,cdf);
  unur_distr_discr_set_pmfsum(distr,158);

  par = unur_dss_new(distr);
  gen = unur_init(par);

  for (i=0;i<1000;i++)
    printf("%d ",unur_sample_discr(gen));
  printf("\n");

  unur_distr_free(distr);
  unur_free(gen);

/*    unur_run_tests(par,UNUR_TEST_CHI2); */
/*    unur_run_tests(par,RUN_TESTS); */

/*    unur_test_chi2(gen,100,0,20,3,stdout); */

  return 0;

}

/*---------------------------------------------------------------------------*/

