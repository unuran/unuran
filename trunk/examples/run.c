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

int main()
{
  int i,j;
  double vec[3];

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  UNUR_DISTR *uv;
  UNUR_GEN *uvgen;

  double mean[5] = { 1.,2.,3.};

  double covar[25] = {  1.0, 0.5, 0.6,
		        0.5, 2.0, 0.7,
		        0.6, 0.7, 3.0};


  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

/*    uv = unur_distr_cauchy(NULL,0); */
/*    uvgen = unur_init( unur_tdr_new( uv )); */
/*    unur_distr_free(uv); */

  distr = unur_distr_multinormal(3,mean,covar);
  par = unur_vmt_new(distr);
/*    unur_vmt_set_marginalgen( par, uvgen ); */
/*    unur_run_tests(par,RUN_TESTS); */
  gen = unur_init(par);

  for (i=0; i<20; i++) {
    unur_sample_vec(gen,vec);
    printf("%g, %g, %g\n",vec[0],vec[1],vec[2]);
  }

  unur_distr_free(distr);
  unur_free(gen);
/*    unur_free(uvgen); */

  exit (0);
}

/*---------------------------------------------------------------------------*/







