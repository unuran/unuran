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
  int i;
  double vec[3];

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  double data[] = {1.1, 2.1, 2.3,
		   0.9, 2.0, 2.0,
		   1.2, 1.8, 2.5,
		   1.5, 2.5, 1.7,
		   0.,  0.5, 0.8,
		   1.8, 2.8, 0. };

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  distr = unur_distr_cvemp_new(3);
  unur_distr_cvemp_set_data( distr, data, 6 );

  par = unur_vempk_new(distr);
  gen = unur_init(par);

  for (i=0; i<20; i++) {
    unur_sample_vec(gen,vec);
    printf("%g, %g, %g\n",vec[0],vec[1],vec[2]);
  }

  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}

/*---------------------------------------------------------------------------*/







