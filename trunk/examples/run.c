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

/*---------------------------------------------------------------------------*/

int main()
{

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  double data[] = {2.,6.,1.,9.,15.,14.,4.,3.,8.,11.,5.,0.,7.};

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);
  
  distr = unur_distr_new(UNUR_DISTR_CEMP);
  unur_distr_cemp_set_data( distr, data, 13 );

  par = unur_empk_new( distr );
  unur_empk_set_kernel(par, UNUR_DISTR_GAUSSIAN);

  gen = unur_init( par );

  unur_test_printsample(gen,10,8);
  

  exit (0);
}

/*---------------------------------------------------------------------------*/







