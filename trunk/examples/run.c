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
  //  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_debug(0);
  unur_set_default_urng(urng);
  
  distr = unur_distr_normal(NULL,0);

  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_gw(par);
  unur_run_tests(par,~0u);
  
  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_ps(par);
  unur_run_tests(par,~0u);
  
  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_ia(par);
  unur_run_tests(par,~0u);
  
  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_gw(par);
  unur_run_tests(par,~0u);
  
  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_ps(par);
  unur_run_tests(par,~0u);
  
  par = unur_tdr_new(distr);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_tdr_set_variant_ia(par);
  unur_run_tests(par,~0u);
  

  exit (0);
}

/*---------------------------------------------------------------------------*/







