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

  double fpar[2] = {2.,2.};
  double cp[5] = {0.1, 0.5,0.6,2.,10.};
  
  UNUR_URNG *urng, *urng_aux;

  urng = prng_new("mt19937(2345)");
  urng_aux = prng_new("mt19937(9987)");
  unur_set_default_urng(urng);
  

  distr = unur_distr_normal(NULL,0);
/*    unur_set_default_debug(1); */

/*    distr = unur_distr_beta(fpar,2); */
/*    unur_distr_cont_set_domain(distr,-UNUR_INFINITY,UNUR_INFINITY); */
/*    unur_distr_cont_set_mode(distr,0.5); */
/*    distr = unur_distr_normal(fpar,2); */

  par = unur_tdr_new(distr);
  unur_set_urng_aux(par,urng_aux);
  unur_tdr_set_variant_ps(par);
/*    unur_tdr_set_cpoints(par,5,cp); */
/*    unur_tdr_set_max_sqhratio(par,0.); */
/*    unur_tdr_set_verify(par,1); */
  unur_tdr_set_c(par,0.);
  unur_run_tests(par,~0);
  
/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ia(par); */
/*    unur_tdr_set_cpoints(par,10,NULL); */
/*    unur_run_tests(par,~0); */
  

  /* choose method */


/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_gw(par); */
/*    unur_tdr_set_cpoints(par,40,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_gw(par); */
/*    unur_tdr_set_c(par,0.); */
/*    unur_tdr_set_cpoints(par,40,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ps(par); */
/*    unur_tdr_set_cpoints(par,4,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ps(par); */
/*    unur_tdr_set_c(par,0.); */
/*    unur_tdr_set_cpoints(par,4,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ia(par); */
/*    unur_tdr_set_cpoints(par,40,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ia(par); */
/*    unur_tdr_set_c(par,0.); */
/*    unur_tdr_set_cpoints(par,40,NULL); */
/*    unur_run_tests(par,~0); */

/*    par = unur_arou_new(distr); */
/*    unur_arou_set_cpoints(par,40,NULL); */
/*    unur_run_tests(par,~0); */
  
  exit (0);
}

/*---------------------------------------------------------------------------*/







