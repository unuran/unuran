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

  distr = unur_distr_normal(NULL,0);
  unur_set_default_debug(1);


  /* choose method */

  par = unur_tdr_new(distr);
  unur_tdr_set_version_gw(par);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_tdr_new(distr);
  unur_tdr_set_version_gw(par);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_tdr_new(distr);
  unur_tdr_set_version_ps(par);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_tdr_new(distr);
  unur_tdr_set_version_ps(par);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_tdr_new(distr);
  unur_tdr_set_version_ia(par);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_tdr_new(distr);
  unur_tdr_set_version_ia(par);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);

  par = unur_arou_new(distr);
  unur_arou_set_cpoints(par,40,NULL);
  unur_run_tests(par,~0);
  
  exit (0);
}

/*---------------------------------------------------------------------------*/







