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

void _unur_distr_corder_debug( struct unur_distr *os, char *genid );
void _unur_distr_cont_debug( struct unur_distr *distr, char *genid );

#define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER)

/*---------------------------------------------------------------------------*/

int main()
{

  UNUR_DISTR *distr, *os;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  double data[] = {2.,6.,1.,9.,15.,14.,4.,3.,8.,11.,5.,0.,7.};

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);
  
  distr = unur_distr_normal(NULL,0);
  os = unur_distr_corder_new(distr,10,5);
  // unur_distr_cont_set_mode(os,0.5);
  // unur_distr_cont_set_mode(os,0.);

  // _unur_distr_cont_debug( os, "junk" );

  unur_distr_cont_set_domain(os,0.,1.);
  //  unur_distr_cont_upd_pdfarea(os);

  par = unur_tdr_new(os);
  unur_run_tests(par,RUN_TESTS);

  printf("end\n");

  /*
  distr = unur_distr_new(UNUR_DISTR_CEMP);
  unur_distr_cemp_set_data( distr, data, 13 );

  par = unur_empk_new( distr );
  unur_empk_set_kernel(par, UNUR_DISTR_GAUSSIAN);

  gen = unur_init( par );

  unur_test_printsample(gen,10,8);
  */

  exit (0);
}

/*---------------------------------------------------------------------------*/







