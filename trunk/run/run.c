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

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  double fpm[] = { 2. };

  unur_set_default_debug(~0u);


/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new( distr ); */
/*    gen = unur_init( par ); */
/*    unur_distr_free(distr); */
/*    unur_free(gen); */

  distr = unur_distr_powerexponential(fpm,1);
  par = unur_tdr_new( distr );
/*    unur_tdr_set_variant_gw(par); */
  gen = unur_init( par );
  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}

/*---------------------------------------------------------------------------*/














