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

#include <mcgibbs.h>


#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

/*---------------------------------------------------------------------------*/

int main()
{
/*   UNUR_DISTR *distr; */
  UNUR_PAR *par;
  UNUR_GEN *gen;

  UNUR_DISTR *normal;

  unur_set_default_debug(~0U);

  normal = unur_distr_multinormal(3,NULL,NULL);

  par = unur_mcgibbs_new(normal);
  unur_set_debug( par, 1u);
/*   unur_mcgibbs_set_thinning(par,3); */
/*   unur_mcgibbs_set_variant_random_direction(par); */

/*   gen = unur_init(par); */

  unur_run_tests(par,~0u);
  
/*   unur_free(gen); */
  unur_distr_free(normal);

  return 0;
}

/*---------------------------------------------------------------------------*/

