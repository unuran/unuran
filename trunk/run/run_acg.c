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

#include <unuran.h>
#include <unuran_acg.h>
#include <unuran_tests.h>

#define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER)

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  unur_set_default_debug(~0u);

  distr = unur_distr_normal(NULL,0);
  par = unur_tdr_new( distr );
  gen = unur_init( par );
  unur_acg( gen, stdout, NULL );
  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}

/*---------------------------------------------------------------------------*/














