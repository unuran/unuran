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
  double fpar[2] = {-1., 5.};

  unur_set_default_debug(~0u);

  distr = unur_distr_normal(fpar,2);
  unur_distr_cont_set_domain(distr,0,1);
  par = unur_tdr_new( distr );
/*    unur_tdr_set_c(par,0.); */
  unur_tdr_set_cpoints(par,4,NULL),
  gen = unur_init( par );
/*    unur_acg_UNURAN( gen, stdout, NULL ); */
/*    unur_acg_C( gen, stdout, NULL ); */
  unur_acg_FORTRAN( gen, stdout, NULL );
/*    unur_acg_JAVA( gen, stdout, NULL ); */
  unur_distr_free(distr);
  unur_free(gen);

  exit (0);
}

/*---------------------------------------------------------------------------*/














