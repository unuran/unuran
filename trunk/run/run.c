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

#define RUN_TESTS       (~0x0u)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;
  struct unur_slist *mlist;

  double time,time_0;

  double unur_urng_mstd(void);

  double fpm[] = {1,5};

  unur_set_default_debug(~0u);
/*    unur_set_stream(stdout); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new(distr); */
/*    unur_tdr_set_variant_ia(par); */

/*    distr = unur_str2distr("normal"); */
/*    par = _unur_str2par(distr,"method=tdr; variant_ia",&mlist); */

  gen = unur_str2gen("normal; domain=(0.,inf); orderstatistics=(10,7) & \
		      method = arou; max_sqhratio = 0.; cpoints=5");


/*    distr = unur_str2distr("normal; domain=(0.,inf); orderstatistics=(10,7)"); */
/*    par = _unur_str2par(distr,"method = arou; max_sqhratio = 0.; cpoints=5",&mlist); */
/*    gen = unur_init(par); */
/*    _unur_slist_free(mlist); */

  for (i=0; i<10; i++)
    printf("%g\n",unur_sample_cont(gen));

  unur_free(gen);

/*    par = unur_unif_new(NULL); */

/*    unur_set_urng(par,unur_urng_mstd); */

/*    gen = unur_init(par); */
/*    unur_srou_chg_domain(gen,0.9,0.95); */
/*    unur_srou_upd_pdfarea(gen); */
/*    unur_srou_upd_mode(gen); */
/*    unur_srou_reinit(gen);   */

/*    unur_test_chi2(gen, 10, 100, 5, 1, stdout); */

/*    unur_run_tests(par,RUN_TESTS); */

/*    _unur_slist_free(mlist); */

  return 0;

}

/*---------------------------------------------------------------------------*/













