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

/* pdf with piecewise linear function as transformed density with T = -1/sqrt */
double pdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
	double y = 1./(fabs(x)+1.);
	return y*y;
}
double dpdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
	double y = 1./(fabs(x)+1.);
	y = 2.*y*y*y;
	return ((x<0.) ? y : - y);
}
double cdf_sqrtlin( double x, UNUR_DISTR *distr )
{ 
	if (x<=0.)
		return 0.5/(1.-x);
	else
		return (1.-0.5/(1.+x));
}

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  double fpm[] = { 2., 3. };

  unur_set_default_debug(~0u);
  unur_set_stream(stdout);

  distr = unur_distr_cont_new();
  unur_distr_cont_set_pdf(distr,pdf_sqrtlin);
  unur_distr_cont_set_dpdf(distr,dpdf_sqrtlin);
  unur_distr_cont_set_cdf(distr,cdf_sqrtlin);
  unur_distr_set_name(distr,"sqrtlin");

  unur_errno = 0;
  par = unur_tdr_new(distr);
  unur_tdr_set_variant_gw(par);
  unur_tdr_set_c(par,-0.5);
  gen = unur_init(par);

  unur_free(gen);


/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new( distr ); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,1u); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,3u); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,UNUR_STDGEN_DEFAULT); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_exponential(NULL,0); */
/*    par = unur_cstd_new( distr ); */
/*    unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*    distr = unur_distr_powerexponential(fpm,1); */
/*    par = unur_tdr_new( distr ); */
/*    unur_tdr_set_variant_gw(par); */
/*    gen = unur_init( par ); */
/*    unur_distr_free(distr); */
/*    unur_free(gen); */

  exit (0);
}

/*---------------------------------------------------------------------------*/














