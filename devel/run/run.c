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
#include <testdistributions.h>

#include <ptx.h>

#define MATHLIB_STANDALONE 1
#include <Rmath.h>

#define RUN_TESTS       (~0x0u)
/* #define RUN_TESTS       UNUR_TEST_SAMPLE */

#define unur_distr_multinormal  unur_distr_multinormal_w_marginals


/*---------------------------------------------------------------------------*/

/* double pdf( double x, const UNUR_DISTR *distr ) */
/* { */
/* /\*   return (1.e32 * pow(1-x*x,-50) * exp((78.043 * x - 107.415)/(1 - x*x))); *\/ */
/*   return exp(log(1.e32)+  50*log(1-x*x) + ((78.043 * x - 107.415)/(1 - x*x))); */
/* } */

double cdft( double x, const UNUR_DISTR *distr )
{
  const double *df;
  unur_distr_cont_get_pdfparams(distr, &df);
  return pt(x,df[0],TRUE,FALSE);
}

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *d1, *d2;
  UNUR_PAR *par;
  UNUR_GEN *gen;

  UNUR_GEN *tgen, *xgen;

  int n=10;
  double u,du;
  double udiff;
  double uin,uout;
  double du1, du2;
  
  
  int i;
  double t,x;
  
  double df1[] = {5.};
  double df2[] = {10.};

  double max_u_error, max_x_error, max_xr_error;
  double max_x1_error, max_x2_error, max_x3_error;

  unur_set_default_debug(~0U);

  d1 = unur_distr_student(df1,1);
  d2 = unur_distr_student(df2,1);

  unur_distr_cont_set_cdf(d1,cdft);
  unur_distr_cont_set_cdf(d2,cdft);

  par = unur_pinv_new(d1);
  unur_pinv_set_u_resolution(par,1.e-15);
  tgen = unur_init(par);

  par = unur_pinv_new(d2);
  unur_pinv_set_u_resolution(par,1.e-12);
  xgen = unur_init(par);

  par = unur_ptx_new(d1,d2);
  unur_ptx_set_u_resolution(par,1e-12);
/*   unur_ptx_set_usecdf(par); */
  gen = unur_init(par);


  n = 33;
  du = 1./n;

  for (i=0; i<n; i++) {
    u = (i+0.5)*du;
    t = unur_quantile(tgen,u);
    x = unur_ptx_eval_approxinvcdf(gen,t);
    udiff = fabs(u-unur_distr_cont_eval_cdf(x,d2));
    printf("u=%10g, t=%13g, x=%13g, udiff=%13g\n",u,t,x,udiff);
  }

/*   for (i=0; i<10000000; i++) { */
/*     t = unur_sample_cont(tgen); */
/* #if 0 */
/*     x = unur_ptx_eval_approxinvcdf(gen,t); */
/* #else */
/*     u = unur_distr_cont_eval_cdf(t,d1); */
/*     x = unur_quantile(xgen,u); */
/* #endif */
/*   } */
  
/*   maxerror = 0.; */
/*   for (i=0; i<1000000; i++) { */
/*     t = unur_sample_cont(tgen); */
/*     x = unur_ptx_eval_approxinvcdf(gen,t); */
/*     u = unur_distr_cont_eval_cdf(t,d1); */
/*     x = fabs(x-unur_quantile(xgen,u)); */
/*     if (x>maxerror) maxerror=x; */
/*   } */
/*   printf("max x-error=%g\n",maxerror); */

  n = 1000000;
  max_u_error = 0.;
  for (i=0; i<n; i++) {
    u = (i+0.5) / n;
    t = unur_quantile(tgen,u);
    x = unur_ptx_eval_approxinvcdf(gen,t);
    u = fabs(unur_distr_cont_eval_cdf(t,d1) -
	     unur_distr_cont_eval_cdf(x,d2));
    if (u>max_u_error) max_u_error=u;
  }
  printf("max u-error=%g\n",max_u_error);

#if 0
  n = 10000000;
  max_u_error = 0.;
  max_x1_error = 0.;
  max_x2_error = 0.;
  max_x3_error = 0.;
  for (i=0; i<n; i++) {
    double x1,x2,dx1,dx2,dx3;
    u = (i+0.5) / n;
    t = unur_quantile(tgen,u);
    uin = unur_distr_cont_eval_cdf(t,d1);

    x = unur_ptx_eval_approxinvcdf(gen,t);
/*     uout = unur_distr_cont_eval_cdf(x,d2); */
/*     du2 = fabs(uout-uin); */
/*     if (du2>max_u_error) max_u_error=du2; */

    x1 = unur_quantile(xgen,uin);
    x2 = qt(uin,df2[0],TRUE,FALSE);
    dx1 = fabs(x1-x);
    dx2 = fabs(x2-x);
    dx3 = fabs(x1-x2);
    if (dx1>max_x1_error) max_x1_error=dx1;
    if (dx2>max_x2_error) max_x2_error=dx2;
    if (dx3>max_x3_error) max_x3_error=dx3;

    if (1 && ( dx1 > 1.e-8 || dx2 > 1.e-8 || dx3 > 1.e-8) )
      printf("u=%10g, t=%13g, x=%13g, dx1=%13g, dx2=%13g, dx3=%13g\n",
	     u,t,x,dx1,dx2,dx3);

  }
  printf("max  u-error = %g\n",max_u_error);
  printf("max x1-error = %g\n",max_x1_error);
  printf("max x2-error = %g\n",max_x2_error);
  printf("max x3-error = %g\n",max_x3_error);
#endif
  
  unur_distr_free(d1);
  unur_distr_free(d2);
  unur_free(tgen);
  unur_free(xgen);
  unur_free(gen);

  return 0;
}

/*---------------------------------------------------------------------------*/

