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

#include <unur_methods.h>
#include <unur_distr.h>
#include <unur_tests.h>

/*---------------------------------------------------------------------------*/

#define Q 0.95                /* parameter for geometric distribution */
#define PSIZE       100       /* size of probability vector */
#define SAMPLE_SIZE 100000    /* sample size for tests */
#define N_INTERVALS 1000


#define RUN_TESTS       (~0x0UL)
/*  #define RUN_TESTS       (~0x0UL & ~UNUR_TEST_CHI2) */
/*  #define RUN_TESTS       (~0x0UL & ~UNUR_TEST_N_URNG) */
/*  #define RUN_TESTS       (~0x0UL & ~UNUR_TEST_SCATTER) */


/* define which tests should run (1) or not (0) */
#define RUN_DAU           0
#define RUN_DIS           0

#define RUN_UTDR          0
#define RUN_AROU          0
#define RUN_TDRSQRT       0
#define RUN_TDRLOG        1
#define RUN_TABL          0

#define RUN_NORMAL        1
#define RUN_GAMMA         1
#define RUN_BETA          1
#define RUN_CAUCHY        1
#define RUN_UNIFORM       1

#define RUN_RECT          0

#define RUN_CSTD          0

/*---------------------------------------------------------------------------*/

double *set_prob(int len);       /* make a probability vector */

/*---------------------------------------------------------------------------*/

int main()
{ 
  double     *prob;    /* probability vector */
  struct unur_par *par;
  double fpar[2];
/*    double stp[10]; */
  double slopes[10];
/*    UNUR_URNG_TYPE urng; */

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Special Generators:\n");
  printf("*******************************************************\n\n");

#if RUN_CSTD == 1

/*    par = unur_cstd_new("  normal (0 , 1 ) "); */
/*    unur_run_tests(par,RUN_TESTS,unur_cdf_normal); */

  par = unur_cstd_new("normal");
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

/*    par = unur_cstd_new("gamma"); */
/*    unur_run_tests(par,RUN_TESTS,unur_cdf_gamma); */

  par = unur_cstd_new("gamma(5)");
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

/*    par = unur_cstd_new("gamma(5,2)"); */
/*    unur_run_tests(par,RUN_TESTS,unur_cdf_gamma); */

/*    par = unur_cstd_new("gamma(5,3,-3)"); */
/*    unur_run_tests(par,RUN_TESTS,unur_cdf_gamma); */

/*    par = unur_cstd_new("exponential"); */
/*    unur_run_tests(par,RUN_TESTS,unur_cdf_exponential); */

  par = unur_cstd_new("exponential(5)");
  unur_run_tests(par,RUN_TESTS,unur_cdf_exponential);


#endif

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Discrete Distributions:\n");
  printf("*******************************************************\n\n");

  /* ------------------------- */
  prob = set_prob(PSIZE);
  /* ------------------------- */


#if RUN_DAU == 1

  printf("DISTRIBUTION:\ttruncated geometric distribution. q = %g, len = %d\n",Q,PSIZE);

  /* get default parameters for new generator */
  par = unur_dau_new(prob,PSIZE);
  /* change defaults */
  unur_set_factor(par,2.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,NULL);

  printf("\n*******************************************************\n\n");

#endif

  /* ------------------------- */

#if RUN_DIS == 1

  printf("DISTRIBUTION:\ttruncated geometric distribution. q = %g, len = %d\n",Q,PSIZE);

  /* get default parameters for new generator */
  par = unur_dis_new(prob,PSIZE);
  /* change defaults */
  /*    unur_set_variant(par,2); */
  /*    unur_set_factor(par,1.); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS,NULL);

  printf("\n*******************************************************\n\n");

#endif

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Continous Distributions:\n");
  printf("*******************************************************\n\n");

  /* ------------------------- */

#if RUN_NORMAL == 1

  printf("DISTRIBUTION:\tstandard normal. (mu = 0., sigma = 1.)\n");

#if RUN_UTDR == 1

  /* get default parameters for new generator */
  par = unur_utdr_new(unur_pdf_normal,0.);
  unur_set_pdf_area(par,unur_area_normal(NULL,0));

  unur_set_pdf_area(par,-3.);


  unur_set_factor(par,3.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  par = unur_arou_new(unur_pdf_normal,unur_dpdf_normal);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,1.);
/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  par = unur_tdr_new(unur_pdf_normal,unur_dpdf_normal);
  unur_set_mode(par,0.);
  unur_set_tdr_c(par,-0.5);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);
/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  par = unur_tdr_new(unur_pdf_normal,unur_dpdf_normal);
  unur_set_mode(par,0.);
  unur_set_tdr_c(par,0.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1
    
  par = unur_tabl_new(unur_pdf_normal);

/*    slopes[0] = slopes[2] = 0.; */
/*    slopes[1] = -50.; */
/*    slopes[3] = 50.; */
/*    unur_set_slopes(par,slopes,2); */

  unur_set_mode(par,0.);
  unur_set_domain(par,-50.,50.);
  unur_set_variant(par,1UL);

/*    unur_set_max_intervals(par,1000); */
/*    unur_set_max_shratio(par,1.); */

  unur_set_pdf_area(par,unur_area_normal(NULL,0));
  unur_set_tabl_c(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_normal);

#endif

#endif

  /* ------------------------- */

#if RUN_GAMMA == 1
#define ALPHA 3.

  fpar[0] = ALPHA;

  printf("\n*******************************************************\n\n");
  printf("DISTRIBUTION:\tgamma. alpha = %g\n",ALPHA);

#if RUN_UTDR == 1

  /* get default parameters for new generator */
  par = unur_utdr_new(unur_pdf_gamma,unur_mode_gamma(fpar,1));
  unur_set_domain(par,0.,UNUR_INFINITY);
  unur_set_pdf_param(par,fpar,1);

/*    unur_set_pdf_area(par,area_gamma(fpar,1)); */
  unur_set_pdf_area(par,ALPHA-1.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(unur_pdf_gamma,unur_dpdf_gamma);
  unur_set_domain(par,-1.,UNUR_INFINITY);
  unur_set_pdf_param(par,fpar,1);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.99);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_gamma,unur_dpdf_gamma);
  unur_set_domain(par,0.,UNUR_INFINITY);
  unur_set_pdf_param(par,fpar,1);
  unur_set_mode(par,unur_mode_gamma(fpar,1));
  unur_set_tdr_c(par,-0.5);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_gamma,unur_dpdf_gamma);
  unur_set_domain(par,0.,UNUR_INFINITY);
  unur_set_pdf_param(par,fpar,1);
  unur_set_mode(par,unur_mode_gamma(fpar,1));
  unur_set_tdr_c(par,0.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1
    
  par = unur_tabl_new(unur_pdf_gamma);
  unur_set_pdf_param(par,fpar,1);
/*    unur_set_domain(par,0.,50.); */
/*    unur_set_mode(par,unur_mode_gamma(fpar,1)); */
/*    unur_set_variant(par,1UL); */
  slopes[0] = slopes[2] = unur_mode_gamma(fpar,1);
  slopes[1] = 0.;
  slopes[3] = 50.;
  unur_set_slopes(par,slopes,2);
/*    unur_set_pdf_area(par,unur_area_gamma(fpar,1)); */
  unur_set_pdf_area(par,ALPHA-1.);
  unur_set_tabl_c(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_gamma);

#endif
#undef ALPHA
#endif

  /* ------------------------- */

#if RUN_BETA == 1
#define A 5.2
#define B 7.9

  fpar[0]=A;
  fpar[1]=B;

  printf("\n*******************************************************\n\n");
  printf("DISTRIBUTION:\tbeta. a = %g, b = %g\n",A,B);

#if RUN_UTDR == 1

  /* get default parameters for new generator */
  par = unur_utdr_new(unur_pdf_beta,unur_mode_beta(fpar,2));
  unur_set_domain(par,0.,1.);
  unur_set_pdf_param(par,fpar,2);

/*    unur_set_pdf_area(par,unur_area_beta(fpar,2)); */
  unur_set_pdf_area(par,0.000217719);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_beta);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(unur_pdf_beta,unur_dpdf_beta);
  unur_set_domain(par,0.,1.);
  unur_set_pdf_param(par,fpar,2);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_beta);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_beta,unur_dpdf_beta);
  unur_set_domain(par,0.,1);
  unur_set_pdf_param(par,fpar,2);
  unur_set_mode(par,unur_mode_beta(fpar,2));
  unur_set_tdr_c(par,-0.5);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_beta);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_beta,unur_dpdf_beta);
  unur_set_domain(par,0.,1);
  unur_set_pdf_param(par,fpar,2);
  unur_set_mode(par,unur_mode_beta(fpar,2));
  unur_set_tdr_c(par,0.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_beta);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1
    
  par = unur_tabl_new(unur_pdf_beta);
  unur_set_pdf_param(par,fpar,2);
  unur_set_domain(par,0.,1.);
  unur_set_mode(par,unur_mode_beta(fpar,2));
/*    unur_set_variant(par,1UL); */
/*    unur_set_pdf_area(par,unur_area_beta(fpar,2)); */
  unur_set_pdf_area(par,0.000217719);
  unur_set_tabl_c(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_beta);

#endif
#undef A
#undef B
#endif

  /* ------------------------- */

#if RUN_CAUCHY == 1
#define theta 0.
#define lambda 1.

  fpar[0] = theta;
  fpar[1] = lambda;

  printf("\n*******************************************************\n\n");
  printf("DISTRIBUTION:\tcauchy. theta = %g, lambda = %g\n",theta,lambda);

#if RUN_UTDR == 1

  /* get default parameters for new generator */
  par = unur_utdr_new(unur_pdf_cauchy,unur_mode_cauchy(fpar,2));
  unur_set_pdf_param(par,fpar,2);
  unur_set_pdf_area(par,unur_area_cauchy(fpar,2));

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_cauchy);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(unur_pdf_cauchy,unur_dpdf_cauchy);
  unur_set_pdf_param(par,fpar,2);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_cauchy);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_cauchy,unur_dpdf_cauchy);
  unur_set_pdf_param(par,fpar,2);
  unur_set_mode(par,unur_mode_cauchy(fpar,2));
  unur_set_tdr_c(par,-0.5);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_cauchy);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(unur_pdf_cauchy,unur_dpdf_cauchy);
  unur_set_pdf_param(par,fpar,2);
  unur_set_mode(par,unur_mode_cauchy(fpar,2));
  unur_set_tdr_c(par,0.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_cauchy);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1
    
  par = unur_tabl_new(unur_pdf_cauchy);
  unur_set_pdf_param(par,fpar,2);
  unur_set_mode(par,0.);
  unur_set_domain(par,-50.,50.);
/*    unur_set_variant(par,1UL); */
  unur_set_pdf_area(par,unur_area_cauchy(fpar,2));
  unur_set_tabl_c(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_cauchy);

  printf("\n*******************************************************\n\n");

#endif
#undef theta
#undef lambda
#endif

  /* ------------------------- */

#if RUN_UNIFORM == 1

  printf("DISTRIBUTION:\tuniform distribution. U(0,1)\n");

#if RUN_UTDR == 1

  /* get default parameters for new generator */
  par = unur_utdr_new(unur_pdf_uniform,0.5);
  unur_set_pdf_area(par,1.);
  unur_set_domain(par,0.,1.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_uniform);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  par = unur_arou_new(unur_pdf_uniform,unur_dpdf_uniform);
  unur_set_domain(par,0.,1.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,1.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_uniform);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  par = unur_tdr_new(unur_pdf_uniform,unur_dpdf_uniform);
/*    unur_set_mode(par,0.5); */
  unur_set_domain(par,-0.5,1.5);
  unur_set_tdr_c(par,-0.5);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,1.);
/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_uniform);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  par = unur_tdr_new(unur_pdf_uniform,unur_dpdf_uniform);
  unur_set_mode(par,0.5);
  unur_set_domain(par,0.,1.);
  unur_set_tdr_c(par,0.);
  unur_set_cpoints(par,30,NULL);
  unur_set_max_shratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_uniform);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1
    
  par = unur_tabl_new(unur_pdf_uniform);
  unur_set_mode(par,0.5);
  unur_set_domain(par,-1.5,2.5);
/*    unur_set_variant(par,1UL); */
  unur_set_pdf_area(par,1.);
  unur_set_tabl_c(par,0.25);

  unur_set_max_shratio(par,1.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,unur_cdf_uniform);

#endif

#endif

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Multivariate Distributions:\n");
  printf("*******************************************************\n\n");

  /* ------------------------- */

#if RUN_RECT == 1

#define DIM         3         /* dimension of hypercube */

  printf("\n*******************************************************\n\n");
  printf("DISTRIBUTION:\tuniform over hypercube. dim = %d\n",DIM);

  par = unur_rect_new(DIM);

  /* run tests */
  unur_run_tests(par,RUN_TESTS,NULL);

#undef DIM

#endif

  /* ------------------------- */

  exit(0);

} /* end of main */

/*---------------------------------------------------------------------------*/

double *set_prob(int len)
/* 
   make a probability vector (need not sum to one)
   (use truncated geometric distribution)
*/
{
  double *prob;
  int i;
  extern double uniform();

  /* allocate memory */
  prob = malloc(len*sizeof(double));
  if (!prob) { fprintf(stderr,"allocation failure\n"); abort(); }
     
  /* main part of geometric distribution */
  prob[0] = 1.;
  for( i=1; i<len; i++ ) 
    prob[i] = prob[i-1] * Q;
/*      prob[i] = uniform(); */

  return prob;

} /* end of set_prob() */

/*---------------------------------------------------------------------------*/

