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

#include <unuran.h>
#include <unuran_distributions.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

#define Q 0.95                /* parameter for geometric distribution */
#define PSIZE       100       /* size of probability vector */
#define SAMPLE_SIZE 100000    /* sample size for tests */
#define N_INTERVALS 1000


/*  #define RUN_TESTS       (~0x0u) */
/*  #define RUN_TESTS       (~0x0u & ~UNUR_TEST_CHI2) */
/*  #define RUN_TESTS       (~0x0u & ~UNUR_TEST_N_URNG) */
#define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER)

/* define which tests should run (1) or not (0) */
#define RUN_DAU           1
#define RUN_DIS           1

#define RUN_UTDR          1
#define RUN_AROU          1
#define RUN_SROU          1
#define RUN_STDR          1
#define RUN_TDRSQRT       1
#define RUN_TDRLOG        1
#define RUN_TABL          1

#define RUN_NORMAL        1
#define RUN_GAMMA         1
#define RUN_BETA          1
#define RUN_CAUCHY        1
#define RUN_UNIFORM       1

#define RUN_RECT          1

#define RUN_CSTD          1

/*---------------------------------------------------------------------------*/

double *set_prob(int len);       /* make a probability vector */

/*---------------------------------------------------------------------------*/

int main()
{ 
  double     *prob;    /* probability vector */
  struct unur_par *par;
/*    struct unur_gen *gen; */
  double fpar[2];
/*    double stp[10]; */
  double slopes[10];
/*    UNUR_URNG_TYPE urng; */

  struct unur_distr *distr_normal;
  struct unur_distr *distr_gamma;
  struct unur_distr *distr_beta;
  struct unur_distr *distr_cauchy;
  struct unur_distr *distr_uniform;

  struct unur_distr *distr_geom;

  struct unur_distr *distr_xxx;

  /* ------------------------- */

  distr_normal = unur_distr_normal(NULL,0);
  //  unur_distr_cont_set_domain(distr_normal,3,UNUR_INFINITY);
  // unur_distr_cont_set_pdfarea(distr_normal,0.01);

  fpar[0] = 3.;
  distr_gamma = unur_distr_gamma(fpar,1);

  fpar[0] = 5.2;
  fpar[1] = 7.9;
/*    fpar[0] = 500.; */
/*    fpar[1] = 790.; */
  distr_beta = unur_distr_beta(fpar,2);

  fpar[0] = 0.;
  fpar[1] = 1.;
  distr_cauchy = unur_distr_cauchy(fpar,2);

  distr_uniform = unur_distr_uniform(NULL,0);
/*    unur_distr_cont_set_mode(distr_uniform,0.); */

  /* ------------------------- */

  prob = set_prob(PSIZE);

  distr_geom = unur_distr_discr_new();
  unur_distr_discr_set_prob(distr_geom,prob,PSIZE);

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Special Generators:\n");
  printf("*******************************************************\n\n");

#if RUN_CSTD == 1

#if 1
  distr_xxx = unur_distr_normal(NULL,0);
  // unur_distr_cont_set_domain(distr_xxx,3,UNUR_INFINITY);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,0);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,1);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,2);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,3);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,4);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,5);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,6);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  distr_xxx = unur_distr_normal(NULL,0);
  par = unur_cstd_new(distr_xxx);
  unur_cstd_set_variant(par,99);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  fpar[0] = 5.;
  distr_xxx = unur_distr_gamma(fpar,1);
  par = unur_cstd_new(distr_xxx);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);

  fpar[0] = 5.;
  fpar[1] = 0.;
  distr_xxx = unur_distr_exponential(fpar,2);
  unur_distr_cont_set_domain(distr_xxx,3,UNUR_INFINITY);
  par = unur_cstd_new(distr_xxx);
  unur_run_tests(par,RUN_TESTS);
  unur_distr_free(distr_xxx);
#endif

/*    fpar[0] = 7.; */
/*    fpar[1] = 5.; */
/*    distr_xxx = unur_distr_beta(fpar,2); */
/*    par = unur_cstd_new(distr_xxx); */
/*    unur_run_tests(par,RUN_TESTS); */
/*    unur_distr_free(distr_xxx); */

  fpar[0] = 0.5;
  fpar[1] = 0.2;
  distr_xxx = unur_distr_beta(fpar,2);
  par = unur_cstd_new(distr_xxx);
  unur_run_tests(par,RUN_TESTS);
/*    gen = unur_init(par); */
/*    { */
/*      int i; */
/*      for (i=0;i<10;i++) */
/*        printf("\n\nx = %g\n\n",unur_sample_cont(gen)); */
/*    } */


#endif

  /* ------------------------- */

  printf("\n*******************************************************\n");
  printf("Discrete Distributions:\n");
  printf("*******************************************************\n\n");

#if RUN_DAU == 1

  printf("DISTRIBUTION:\ttruncated geometric distribution. q = %g, len = %d\n",Q,PSIZE);

  /* get default parameters for new generator */
  par = unur_dau_new(distr_geom);

  /* change defaults */
  unur_dau_set_urnfactor(par,2.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n*******************************************************\n\n");

#endif

  /* ------------------------- */

#if RUN_DIS == 1

  printf("DISTRIBUTION:\ttruncated geometric distribution. q = %g, len = %d\n",Q,PSIZE);

  /* get default parameters for new generator */
  par = unur_dis_new(distr_geom);

  /* change defaults */
  unur_dis_set_variant(par,1);
  /*    unur_set_factor(par,1.); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  par = unur_utdr_new(distr_normal);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  par = unur_arou_new(distr_normal);
  unur_arou_set_cpoints(par,30,NULL);
  unur_arou_set_max_sqhratio(par,1.);
  unur_arou_set_usecenter(par,0);
/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_SROU == 1

  par = unur_srou_new(distr_normal);

  unur_srou_set_Fmode(par,unur_distr_cont_cdf(distr_normal,unur_distr_cont_get_mode(distr_normal)));
  unur_srou_set_usesqueeze(par,1);
/*    unur_srou_set_usemirror(par,1); */
  unur_srou_set_verify(par,1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_STDR == 1

  par = unur_stdr_new(distr_normal);

  unur_stdr_set_Fmode(par,unur_distr_cont_cdf(distr_normal,unur_distr_cont_get_mode(distr_normal)));
  unur_stdr_set_usesqueeze(par,1);
/*    unur_stdr_set_verify(par,1); */
  

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  par = unur_tdr_new(distr_normal);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);
/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  par = unur_tdr_new(distr_normal);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1

  par = unur_tabl_new(distr_normal);
  unur_tabl_set_boundary(par,-50.,50.);

/*    slopes[0] = slopes[2] = 0.; */
/*    slopes[1] = -50.; */
/*    slopes[3] = 50.; */

/*    unur_tabl_set_variant(par,1UL); */
  unur_tabl_set_areafraction(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  par = unur_utdr_new(distr_gamma);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(distr_gamma);
  unur_arou_set_cpoints(par,30,NULL);
  unur_arou_set_max_sqhratio(par,0.99);
  unur_arou_set_usecenter(par,0);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_SROU == 1

  par = unur_srou_new(distr_gamma);

  unur_srou_set_Fmode(par,unur_distr_cont_cdf(distr_gamma,unur_distr_cont_get_mode(distr_gamma)));
  unur_srou_set_usesqueeze(par,1);
/*    unur_srou_set_usemirror(par,1); */
  unur_srou_set_verify(par,1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_STDR == 1

  par = unur_stdr_new(distr_gamma);

  unur_stdr_set_Fmode(par,unur_distr_cont_cdf(distr_gamma,unur_distr_cont_get_mode(distr_gamma)));
  unur_stdr_set_usesqueeze(par,1);
/*    unur_stdr_set_verify(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_gamma);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_gamma);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1

  fpar[0] = ALPHA;
  par = unur_tabl_new(distr_gamma);
  unur_tabl_set_boundary(par,0.,50.);
  unur_tabl_set_variant(par,1UL);
  slopes[0] = slopes[2] = unur_distr_cont_get_mode(distr_gamma);
  slopes[1] = 0.;
  slopes[3] = 50.;
  unur_tabl_set_slopes(par,slopes,2);
  unur_tabl_set_areafraction(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  par = unur_utdr_new(distr_beta);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(distr_beta);
  unur_arou_set_cpoints(par,30,NULL);
  unur_arou_set_max_sqhratio(par,0.);
  unur_arou_set_usecenter(par,0);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_SROU == 1

  par = unur_srou_new(distr_beta);

  unur_srou_set_Fmode(par,unur_distr_cont_cdf(distr_beta,unur_distr_cont_get_mode(distr_beta)));
  unur_srou_set_usesqueeze(par,1);
/*    unur_srou_set_usemirror(par,1); */
  unur_srou_set_verify(par,1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_STDR == 1

  par = unur_stdr_new(distr_beta);
  unur_stdr_set_Fmode(par,unur_distr_cont_cdf(distr_beta,unur_distr_cont_get_mode(distr_beta)));
  unur_stdr_set_usesqueeze(par,1);
/*    unur_stdr_set_verify(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_beta);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_beta);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1

  par = unur_tabl_new(distr_beta);

/*    unur_set_variant(par,1UL); */
  unur_tabl_set_areafraction(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  par = unur_utdr_new(distr_cauchy);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  /* get default parameters for new generator */
  par = unur_arou_new(distr_cauchy);
  unur_arou_set_cpoints(par,30,NULL);
  unur_arou_set_max_sqhratio(par,0.);
  unur_arou_set_usecenter(par,0);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_SROU == 1

  par = unur_srou_new(distr_cauchy);

  unur_srou_set_Fmode(par,unur_distr_cont_cdf(distr_cauchy,unur_distr_cont_get_mode(distr_cauchy)));
  unur_srou_set_usesqueeze(par,1);
/*    unur_srou_set_usemirror(par,1); */
  unur_srou_set_verify(par,1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_STDR == 1

  par = unur_stdr_new(distr_cauchy);

  unur_stdr_set_Fmode(par,unur_distr_cont_cdf(distr_cauchy,unur_distr_cont_get_mode(distr_cauchy)));
  unur_stdr_set_usesqueeze(par,1);
/*    unur_stdr_set_verify(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_cauchy);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  /* get default parameters for new generator */
  par = unur_tdr_new(distr_cauchy);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1

  par = unur_tabl_new(distr_cauchy);
  unur_tabl_set_boundary(par,-1000.,1000.);
    
/*    unur_set_variant(par,1UL); */
  unur_tabl_set_areafraction(par,0.1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  par = unur_utdr_new(distr_uniform);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_AROU == 1

  par = unur_arou_new(distr_uniform);
  unur_arou_set_cpoints(par,30,NULL);
  unur_arou_set_max_sqhratio(par,1.);
  unur_arou_set_usecenter(par,0);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_SROU == 1

  par = unur_srou_new(distr_uniform);

  unur_srou_set_Fmode(par,unur_distr_cont_cdf(distr_uniform,unur_distr_cont_get_mode(distr_uniform)));
  unur_srou_set_usesqueeze(par,1);
/*    unur_srou_set_usemirror(par,1); */
  unur_srou_set_verify(par,1);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_STDR == 1

  par = unur_stdr_new(distr_uniform);

  unur_stdr_set_Fmode(par,unur_distr_cont_cdf(distr_uniform,unur_distr_cont_get_mode(distr_uniform)));
  unur_stdr_set_usesqueeze(par,1);
/*    unur_stdr_set_verify(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRSQRT == 1

  par = unur_tdr_new(distr_uniform);
  unur_tdr_set_c(par,-0.5);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,1.);
  unur_tdr_set_usemode(par,0);

/*    unur_set_debug(par,1); */

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TDRLOG == 1

  par = unur_tdr_new(distr_uniform);
  unur_tdr_set_c(par,0.);
  unur_tdr_set_cpoints(par,30,NULL);
  unur_tdr_set_max_sqhratio(par,0.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

  printf("\n------------------------------------------------------------\n\n");

#endif

#if RUN_TABL == 1

  par = unur_tabl_new(distr_uniform);
  unur_tabl_set_boundary(par,-1.5,2.5);
  unur_tabl_set_areafraction(par,0.25);
  unur_tabl_set_max_sqhratio(par,1.);

  /* run tests */
  unur_run_tests(par,RUN_TESTS);

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
  unur_run_tests(par,RUN_TESTS);

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

