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
  int i,k;
  double vec[3],sum=0.,x,hx,cx;

  double fpm[10]={500.,0.9};
  double fpmh[10]={10000.,9000.,500.};

  UNUR_DISTR *distr;    /* distribution */
  UNUR_DISTR *distrh;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG *urng;      

  distr = unur_distr_cauchy(NULL,0);
  unur_distr_cont_set_domain(distr,0.1,1.);
  unur_distr_cont_upd_mode(distr);
  unur_distr_cont_upd_pdfarea(distr);

  par = unur_cstd_new(distr);
  fpm[0] = 2.;
  fpm[1] = 5.;
  gen = unur_init(par);
  unur_cstd_chg_pdfparams(gen,fpm,2);

  unur_test_chi2( gen, 100, 0, 0, 1 );

  unur_free(gen);


#if 0


  distrh = unur_distr_hypergeometric(fpmh,3);

  par = unur_dstd_new(distrh);
  if (!unur_dstd_set_variant(par,1)) { par = NULL; }
  gen = unur_init(par);
  for (i=0; i<500; i++)
  { k=unur_sample_discr(gen);
    printf("%d\n",k);
  }


 
  distr = unur_distr_binomial(fpm,2);

  for(i=0;i<=500;i++)
  { x= unur_distr_discr_eval_pmf(i,distr); 
    sum+= x;
    cx= unur_distr_discr_eval_cdf(i,distr); 
    hx= unur_distr_discr_eval_pmf(i,distrh); 
    printf("pmf(%d)=%e h: %e diff %e sum= %e; cdf(%d)=%e\n",i,x,hx,x-hx,sum,i,cx);
  }
  printf("mode hyp: %d \n",unur_distr_discr_get_mode(distrh)); 


   unur_test_chi2( gen, 100, 0,  0, 3 );
#endif  

  /*

  urng = prng_new("mt19937(2345)");
  if (!urng) exit(-1);

  unur_set_default_urng(urng);

  distr = unur_distr_normal(NULL,0);
  unur_distr_cont_set_domain(distr,1,20);

  par = unur_tdr_new(distr);
  unur_tdr_set_variant_ia(par);

  gen = unur_init(par);

  for (i=0; i<10; i++)
    unur_sample_cont(gen);

   unur_tdr_chg_truncated(gen,5.,20000.);

   unur_test_chi2( gen, 100, 0,  0, 2 );


  // unur_run_tests( par, RUN_TESTS) ;

  */


  unur_distr_free(distr);
  exit (0);
}

/*---------------------------------------------------------------------------*/














