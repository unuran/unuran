/* my fourth UNURAN program: test4.c*/


#include <unuran.h>
int main()
{
  int    i;
  double x;
  double params[2] = {10.0, 0.5};

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG  *ug;       /* uniform generator */

  /* choose kind of uniform generator: Mersenne Twister */
  ug =  prng_new("mt19937(1237)");

  /* choose a implemented distribution: Gaussian */
  distr = unur_distr_normal(params, 2);
  //unur_distr_cont_set_domain(distr, -1, 0);

  /* choose method */
  par = unur_ninv_new(distr);
  //unur_ninv_use_newton(par);

  /* Set uniform generator in parameter object */   
  unur_ninv_set_table(par, 50);
  unur_set_urng(par, ug);    

  /* make generator object */
  gen = unur_init(par);
  
  /* sample: print 100 random numbers */
  for (i=0; i<6; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }
  unur_ninv_chg_max_iter(gen,10);
  unur_ninv_chg_table(gen, 100);

  fprintf(stderr, "vorherr\n");
  unur_ninv_chg_table(gen, 500);
  fprintf(stderr, "nachherr\n");
 
  for (i=0; i<8; i++) {
    x = unur_sample_cont(gen);
    printf("1--%f\n",x);
  }

  for (i=0; i<8; i++) {
    x = unur_sample_cont(gen);
    printf("2--%f\n",x);
  }

  for (i=0; i<8; i++) {
    x = unur_sample_cont(gen);
  }

  printf("Und nun mit veraenderten Parametern:\n");
  params[0] = 0.;
  params[1] = 1.;
  unur_ninv_chg_pdfparams(gen, params, 2);
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    printf("3--%f\n",x);
  }

  printf("vorher\n");
  unur_ninv_chg_table(gen, 1000);
  printf("nachher\n");
  
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    printf("4--%f\n",x);
  }


  return 0;
}













