/* my fourth UNURAN program: test4.c*/

#include <unuran.h>

int main()
{
  int    i;
  double x;

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */
  UNUR_URNG  *ug;       /* uniform generator */

  /* choose kind of uniform generator: Mersenne Twister */
  ug =  prng_new("mt19937(1237)");

  /* choose a implemented distribution: Gaussian */
  distr = unur_distr_normal(NULL, 0);
  unur_distr_cont_set_domain(distr, -1, 0);

  /* choose method */
  par = unur_ninv_new(distr);
  unur_ninv_use_newton(par);

  /* Set uniform generator in parameter object */   
  unur_set_urng(par, ug);    

  /* make generator object */
  gen = unur_init(par);

  /* sample: print 100 random numbers */
  for (i=0; i<6; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }
  unur_ninv_chg_x_resolution(gen, 1e-2); 
  unur_ninv_chg_max_iter(gen,10);
  unur_ninv_chg_domain(gen, 0.,1.);
  unur_ninv_chg_start(gen, 1., 1.);

  for (i=0; i<6; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }


  return 0;
}













