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

  /* choose method */
  par = unur_arou_new(distr);

  /* Set uniform generator in parameter object */
  unur_set_urng(par, ug);

  /* make generator object */
  gen = unur_init(par);

  /* sample: print 100 random numbers */
  for (i=0; i<100; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }

  /* destroy generator object */
  unur_free(gen);

  return 0;
}
