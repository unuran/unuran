/* my third UNURAN program example3.c */

#include <unuran.h>

int main()
{
  int    i;
  double x;
  double params[2] = {10.0, 0.5};

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */

  /* choose Gaussian distribution with parameters 10.0 and 0.5 */
  distr = unur_distr_normal(params, 2);

  /* choose method */
  par = unur_ninv_new(distr);

  /* change a parameter of the used method */
  unur_ninv_use_newton(par);

  /* make generator object */
  gen = unur_init(par);

  /* sample: print 100 random numbers */
  for (i=0; i<100; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }

  /* destroy distribution- and generator object */
  unur_distr_free(distr);
  unur_free(gen);

  return 0;
}
