/* my first UNURAN program: example1.c*/

#include <unuran.h>

int main()
{
  int    i;
  double x;

  UNUR_DISTR *distr;    /* distribution */
  UNUR_PAR   *par;      /* parameter */
  UNUR_GEN   *gen;      /* generator */

  /* choose a implemented distribution: Gaussian */
  distr = unur_distr_normal(NULL, 0);

  /* choose method */
  par = unur_arou_new(distr);

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
