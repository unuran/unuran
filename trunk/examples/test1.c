/* my first UNURAN program */ 

#include <unuran.h>
#include <unuran_distributions.h>
#include <unuran_tests.h>



int main()
{
  double x;
  int    i;

  struct unur_distr *distr;  /* distribution */ 
  struct unur_par *par;      /* parameter */
  struct unur_gen *gen;      /* generator */
/* mache neues Verteilungsobjekt unur_distr_cont_new (methods.distr.c)   */
/*     setze pointer auf pdf und cdf  */
 
  /* make distr. object: truncated Gaussian */
  distr = unur_distr_normal(NULL,0);
  unur_distr_cont_set_domain(distr, -2.,2.);

  /* choose method and set parameters */
  par = unur_arou_new(distr);
  unur_arou_set_max_sqhratio(par,0.99);

  /* make generator object */
  gen = unur_init(par);

  /* sample */
  for (i=0;i<100;i++) {
  x = unur_sample_cont(gen);
  printf("%f\n",x);
  }
 


  /* destroy generator object */
  unur_free(gen);

  exit (0);
}






