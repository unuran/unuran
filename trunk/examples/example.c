
#include <unuran.h>
#include <unuran_distributions.h>


main ()
{
  struct unur_distr *distr;  /* distribution object */
  struct unur_par   *par;    /* parameter object */
  struct unur_gen   *gen;    /* generator object */

  double distr_params[2] = {2.,3.};  /* parameters for distribution */

  int i;
  double x;

  /* make a distribution object: 
     use own of the standard distributions */
  distr = unur_distr_beta(distr_params,2);
  if (!distr) exit(0);

  /* make parameter object for method TDR */
  par = unur_tdr_new(distr);
  if (!par) exit(0);

  /* change default parameters for method */
  unur_tdr_set_c(par,-0.5);

  /* make generator object */
  gen = unur_init(par);
  if (!gen) exit(0);

  /* sample */
  for (i=0; i<100; i++) {
    x = unur_sample_cont(gen);
    printf("x = %g\n",x);
  }

  /* destroy generator object */
  unur_free(gen);

}
