/* UNURAN example for discrete distributions             */


/* This program gives a simple example how to use
   discrete distributions with UNURAN.            */


#include <unuran.h>
int main()
{
  int    i;
  double x;
  double params[1] = {0.3};
  double probvec[10] = {1.0, 2.0, 3.0, 4.0, 5.0,\
                        6.0, 7.0, 8.0, 4.0, 3.0};

  UNUR_DISTR *distr1, *distr2;    /* distribution */
  UNUR_PAR   *par1, *par2;         /* parameter */
  UNUR_GEN   *gen1, *gen2;         /* generator */


  distr1 = unur_distr_geometric(params, 1);
  unur_distr_discr_set_domain(distr1, 0, 6);
  /* needs mode */
  unur_distr_discr_set_mode(distr1, 0);
  par1 = unur_dari_new(distr1);
  gen1 = unur_init(par1);

  distr2 = unur_distr_discr_new();
  unur_distr_discr_set_prob(distr2, probvec, 10);
  par2 = unur_dgt_new(distr2);
  gen2 = unur_init(par2);


  printf("Test\n");

  for (i=0; i<10; i++){
    unur_sample_discr(gen1);
    unur_sample_discr(gen2);
  }

  unur_distr_free(distr1);
  unur_distr_free(distr2);
  unur_free(gen1);
  unur_free(gen2);


  return 0;
}






















