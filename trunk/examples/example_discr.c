/* UNURAN example for discrete distributions                       */

/* This program gives a simple example how to use
   discrete distributions with UNURAN.                             */

#include <unuran.h>

int main()
{
  int    i;
  double param = 0.3;
  double probvec[10] = {1.0, 2.0, 3.0, 4.0, 5.0,\
                        6.0, 7.0, 8.0, 4.0, 3.0};
  
  UNUR_DISTR *distr1, *distr2;    /* distribution                  */
  UNUR_PAR   *par1, *par2;        /* parameter                     */
  UNUR_GEN   *gen1, *gen2;        /* generator                     */

  /* First generator will sample from a truncated
     geometric distribution (1 parameter) using the
     method DARI (which needs the mode)                            */
  distr1 = unur_distr_geometric(&param, 1);
  unur_distr_discr_set_mode(distr1, 0);
  par1 = unur_dari_new(distr1);
  gen1 = unur_init(par1);
  /* Always perform this check                                     */
  if ( gen1 == NULL ){
     fprintf(stderr, "Error creating first generation object\n");
     return (1);
  }

  /* Second generator will sample from a probability
     vector using the method DAU -- the urnfactor
     ist set to 1.5 (standard value would be 1.0                   */
  distr2 = unur_distr_discr_new();
  /* probability vectors need not to sum up to 1                   */
  unur_distr_discr_set_pv(distr2, probvec, 10);
  par2 = unur_dau_new(distr2);
  unur_dau_set_urnfactor(par2, 1.5);
  gen2 = unur_init(par2);
  /* Always perform this check                                     */
  if ( gen2 == NULL ){
     fprintf(stderr, "Error creating second generation object\n");
     return (1);
  }
  
  /* print some random integers                                    */
  for (i=0; i<10; i++){
    printf("number %d: %d\n", i*2,   unur_sample_discr(gen1) );
    printf("number %d: %d\n", i*2+1, unur_sample_discr(gen2) );
  }

  /* clear allocated memory                                        */
  unur_distr_free(distr1);
  unur_distr_free(distr2);
  unur_free(gen1);
  unur_free(gen2);


  return 0;
}

