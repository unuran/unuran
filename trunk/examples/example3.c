/* my third UNURAN program example3.c                         */

#include <unuran.h>

int main()
{
  int    i;
  double x;
  double params[2] = {10.0, 0.5};

  UNUR_DISTR *distr;    /* distribution object                */
  UNUR_PAR   *par;      /* parameter object                   */
  UNUR_GEN   *gen;      /* generator object                   */

  /* choose Gaussian distribution with the 2 parmeters
     stored in the array params  -> N(10,5)                   */
  distr = unur_distr_normal(params, 2);

  /* choose method: TABL -- an acceptance/rejection method
     using piecewise constant hats and sqeezes                */
  par = unur_tabl_new(distr);

  /* change a parameter of the used method:
     set upper bound allowed for the ratio of the areas
     below sqeeze and hat                                     */
  unur_tabl_set_max_sqhratio(par, 0.8);

  /* create generator object -- destroy parameter object      */
  gen = unur_init(par);
  if (gen == NULL){
     fprintf(stderr, "Error creating generator object\n");
     return 1;     
  }

  /* sample: print mean of 100 random numbers                 */
  for (i=0, x=0; i<100; i++) 
     x += unur_sample_cont(gen);

  printf("Mean value of 100 random numbers: %f\n",x/100);


  /* destroy distribution- and generator object               */
  unur_distr_free(distr);
  unur_free(gen);

  return 0;
}
