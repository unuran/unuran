/* my fourth UNURAN program: example4.c                    */

/* Invoking the PRNG-library to generate uniformly
   distributed random numbers                              */

#include <unuran.h>

int main()
{
  int i;
  double params[2] = {1.0, 2.0};

  UNUR_DISTR *distr;    /* distribution object             */
  UNUR_PAR   *par;      /* parameter object                */
  UNUR_GEN   *gen;      /* generator object                */
  /* another UNURAN structure containig information
     about the uniform generator                           */
  UNUR_URNG  *ug;       /* uniform generator object        */


  /* Choose kind of uniform generator: Mersenne Twister
     For others see documentation of PRNG                  */
  ug =  prng_new("mt19937(1237)");

  /* choose a (standard) distribution: Beta(1,2)
     the 2 parameters are stored in params                 */
  distr = unur_distr_beta(params, 2);

  /* choos method AROU for sampling                        */
  par = unur_arou_new(distr);

  /* Set uniform generator in parameter object             */
  unur_set_urng(par, ug);

  /* make generator object                                 */
  gen = unur_init(par);

  /* destroy distribution object                           */
  unur_distr_free(distr);

  /* sample: print 100 random numbers                      */
  for (i=0; i<100; i++) 
    printf( "%f\n", unur_sample_cont(gen) );

  /* destroy generator object                              */
  unur_free(gen);

  return 0;
}
