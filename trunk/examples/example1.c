/* my first UNURAN program: example1.c                      */

/* load the UNURAN library                                  */
#include <unuran.h>

int main()
{
  int    i;   /* loop variable                              */
  double x;   /* will hold the random number                */

  /* defining the three UNURAN objects                      */
  UNUR_DISTR *distr;    /* distribution object              */
  UNUR_PAR   *par;      /* parameter object                 */
  UNUR_GEN   *gen;      /* generator object                 */

  /* choose a implemented distribution: Gaussian            */
  /* 0 Parameters, therefore N(0,1)                         */
  distr = unur_distr_normal(NULL, 0);

  /* choose a method: AROU                                  */
  /* for other methods just replace "arou"  with
     the respective name (in lower case letters)            */
  par = unur_arou_new(distr);

  /* the distribution object shouldn't be changed here      */

  /* create generator object                                */
  /* destroys the parameter object                          */
  gen = unur_init(par);

  /* it is important to check if creation of
     the generator object was successfull                   */ 
  if (gen == NULL){
     fprintf(stderr, "Error creating generator object\n");
     return 1;     
  }

  /* It is possible to reuse the distribution object
     to create parameter and generator objects but
     if not, it can be destroyed.                           */     
  unur_distr_free(distr);

  /* sampling: print 100 random numbers                     */
  for (i=0; i<100; i++) {
     x = unur_sample_cont(gen);
     printf("%f\n",x);
  }

  /* destroy generator object                               */
  unur_free(gen);

  return 0;
}
