/* ------------------------------------------------------------- */
/* File: example_vmt_0.c                                         */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

/* Sample from multinormal distribution using method VMT.        */

/* ------------------------------------------------------------- */

#define dim (3)

int main()
{
  int    i;
  double X[dim];

  /* Mean and covariance matrix for multinormal distribution.    */
  double mean[]  = { 1., 2., 3. };

  double covar[] = { 2., 2., 1., 
		     2., 4., 3.,
		     1., 3., 3. };

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  /* Create a multinormal distribution object.                   */
  distr = unur_distr_multinormal( dim, mean, covar );

  /* Choose a method: VMT.                                       */
  par = unur_vmt_new(distr);

  /* Create the generator object.                                */
  gen = unur_init(par);

  /* It is important to check if the creation of the generator   */
  /* object was successful. Otherwise `gen' is the NULL pointer  */ 
  /* and would cause a segmentation fault if used for sampling.  */
  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  /* It is possible to reuse the distribution object to create   */
  /* another generator object. If you do not need it any more,   */
  /* it should be destroyed to free memory.                      */
  unur_distr_free(distr);

  /* Now you can use the generator object `gen' to sample from   */
  /* the distribution. Eg.:                                      */
  for (i=0; i<10; i++) {
    unur_sample_vec(gen, X);
    printf("(%f,%f,%f)\n", X[0], X[1], X[2]);
  }
 
  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);

  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */
