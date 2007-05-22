/* ------------------------------------------------------------- */
/* File: example.c                                               */
/* ------------------------------------------------------------- */

/* Include UNURAN header files.                                  */
#include <unuran.h>
#include <unuran_urng_rngstreams.h>

/* ------------------------------------------------------------- */

int main(void)
{
  int    i;          /* loop variable                            */
  double x;          /* will hold the random number              */

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  /* Declare objects for uniform random number generators.       */
  UNUR_URNG  *urng_global;    /* uniform RN generator object     */
  UNUR_URNG  *urng_local;     /* uniform RN generator object     */

  /* -- Prepare for using RNGSTREAMS library ------------------- */

  /* The RNGSTREAMS library sets a package seed.                 */
  unsigned long seed[] = {111u, 222u, 333u, 444u, 555u, 666u};
  RngStream_SetPackageSeed(seed);

  /* Make a object for uniform random number generator.          */
  /* For details see                                             */
  /* http://statmath.wu-wien.ac.at/software/RngStreams/          */
  urng_global = unur_urng_rngstream_new("urng-global");
  if (urng_global == NULL) exit (EXIT_FAILURE);

  /* Set default URNG.                                           */
  unur_set_default_urng( urng_default );

  /* -- Example 1 ---------------------------------------------- */

  /* Create generator for the Beta distribution with shape       */
  /* parameters 2 and 3. Use method "AUTOmatic".                 */
  gen = unur_str2gen("beta(2,3)");

  /* It is important to check if the creation of the generator   */
  /* object was successful. Otherwise `gen' is the NULL pointer  */ 
  /* and would cause a segmentation fault if used for sampling.  */
  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  /* Now you can use the generator object `gen' to sample from   */
  /* the distribution. Eg.:                                      */
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }

  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);

  /* -- Example 2 ---------------------------------------------- */

  /* Create generator for the t distribution with 3 degrees of   */
  /* freedom. Use method "TDR" (Transformed Density Rejection).  */
  gen = unur_str2gen("student(3) & method=TDR");
  
  /* Now we want to switch to a different (independent) stream   */
  /* of uniform random numbers.                                  */
  urng_local = unur_urng_rngstream_new("urng-local");
  if (urng2 == NULL) exit (EXIT_FAILURE);
  unur_chg_urng( gen, urng-local );

  /* Draw sample ...                                             */
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }

  /* Destroy generator object (when it is not required any more) */
  unur_free(gen);

  /* We also should destroy the uniform random number generators.*/
  unur_urng_free(urng-local);
  unur_urng_free(NULL);

  exit (EXIT_SUCCESS);
} /* end of main() */

/* ------------------------------------------------------------- */
