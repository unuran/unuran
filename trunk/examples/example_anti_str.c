/* ------------------------------------------------------------- */
/* File: example_anti_str.c                                      */
/* ------------------------------------------------------------- */
/* String API.                                                   */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

/* Example how to sample from two streams of antithetic random   */
/* variates from Gaussian N(2,5) and Gamma(4) distribution, resp.*/

/* ------------------------------------------------------------- */

/* This example makes use of the PRNG library (see               */
/* http://statistik.wu-wien.ac.at/prng/) for generating          */
/* uniform random numbers.                                       */
/* To compile this example you must have set                     */
/* #define  UNUR_URNG_TYPE  UNUR_URNG_PRNG                       */
/* in `src/unuran_config.h'.                                     */

/* It also works with necessary modifications with other uniform */
/* random number generators.                                     */

/* ------------------------------------------------------------- */

int main()
{
#if UNUR_URNG_TYPE == UNUR_URNG_PRNG

  int    i;          /* loop variable                            */
  double xn, xg;     /* will hold the random number              */

  /* Declare UNURAN generator objects.                           */
  UNUR_GEN   *gen_normal, *gen_gamma;

  /* PRNG only:                                                  */
  /* Make a object for uniform random number generator.          */
  /* For details see http://statistik.wu-wien.ac.at/prng/.       */

  /* Create the first generator: Gaussian N(2,5)                 */
  gen_normal = unur_str2gen("normal(2,5) & method=tdr; variant_ps & urng=mt19937(1237)");
  if (gen_normal == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }
  /* Set auxilliary uniform random number generator.             */
  /* We use the default generator.                               */
  unur_chgto_urng_aux_default(gen_normal);

  /* The second generator: Gamma(4) with antithetic variates.    */
  gen_gamma = unur_str2gen("gamma(4) & method=tdr; variant_ps & urng=anti(mt19937(1237))");
  if (gen_gamma == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }
  unur_chgto_urng_aux_default(gen_gamma);


  /* Now we can sample pairs of negatively correlated random     */
  /* variates. E.g.:                                             */
  for (i=0; i<10; i++) {
    xn = unur_sample_cont(gen_normal);
    xg = unur_sample_cont(gen_gamma);
    printf("%g, %g\n",xn,xg);
  }


  /* When you do not need the generator objects any more, you    */
  /* can destroy it.                                             */

  /* But first we have to destroy the uniform random number      */
  /* generators.                                                 */
  prng_free(unur_get_urng(gen_normal));
  prng_free(unur_get_urng(gen_gamma));

  unur_free(gen_normal);
  unur_free(gen_gamma);


  exit (EXIT_SUCCESS);

#else
  printf("You must use the PRNG library to run this example!\n\n");
  exit (EXIT_FAILURE);
#endif

} /* end of main() */

/* ------------------------------------------------------------- */

