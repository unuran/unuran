/* ------------------------------------------------------------- */
/* File: example4_str.c                                          */
/* ------------------------------------------------------------- */
/* String API.                                                   */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

/* ------------------------------------------------------------- */

/* This example makes use of the PRNG library (see               */
/* http://statistik.wu-wien.ac.at/prng/) for generating          */
/* uniform random numbers.                                       */
/* To compile this example you must have set                     */
/* #define  UNUR_URNG_TYPE  UNUR_URNG_POINTER                    */
/* in `src/unuran_config.h'.                                     */

/* It also works with necessary modifications with other uniform */
/* random number generators.                                     */

/* ------------------------------------------------------------- */

int main()
{
#if UNUR_URNG_TYPE == UNUR_URNG_PRNG

  int    i;          /* loop variable                            */
  double x;          /* will hold the random number              */

  /* Declare UNURAN generator object.                            */
  UNUR_GEN   *gen;      /* generator object                      */

  /* Declare objects for uniform random number generators.       */
  UNUR_URNG  *urng1, *urng2;    /* uniform generator objects     */

  /* Create the generator object.                                */
  /* Use a predefined standard distribution:                     */
  /*   Beta with parameters 2 and 3.                             */
  /* Choose a method: TDR.                                       */
  /* Use the Mersenne Twister for unifrom random number          */
  /*   generator (requires PRNG library).                        */
  gen = unur_str2gen("beta(2,3) & method=tdr & urng = mt19937(1237)");

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

  /* Now we want to switch to a different uniform random number  */
  /* generator.                                                  */
  /* Now we use an ICG (Inversive Congruental Generator).        */
  urng2 = prng_new("icg(2147483647,1,1,0)");
  if (urng2 == NULL) exit (EXIT_FAILURE);

  /* Change uniform random number generator.                     */
  /* Notice however that we should save the pointer to uniform   */
  /* random number generator in the generator object.            */
  urng1 = unur_chg_urng( gen, urng2 );

  /* ... and sample again.                                       */
  for (i=0; i<10; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);
  }

  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);

  /* We also should destroy the uniform random number generators.*/
  prng_free(urng1);
  prng_free(urng2);

  exit (EXIT_SUCCESS);

#else
  printf("You must use the PRNG library to run this example!\n\n");
#endif

} /* end of main() */

/* ------------------------------------------------------------- */

