/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Examples                                                                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

/*  #include <stdio.h> */
/*  #include <malloc.h> */

/*  #include <unuran.h> */
/*  #include <unuran_tests.h> */

/*  #define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER) */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/*  int main() */
/*  { */
/*    UNUR_DISTR *distr; */
/*    UNUR_PAR *par; */
/*    UNUR_GEN *gen; */
/*    int i; */

/*    double fpm[] = { 2., 3. }; */

/*    unur_set_default_debug(~0u); */
/*    unur_set_stream(stdout); */

#define N 1000

/* ------------------------------------------------------------- */
/* File: example2.c                                              */
/* ------------------------------------------------------------- */

/* Include UNURAN header file.                                   */
#include <unuran.h>

static int whcdfcount=0;
/* ------------------------------------------------------------- */

/* In this example we build a distribution object from scratch   */
/* and sample from this distribution.                            */
/*                                                               */
/* We use method TDR (Transformed Density Rejection) which       */
/* required a PDF and the derivative of the PDF.                 */

/* ------------------------------------------------------------- */


/* The CDF of our distribution:                                  */
double mycdf( double x, UNUR_DISTR *distr )
     /* The second argument (`distr') can be used for parameters */
     /* for the PDF. (We do not use parameters in our example.)  */
{ 
  whcdfcount++;
  return 1.-exp(-x);


} /* end of mycdf() */


/* ------------------------------------------------------------- */

int main()
{
  int    i;     /* loop variable                                 */
  double x;     /* will hold the random number                   */

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */


  unur_set_default_debug(~0u); 

  /* Create a new distribution object from scratch.              */
  /* It is a continuous distribution, and we need a PDF and the  */
  /* derivative of the PDF. Moreover we set the domain.          */

  /* Get empty distribution object for a continuous distribution */
  distr = unur_distr_cont_new();

  /* Assign the PDF and dPDF (defined above).                    */
  unur_distr_cont_set_cdf( distr, mycdf );
  /*  unur_distr_cont_set_dpdf( distr, mydpdf );*/

  /* Set the domain of the distribution (optional for TDR).      */
  unur_distr_cont_set_domain( distr, 0, 1.e20 );

  /* Choose a method: TDR.                                       */
  par = unur_ninv_new(distr);
    unur_ninv_set_start(par,1.,10.);
/*wenn man das set_start hier weglaesst funktioniert alles
bestens*/

  gen = unur_init(par);

  if (gen == NULL) {
     fprintf(stderr, "ERROR: cannot create generator object\n");
     exit (EXIT_FAILURE);
  }

  unur_distr_free(distr);

  printf("for set-up %d evaluations of the cdf \n",whcdfcount);
  whcdfcount=0;
  
  /* Now you can use the generator object `gen' to sample from   */
  /* the distribution. Eg.:                                      */
  for (i=0; i<N; i++) {
    x = unur_sample_cont(gen);
    printf("%f\n",x);

  }


  /* When you do not need the generator object any more, you     */
  /* can destroy it.                                             */
  unur_free(gen);



  printf("whcdfcount %f \n",(double)whcdfcount/N);


  exit (EXIT_SUCCESS);

} /* end of main() */

/* ------------------------------------------------------------- */


/*    distr = unur_distr_normal(NULL,0); */
/*    par = unur_tdr_new( distr ); */
/*    unur_run_tests(par,UNUR_TEST_TIME); */
/*    unur_distr_free(distr); */

/*---------------------------------------------------------------------------*/














