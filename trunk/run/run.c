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

#include <stdio.h>
#include <malloc.h>

#include <unuran.h>
#include <unuran_tests.h>

#define RUN_TESTS       (~0x0u & ~UNUR_TEST_SCATTER)

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

/* ------------------------------------------------------------- */

int main()
{
  int    i;     /* loop variable                                 */
  double x;     /* will hold the random number                   */
  double cdfbmode;


  double moments[10];

  /* Declare the three UNURAN objects.                           */
  UNUR_DISTR *distr;    /* distribution object                   */
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  double fpar[] = {50};

  unur_set_default_debug(~0u);

  distr = unur_distr_normal(NULL,0);
  par = unur_gsrou_new(distr);
/*    gen = unur_init(par); */
  
  unur_run_tests(par,RUN_TESTS);

/*    gen = unur_init(par); */

/*    for (i=0; i<100; i++) */
/*      unur_sample_discr(gen); */
/*      printf("%d\n",unur_sample_discr(gen)); */

/*    unur_test_count_urn(gen,10000,1,stdout); */

/*    unur_test_moments( gen, moments, 2, 24000, 0, stdout ); */
/*    printf("mean = %g\n",moments[1]); */

  return 0;

}

/*---------------------------------------------------------------------------*/













