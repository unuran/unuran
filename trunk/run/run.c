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

#define RUN_TESTS       (~0x0u)

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int main()
{
  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  double *data;

  unur_set_default_debug(~0u);

  data = malloc(10000*sizeof(double));
  gen = unur_str2gen("normal & method=cstd");
  for (i=0;i<10000;i++)
    data[i] = unur_sample_cont(gen);
  unur_free(gen);

  unur_set_default_debug(~0u);
  distr = unur_distr_cemp_new();
  unur_distr_cemp_set_data(distr,data,10000);
  free(data);

  par = unur_empk_new(distr);
  
  unur_run_tests(par,RUN_TESTS);

  return 0;

}

/*---------------------------------------------------------------------------*/













