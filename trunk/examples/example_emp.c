/* UNURAN program example_emp.c                                    */

#include <unuran.h>

int main()
{
  int    i;
  double data[15] = { -0.1,  0.05, -0.5,   0.08,  0.13,\
		      -0.21,-0.44, -0.43, -0.33, -0.3, \
		       0.18, 0.2,  -0.37, -0.29, -0.9 };

  UNUR_DISTR *distr;  /* distribution                              */
  UNUR_PAR   *par;    /* parameter                                 */
  UNUR_GEN   *gen;    /* generator                                 */

  /* create distribution object and set empirical sample           */
  distr = unur_distr_cemp_new();
  unur_distr_cemp_set_data(distr, data, 15); 

  par = unur_empk_new(distr);

  /* set set the smooting factor in paramter object                */
  unur_empk_set_smoothing(par, 0.8);
  
  /* createte generator object to enable sampling                  */
  gen = unur_init(par);
  /* Always perform this check                                     */
  if ( gen == NULL ){
     fprintf(stderr, "Error creating generation object\n");
     return (1);
  }

  /* !!! sampling from continuous empirical univariate
     distributions is also done by unur_sample_cont() !!!          */  
  for (i=0; i<10; i++)
     printf( "%f\n", unur_sample_cont(gen) );

  /* change the smoothing factor in the generator object
     during sampling                                               */
  unur_empk_chg_smoothing(gen, 0.0);
  for (i=0; i<10; i++)
     printf( "%g\n", unur_sample_cont(gen) );

  /* destroy distribution object and generator object              */
  unur_distr_free(distr);
  unur_free(gen);

  return (0);
}




