/* UNURAN program example_vemp.c                                   */

#include <unuran.h>

int main()
{
  int    i;
  double data[15][2] = { {1.,1.},{-1.,1.},{1.,-1.},{-1.,-1.} };       

  double result[2];

  UNUR_DISTR *distr;  /* distribution                              */
  UNUR_PAR   *par;    /* parameter                                 */
  UNUR_GEN   *gen;    /* generator                                 */

  /* create distribution object and set empirical sample           */
  distr = unur_distr_cvemp_new(2);
  unur_distr_cvemp_set_data(distr, &data[0][0], 4); 

  par = unur_vempk_new(distr);
  /* resample from the four given points                           */
  unur_vempk_set_smoothing(par, 0.);

  /* create generator object to enable sampling                    */
  gen = unur_init(par);
  /* Always perform this check                                     */
  if ( gen == NULL ){
     fprintf(stderr, "Error creating generation object\n");
     return (1);
  }

  for (i=0; i<10; i++){
    unur_sample_vec(gen, &result[0]);
    printf("(%f,%f)\n", result[0], result[1]);
  }
 
  /* destroy distribution object and generator object              */
  unur_distr_free(distr);
  unur_free(gen);

  return (0);
}





