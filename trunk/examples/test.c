/*   my second UNURAN program test2.c */ 

#include <unuran.h>
#include <unuran_tests.h>


double trig (double x)
{  return ( (x<0.) ? 0. : (x>1.) ? 1. : x ); }  
double pdftrig (double x)
{  return ( (x<0. || x>1.) ? 0. : 1.); }
double dpdftrig (double x)
{  return 0.; }


int main()
{
  double x, x1;
  int    i;
  double q0, q1, q2, q3, q4;

  UNUR_DISTR  *distr,*distr1;    /* distribution */ 
  UNUR_PAR    *par,  *par1;      /* parameter */
  UNUR_GEN    *gen,  *gen1;      /* generator */
/* mache neues Verteilungsobjekt unur_distr_cont_new (methods.distr.c)   */
/*     setze pointer auf pdf und cdf  */
 
  /* make distr. object: truncated Gaussian */
    distr  = unur_distr_normal(NULL,0);
    distr1 = unur_distr_normal(NULL,0);
 //distr1 = unur_distr_cont_new();
 
  /* choose method and set parameters */
    par  = unur_ninv_new(distr);
    par1 = unur_ninv_new(distr1);
    unur_set_debug(par,  0);
    unur_set_debug(par1, 0);


  /* make generator object */
    gen  = unur_init(par);
    gen1 = unur_init(par1);
    if (gen1 == NULL){
      printf("NULL Pointer bei Generatorobjekt\n");
      exit(1);
    }

  /* sample */
 
  for (i=0;i<100;i++) {
  x  = unur_sample_cont(gen);
  x1 = unur_sample_cont(gen1);
  printf("%f %f \n",x, x1);
  }
 

  /* test of PP-algorithm in quantile.c */  
  q0 = 0.;
  q1 = 0.;
  q2 = 0.;
  q3 = 0.;
  q4 = 0.; 
  unur_test_quartiles( gen1, &q0 ,&q1, &q2, &q3, &q4, 1000 );

  unur_test_correlation(gen, gen1, 10000);

  /* destroy generator object */
  // unur_free(gen);
  unur_free(gen1);

  exit (0);
}















