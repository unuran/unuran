/* my second UNURAN program test2.c */ 

#include <unuran.h>
#include <unuran_distributions.h>
#include <unuran_tests.h>

double trig (double x)
{  return ( (x<0.) ? 0. : (x>1.) ? 1. : x ); }  
double pdftrig (double x)
{  return ( (x<0. || x>1.) ? 0. : 1.); }
double dpdftrig (double x)
{  return 0.; }


int main()
{
  double x;
  int    i;

  struct unur_distr *distr,*distr1;  /* distribution */ 
  struct unur_par *par,*par1;      /* parameter */
  struct unur_gen *gen,*gen1;      /* generator */
/* mache neues Verteilungsobjekt unur_distr_cont_new (methods.distr.c)   */
/*     setze pointer auf pdf und cdf  */
 
  /* make distr. object: truncated Gaussian */
  //  distr = unur_distr_normal(NULL,0);
  distr1 = unur_distr_cont_new();
  unur_distr_cont_set_cdf(distr1, trig );
  unur_distr_cont_set_pdf(distr1, pdftrig );  
  unur_distr_cont_set_dpdf(distr1,dpdftrig );  
  // unur_distr_cont_set_mode(distr1,0.5 );

  // unur_distr_cont_set_domain(distr, -2.,2.);

  /* choose method and set parameters */
  //  par = unur_arou_new(distr);
  par1 = unur_tabl_new(distr1);
  //  unur_ninv_use_newton(par1);
  // unur_arou_set_max_sqhratio(par,0.99);

  /* make generator object */
  //  gen = unur_init(par);
  gen1 = unur_init(par1);
  if (gen1 ==NULL){
    printf("NULL Pointer bei Generatorobjekt\n");
    exit(1);
  }

  /* sample */
  for (i=0;i<100;i++) {
  x = unur_sample_cont(gen1);
  printf("%f\n",x);
  }
 


  /* destroy generator object */
  // unur_free(gen);
  unur_free(gen1);

  exit (0);
}











