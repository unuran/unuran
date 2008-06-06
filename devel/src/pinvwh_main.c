

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <unuran.h>
#include <pinv.h>

int check_inversion_unuran(UNUR_GEN *gen,double uerror,double (*cdf)(double x),int printyn);


/*auch getestet mit f eine quasi-dichte mit Integral nicht 1!! Keine Probleme */

/* double f_normal(double x){//PDF normal distribution */
/* return exp(-x*x*0.5)/sqrt(2*M_PI); */
/* } */
/*die faktoren in den Dichten sind nur, um zu testen, dass die Methode auch mit
  nicht normierten Quasi-dichten funktioniert.
*/

double cdf_normal(double x)//CDF normal distribution
{ return 0.5*(1+erf(x*.7071067811865475244));
}



/* double f_gamma2(double x){//PDF gamma 2 distribution */
/* return x*exp(-x); */
/* } */

double cdf_gamma2(double x)//CDF gamma 2 distribution
{ return 1.-exp(-x)*(1.+x);
}




/* double f_cauchy(double x){//PDF half Cauchy distribution */
/* return 1./(M_PI*(1.+x*x)); */
/* } */


double cdf_cauchy(double x)//CDF half Cauchy distribution
{ return 0.5+atan(x)/M_PI;
}



////////////////////////////////////////////////////////////////////////////////////


int main(){


  UNUR_DISTR *normal, *cauchy, *gamma2;
  double fpm[1] = {2.};

  UNUR_PAR *par;
  UNUR_GEN *gen;

  double uerror=1.e-11,x0=1.;
  int g;  
 
  normal = unur_distr_normal(NULL,0);
  unur_distr_cont_set_center(normal,x0);

  cauchy = unur_distr_cauchy(NULL,0);
  unur_distr_cont_set_center(cauchy,x0);

  gamma2 = unur_distr_gamma(fpm,1);
  unur_distr_cont_set_center(gamma2,x0);



  for(g=3;g<12;g++){
    for(uerror=1.e-8;uerror>2.e-13;uerror*=0.1){
      
      /* --------------------------------------------- */
      printf("\n Normal distribution\n");

      par = unur_pinv_new(normal);
      unur_pinv_set_order(par,g);
      unur_pinv_set_u_resolution(par,uerror);
      unur_pinv_set_boundary(par,-1.e100, 1.e100);
      unur_pinv_set_searchboundary(par,TRUE,TRUE);
      gen = unur_init(par);
      check_inversion_unuran(gen,uerror,cdf_normal,0);
      unur_free(gen);
      
      
      /* --------------------------------------------- */
      printf("\n Gamma(2) distribution\n");
      
      par = unur_pinv_new(gamma2);
      unur_pinv_set_order(par,g);
      unur_pinv_set_u_resolution(par,uerror);
      unur_pinv_set_boundary(par,0, 1.e100);
      unur_pinv_set_searchboundary(par,FALSE,TRUE);
      gen = unur_init(par);
      check_inversion_unuran(gen,uerror,cdf_gamma2,0);
      unur_free(gen);
      
      
      /* --------------------------------------------- */
      printf("\n Cauchy distribution\n");

      par = unur_pinv_new(cauchy);
      unur_pinv_set_order(par,g);
      unur_pinv_set_u_resolution(par,uerror);
      unur_pinv_set_boundary(par,-1.e100, 1.e100);
      unur_pinv_set_searchboundary(par,TRUE,TRUE);
      gen = unur_init(par);
      check_inversion_unuran(gen,uerror,cdf_cauchy,0);
      unur_free(gen);
      
    }
  }

  unur_distr_free(normal);
  unur_distr_free(cauchy);
  unur_distr_free(gamma2);

  return 0;
}
