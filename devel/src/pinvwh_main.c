

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <pinvwh.h>
#include <unuran.h>




/*auch getestet mit f eine quasi-dichte mit Integral nicht 1!! Keine Probleme */

double f_normal(double x){//PDF normal distribution
return exp(-x*x*0.5)/sqrt(2*M_PI);
}
/*die faktoren in den Dichten sind nur, um zu testen, dass die Methode auch mit
  nicht normierten Quasi-dichten funktioniert.
*/

double cdf_normal(double x)//CDF normal distribution
{ return 0.5*(1+erf(x*.7071067811865475244));
}



double f_gamma2(double x){//PDF gamma 2 distribution
return x*exp(-x);
}

double cdf_gamma2(double x)//CDF gamma 2 distribution
{ return 1.-exp(-x)*(1.+x);
}




double f_cauchy(double x){//PDF half Cauchy distribution
return 1./(M_PI*(1.+x*x));
}


double cdf_cauchy(double x)//CDF half Cauchy distribution
{ return 0.5+atan(x)/M_PI;
}



////////////////////////////////////////////////////////////////////////////////////

void printvec(int n,double v[]){
  int i;
  printf("(");
  for(i=0;i<n;i++) printf(" %g ,",v[i]);
  printf(")\n");
}



/****************************************/
int check_inversion(struct genobject *geno,double uerror,double (*cdf)(double x),int printyn){
  /* checks the inversion for a generator object by calculating the u-error using the CDF
     checks many u-values especially in the etreme tails
  */

  double u,x,uerr,maxerror=0.,maxu;

  uerror*=1.1;
/* "*1.1" is necessary for this controll as the left cut-off tail has about 10% of the uerror */
  maxu=1.;
  for(u=1.e-10;u<maxu;u+=1./33211){
    x=quantile(u,geno);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror) printf("%g %g %g\n",u,x,uerr);
  }
  for(u=1.e-13;u<1.e-5;u+=1.e-8){//check left tail
    x=quantile(u,geno);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror)printf("%g %g %g\n",u,x,uerr);
  }
  for(u=maxu-1.e-13;u>maxu-1.e-5;u-=1.e-8){//check right tail
    x=quantile(u,geno);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror)printf("%g %g %g\n",u,x,uerr);
  }
  printf("aimed uerror: %g ;  observed maxerror : %g\n",uerror,maxerror);
  return 0;
}

/**********************************************/





int main(){


  UNUR_DISTR *normal, *cauchy, *gamma2;
  double fpm[1] = {2.};

  double uerror=1.e-11,x0=1.;
  int g;  
  struct genobject *geno;
  
 
  normal = unur_distr_normal(NULL,0);
  cauchy = unur_distr_cauchy(NULL,0);
  gamma2 = unur_distr_gamma(fpm,1);


 for(g=3;g<12;g++){
   for(uerror=1.e-8;uerror>2.e-13;uerror*=0.1){
     printf("\n Normal distribution\n");
     geno = pinvsetup(f_normal,g, uerror, x0, 1, 1,-1.e100, 1.e100);
     check_inversion(geno,uerror,cdf_normal,0);
     free_genobject(geno);



     printf("\n Gamma(2) distribution\n");
     geno = pinvsetup(f_gamma2,g, uerror, x0, 0, 1,0, 1.e100);//fuer Gamma 2
     check_inversion(geno,uerror,cdf_gamma2,0);
     free_genobject(geno);



     printf("\n Cauchy distribution\n");
     geno = pinvsetup(f_cauchy,g, uerror, x0, 1, 1,-1.e100, 1.e100);
     check_inversion(geno,uerror,cdf_cauchy,0);
     free_genobject(geno);
   }
 }
 return 0;
}
