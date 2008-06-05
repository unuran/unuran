#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <pinvwh.h>





double lobato5(double x, double h, double fx, double *fxph, double (*f)(double x))
/************************************
 * Numerical Integration of the interval (x,x+h)
 * using Gauss-Lobato integration with 5 points.
 * fx ... f(x) to save calls to f()
 * *fxph ... f(x+h) 
 ************************************/
#define W1 0.17267316464601146  //= 0.5-sqrt(3/28)
#define W2 (1.-W1)
{ double ifx,ifxph;

 ifxph = (*f)(x+h);
 if(fxph!=NULL){ 
   ifx = fx;
   *fxph = ifxph;
 }
 else ifx = (*f)(x);
 return (9*(ifx+(ifxph))+49.*((*f)(x+h*W1)+(*f)(x+h*W2))+64.*(*f)(x+h/2.))*h/180.;
}
#undef W1
#undef W2


double nint_12(double a,double b,double *res_relerror, double (*f)(double x)){
  //with 1 and two intervals reports the resulting (relative-error) in res_relerror 
  int i;
  double res,reso;

  reso = lobato5(a,(b-a),0,NULL,f);
  i=1;
  res = lobato5(a,(b-a)*0.5,0,NULL,f)+lobato5((b+a)*0.5,(b-a)*0.5,0,NULL,f);
  *res_relerror = fabs(res-reso)/res; 
  return res;
}

double nint_monoton_dens(double a,double b,double step,double crit, double (*f)(double x)){
  // step <= b-1 !!!
// numerical integration with a variable step-size for a monoton density, starting with "step"
// crit ... maximal accepted relative error

//TODO: Kann man verbessern, wenn man in den Tails den Error relativ zum Gesamtintegral rechnet

  double sum=0.,sumi,x,error;
  int i=1;

  x=a;
  while(x<b){
    sumi = nint_12(x,x+step,&error,f);
    //printf("%d: x%g step%g sumi%g error%g\n",i,x,step,sumi,error);
    if(error > crit){ 
      step = step*pow(0.5*crit/error,1./9.);
    }else{
      i++;
      sum += sumi;
      x = x+step;
      //      step = step*1.2*pow(crit/error,1./9.);
      step = step*pow(crit/error,1./9.);
      if(x+step >b){
        step = (b-x)*(1+5.e-16);
      }
    }    
  }
  return sum;
}


////////////////////////////////////////////////////////////////////////////////////

double evalnewtoninterpol(double u,int g,double ui[],double zi[]){
  /* zi = pi;  ui=yi */
  double q,chi;
  int k;
  
  q = u;
  chi = zi[g];
  for(k=g-1;k>=1;k--) chi = chi*(q-ui[k])+zi[k];
  return( chi*q);
}


int newtoninterpol(double x0, double h,int g,double ui[],double zi[],double *x,double (*f)(double x)){
  /*calculates ui and zi arrays, xi pointer may be NULL */
  double xi, dxi, temp,// zi[20], ui[20]={0.},
phi;/*20 statt g+1 */
  int i,k;
  phi = M_PI*0.5/(g+1);
  ui[0]=0.;
  if(x!=NULL)x[0]=xi;
  for(i=1; i<=g; i++){
    xi = x0 + h*sin((i-1)*phi)*sin(i*phi)/cos(phi);
    dxi = h*sin(2*i*phi)*tan(phi);
    if(x!=NULL)x[i]=xi+dxi;
   temp = lobato5(xi, dxi,0,NULL,f);
   if(temp<1.e-50){
     printf("ERROR!! Newtoninterpolation interval too short. or density 0 EXITINGNG\n");
     exit(1);
   }
    zi[i]=dxi/temp;
    ui[i]=ui[i-1]+temp; /* ui[0] initialisation??? */
  } // ui[g] is the probability of the interval
  for(k=2; k<=g; k++)
    for(i=g; i>=k; i--) zi[i]=(zi[i]-zi[i-1])/(ui[i]-ui[i-k]);

  //  printf("x:");
  //    printvec(g+1,x);
  //printf("u:");
  //  printvec(g+1,ui);

  return 1;
}


int tstpt(int g,double ui[],double utest[]){
  int k,j,i;
  double sum, qsum,x;
  utest[0]=0.;
  for(k=1; k<=g; k++){
    x = 0.5*(ui[k-1]+ui[k]);
    for(j=1; j<=2; j++){
      sum = 1./x;
      qsum = sum*sum;
      for(i=1; i<=g; i++){
	sum = sum + 1./(x-ui[i]);
	qsum = qsum + 1./((x-ui[i])*(x-ui[i]));
      }
      x +=sum/qsum;
    }
    utest[k] = x;
  }
  return 1;
}

double maxerrornewton(int g,double ui[],double zi[],double x0,double xval[],double (*f)(double x)){
  double maxerror=0.,uerror,*testu,uarr[21],x;
  int i,n;
    n=g;
    testu=uarr;
    tstpt(g,ui,testu);
    //printvec(n+1,testu);
    //printvec(n+1,ui);
    //if(xval!=NULL)printvec(3,xval);
  for(i=0;i<n;i++){
    //Naechste Zeile: TODO Verbesserung moeglich, wenn man die xi verwendet und so kuerzere
    // intervalle fuer LObato integration bekommt
    x=evalnewtoninterpol(testu[i+1],g,ui,zi);
    if(i==0||xval==NULL) uerror=fabs(lobato5(x0,x ,0.,NULL,f)-testu[i+1]);
    else uerror=fabs(ui[i]+lobato5(xval[i],x+x0-xval[i],0.,NULL,f)-testu[i+1]);
    if(uerror>maxerror) maxerror=uerror;
    //    printf("%d:u max %g %g\n",i,uerror,maxerror);
  }
  //printf("maxerror %g\n",maxerror);
  return maxerror;
}



struct genobject *init_genobject(int g,int maxint){
 struct genobject *p;
 p = malloc(sizeof(struct genobject));
 p->g = g;
 p->ni = 0;
 p->iv = malloc(sizeof(struct siv)*maxint);
 return p;
}

int free_genobject(struct genobject *p){
  int i;
  for(i=0;i<=p->ni;i++){ 
    free(p->iv[i].ui);
    free(p->iv[i].zi);
  }
  free(p->iv);
  p->iv=NULL;
  free(p->gt);
  p->gt=NULL;
  free(p);
  return 0;
}



struct genobject *setup(double (*f)(double x),int g, double a, double b, double hh, double uerror){
  /*
    f ... PDF
    g ... order of polynomial
    a ...
    b ...
    hh ...
    uerror ... u-error
  */
  double maxerror,h=hh,*xval;
  int i,j,cont,countextracalc=0,maxint=10000;
  struct genobject *geno;
  xval=malloc(sizeof(double)*(g+1));
  geno = init_genobject(g,maxint);
  geno -> iv[0].ui = malloc(sizeof(double)*(g+1));
  geno -> iv[0].zi = malloc(sizeof(double)*(g+1));
  geno -> iv[0].xi = a;
  geno -> iv[0].cdfi = 0.;//cdfi holds cdf value at the left border of the interval
  cont=1;
  i=0;
  while(cont){
    if(geno->iv[i].xi+h >b){
      h = b - geno->iv[i].xi;
      cont=0;
    }
    newtoninterpol(geno->iv[i].xi,h,g,geno->iv[i].ui,geno->iv[i].zi,xval,f);
    maxerror = maxerrornewton(g,geno->iv[i].ui,geno->iv[i].zi,geno->iv[i].xi,xval,f);
    if(maxerror > uerror){ 
      countextracalc++;
      h*= 0.9;
     if(maxerror>4.*uerror) h*=0.9;
 
    }else{
	geno->iv[i+1].ui = malloc(sizeof(double)*(g+1));
  	geno->iv[i+1].zi = malloc(sizeof(double)*(g+1));
  	geno->iv[i+1].xi = geno->iv[i].xi+h;
  	geno->iv[i+1].cdfi = geno->iv[i].cdfi +(geno -> iv)[i].ui[g];//cdfi holds cdf value at the left border of the interval
        if(maxerror < 0.3*uerror) h*=1.2;
        if(maxerror < 0.1*uerror) h*=2.;
        i++;
   }
   if(i>maxint){
     printf("error setup(); i>maxint; EXITING\n");
     exit(1);
   }
 }
 geno->ni = i;
 geno->iv = realloc(geno->iv,sizeof(struct siv)*(geno->ni+1));

 free(xval);

 geno->umax = geno->iv[geno ->ni].cdfi;


 geno->C = geno->ni;//size of guide-table
 geno->gt = malloc(sizeof(int)*geno->C);  
 geno->gt[0] = 0;
 i=0;
 for(j=1; j<geno->C; j++){
   while(j/(double)geno->C > geno->iv[i+1].cdfi/geno->umax) i++;
   /* "/geno->umax" above is necessary, as we need the guide table for u in (0,umax) */
   geno->gt[j] = i;
 }
 printf("Set-up finished: g=%d,  Number of intervals = %d,\n         additional calculated interpolations=%d\n",g,geno->ni,countextracalc);
 printf("u in (0,%.18g)   1-umax%g\n",geno->umax,1-geno->umax);
 return geno;
} 



double searchborder(double x0, double step,double border,double (*f)(double x)){
  // x0 starting point with f(x0) not small; may but need not be the mode
  // step first step size, includes direction
  double fx0=f(x0),fx=fx0,x=x0,xa;
  int i;

  for(i=0;i<100 && f(x)>fx0*1.e-13;i++){
    xa=x;
    if(i>10) step*=2;
    if((x-border)*(x+step-border)>0.) x+= step;
    else x=(x+border)*0.5;
    fx=f(x);
    //printf("%d:x%g fx%g\n",i,x,fx);
  }

  do{
    x=(x+xa)*0.5;
    fx=f(x);
    //printf("%d:x%g fx%g\n",i,x,fx);
  }while(fx<fx0*1.e-13);
 


  return x;
}

/************************************/
double tail(double x, double d,double (*f)(double x)){
  /********************************
     calculates approximate tail area = f(x)^2/((lc_f(x)+1)*abs(f'(x)) 
     x...cut off point
     d... step length for numeric differentiation
  **********************************/
  double ff,fp,fm,cplus1;
  ff = f(x);
  fp = f(x+d);
  fm = f(x-d);
  if(fm-2.*ff+fp < 0.||fm<1.e-100||ff<1.e-100||fp<1.e-100){
    printf("warning possible problem in function tail() !!!\n");
  }
  cplus1 = fp/(fp-ff) + fm/(fm-ff);
  return ff*ff/(cplus1*fabs((fp-fm)/(2.*d)));
}




double cut(double w,double dw, double crit,double (*f)(double x)){
  /**********************
    calculates starting from w the left (dw<0) or right (dw>0) cut off point.
    crit... u-error criterium for tail cut off 
    dw ... initial step-length
    the are outside is approximately crit/10
  ****************/

  double H,rezH,y,d,rezy,yplus,rezys,corr;
  int j,k,cont;

  H = crit;
  rezH=1./H;
  cont=1;
  for(j=1;j<1000&&cont;j++){
    y=tail(w,dw/64.,f);
    if(y<H) cont=0;
    else{
      w+=dw;
      if(j>32) dw*=1.5;
    }
  }
  if(cont){
    printf("error in cut(), first loop; exiting!!!\n");
    exit(1);
  }

  d=dw/64.;
  for(k=0;k<50;k++){
    rezy=1./y;
    yplus = tail(w+d,d,f);
    if(yplus<0){
      printf("error in cut(), yplus negative; exiting!!!\n");
      exit(1);
    }
    rezys = (1./yplus-rezy)/d;
    corr = -(rezy-rezH)/rezys;
    if(fabs(rezy/rezH-1.)<1.e-7) return w;
    w+=corr;
    y=tail(w,d,f);
    if(y<0){
      printf("error in cut(), y negative; exiting!!!\n");
      exit(1);
    }
  }
    printf("error in cut(), second loop; exiting!!!\n");
    exit(1);


}

/**********************************************************************/

double quantile(double u, struct genobject *geno){
  /* evaluates the inverse CDF of the generator object geno for a given
    u-value in [0,1];
    for u>1 the right border of the domain is returned.
    for u=0 the left border of the domain is returned.
    for u<0 seg-fault is possible
  */
  
int i;
 double x,un;

 if(u>=1.)u=(1-3.e-16);
 un =u*geno->umax;
 for(i= geno -> gt[(int)(u * geno ->C)]; geno->iv[i+1].cdfi < un ; i++);
 //for(i= 0; geno->iv[i+1].cdfi < un ; i++);
 //printf("i-gti%d\n",i-geno -> gt[(int)(u * geno ->C)]);

 x=evalnewtoninterpol(un- geno->iv[i].cdfi,geno->g,geno->iv[i].ui,geno->iv[i].zi);


 return (geno -> iv)[i].xi+x;

}






/******************************************************/

struct genobject *pinvsetup( double (*f)(double x),int g, double uerror, double x0, int asearch, int bsearch,double a, double b){
/***********************************************
  starts the set-up and returns a pointer to the generator object

  f      ... pointer to the pdf)
  g      ... degree of polynomial
  uerror ... maximal accepted u-error
  x0     ... value with f(x0) not small
  asearch... 1 ... search for a cut off value >= a
  bsearch... 1 ... search for a cut off value <= b
  a   ...left domain border, can be set very small eg. -1.e100 for asearch==1
  b   ...right domain border, can be set very large eg. 1.e100 for bsearch==1

  asearch==0 means that the direct a-value is left unchanged
  bsearch==0 means that the direct a-value is left unchanged
  note that this may lead to numericl problems in the set-up
  if f(x) is very small close to the border.

 So the algorithm should run more stable for  asearch and bsearch =1.
 
 asearch=0 is useful eg. for the Gamma(2) distribution where the left border 0
 is fixed and should not be searched for, as f(x) has no left tail

 ***********************************************/

  double area,areal,tailcutfact;

  if(asearch) a = searchborder(x0, -1, a,f);
  if(bsearch) b = searchborder(x0, 1, b,f);
  area = nint_monoton_dens(x0,b,1.,1.e-8,f);
  areal = nint_monoton_dens(a,x0,1.,1.e-8,f);
  area+=areal;
  printf("after searchborder: a=%g  b=%g area=%g\n",a,b,area);
//  printf("a=%g  b=%g area=%g integralerror%g\n",a,b,area,area/(cdf(b)-cdf(a))-1.);

 
  tailcutfact= 0.1;
  if(uerror<=9.e-13) tailcutfact=0.5;
/* above command necessary for Cauchy distribution where cut has problems with very small values*/
 
  if(asearch) a = cut(a,(a-b)/128,uerror*area*tailcutfact,f);
  if(bsearch) b = cut(b,(b-a)/128,uerror*area*tailcutfact,f);
  printf("after cut: a=%g b=%g\n",a,b);

  return setup(f,g,a,b,(b-a)/128,uerror*area);
}

/****************************************/
