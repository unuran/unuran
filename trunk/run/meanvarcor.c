/***************************************************************/
/*utility for computing mean, variance and correlation ***********/
/* this file meanvarcor.c */

/*stddev and correlation for parent population and not for sample 
  (sums divided through n not through n-1  !!!!!!!!!*/

#define SQ(a) ((a)*(a))

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct{
  double mx,my,sx,sy,sxy;
  int i,cookie;
} MEANCOR;

MEANCOR *init_meancor()
{ 
  MEANCOR *p;
  p = malloc(sizeof(MEANCOR)); 
  p->mx = 0.;
  p->sx = 0.;
  p->my = 0.;
  p->sy = 0.;
  p->sxy = 0.;
  p->i = 0;
  p->cookie = -9999;
  return p;
}

void free_meancor(MEANCOR *p)
{ free(p);
}

void update_meancor(MEANCOR *p, double xi, double yi)
{ double dx,dy;
  if(p==NULL || p->cookie != -9999){
    printf("error update_meancor; cookie check failed!! exiting\n");
    exit(1);
  }
  p->i++;
  dx= (xi - p->mx)/(double)p->i;
  p->mx+=dx;
  p->sx+= p->i*(p->i-1.)*dx*dx;
  dy= (yi - p->my)/(double)p->i;
  p->my+=dy;
  p->sy+= p->i*(p->i-1.)*dy*dy;
  p->sxy+= p->i*(p->i-1.)*dx*dy;
}

void printresult_meancor(MEANCOR *p)
{
  if(p==NULL || p->cookie != -9999){
    printf("error printresult_meancor; cookie check failed!! exiting\n");
    exit(1);
  }
  printf("correlation: %g \n",p->sxy/sqrt(p->sx*p->sy));
}

/********************************************************/
/************for computing quantiles***********/
void swap(double sample[], int i, int j)
{
 double temp;
 temp = sample[i];
 sample[i] = sample[j];
 sample[j] = temp;
}
void qksort(double sample[],int left, int right)
{
 int i, last;
 if (left >= right) return;
 swap(sample, left, (left+right)/2);
 last = left;
 for (i = left+1; i <= right; i++)
     if (sample[i] < sample[left])
         swap(sample,++last, i);
 swap(sample, left, last);
 qksort(sample, left, last-1);
 qksort(sample,last+1,right);
}


/******computing quantiles***********/
typedef struct{
  double window_max,window_min,*data;
  int i,i_window,i_below_window,cookie,n;
} QUANTILE;

QUANTILE *init_quantile(int n/*maximal size of array*/)
{ 
  QUANTILE *p;
  p = malloc(sizeof(QUANTILE)); 
  p->window_max = 1.e99;
  p->window_min = -1.e99;
  p->i = 0;
  p->i_window = 0;
  p->i_below_window = 0;
  if(n>1)
    p->n = n;
  else{
    printf("init_quantile n=%d; that's too small! NULL-pointer returned\n",n);
    return NULL;
  }
  p->data = malloc(sizeof(double)*p->n);
  p->cookie = -9191;
  return p;
}

void set_window_quantile(QUANTILE *p,double window_min,double window_max)
{ 
  if(p==NULL || p->cookie != -9191){
    printf("error set_window_quantile; cookie check failed!! window not changed\n");
  }
  else{
    if(window_min < window_max){
      p->window_min = window_min;
      p->window_max = window_max;
    }    
    else{
      printf("ERROR in set_window_quantile:\n");
      printf("window_min %f >= window_max %f window not changed!!\n",window_min,window_max);
    }
  }
}

void free_quantile(QUANTILE *p)
{ 
  if(p==NULL || p->cookie != -9191){
    printf("error free_quantile; cookie check failed!! nothing freed\n");
  }
  else{
    if(p->data) free(p->data);
    p->data = NULL;
    free(p);
  }
}

void update_quantile(QUANTILE *p, double xi)
{ 
  if(p==NULL || p->cookie != -9191){
    printf("error update_quantile; cookie check failed!! exiting\n");
    exit(1);
  }
  p->i++;
  if(p->window_min < xi && xi < p->window_max){
    if(p->i_window == p->n){
      printf("ERROR update_quantile; too many data in window, %f not stored\n",xi);
    }      
    p->data[p->i_window] = xi;
    p->i_window++;
  }
  if(p->window_min > xi)  p->i_below_window++;
}


void printresult_quantile(QUANTILE *p)
#define NQUANTILE 3
{ double q[NQUANTILE]={0.025,0.5,0.975};
 int i,nquantil;
  if(p==NULL || p->cookie != -9191){
    printf("error printresult_quantile; cookie check failed!!\n");
  }
  qksort(p->data, 0, p->i_window - 1);

  for(i=0; i<NQUANTILE; i++){
    nquantil = (int)(q[i]*p->i);
    nquantil = nquantil - p->i_below_window;
    if(nquantil < p->i_window && nquantil>= 0)
      printf("%g-p qu %g ",q[i]*100.,p->data[nquantil]);
    else printf("point for %f-percent-quantile not in window!!\n",q[i]*100.);
  }
    printf("\n");
}


void result_quantile(QUANTILE *p,int nq, double q[], double res[])
     /*nq ... number of quantiles computed
       q[nq]... array of quantile values that will be computed
       res[nq]... vector holding result
       sorts the data,computes the quantiles and puts them into res[]*/
{ 
 int i,nquantil;
  if(p==NULL || p->cookie != -9191){
    printf("error result_quantile; cookie check failed!!\n");
  }
  qksort(p->data, 0, p->i_window - 1);

  for(i=0; i<nq; i++){
    nquantil = (int)(q[i]*p->i);
    nquantil = nquantil - p->i_below_window;
    if(nquantil < p->i_window && nquantil>= 0) res[i]=p->data[nquantil];
    else {
      printf("error result_quantile: point for %f-percent-quantile not in window!! value set to 1.e100\n",q[i]*100.);
      res[i] = 1.e100;
    }
  }
}

double get_quantile(QUANTILE *p, double q) {
  double res;
  result_quantile(p, 1, &q, &res);
  return res;
}

/******computing mean and variance***********/
typedef struct{
  double mean,var,max,min;
  int i,cookie;
} MEANVAR;

MEANVAR *init_meanvar()
{ 
  MEANVAR *p;
  p = malloc(sizeof(MEANVAR)); 
  p->mean = 0.;
  p->max = -1.e99;
  p->min = 1.e99;
  p->var = 0.;
  p->i = 0;
  p->cookie = -99;
  return p;
}

void free_meanvar(MEANVAR *p)
{ free(p);
}

void update_meanvar(MEANVAR *p, double xi)
{ double dx;
  if(p==NULL || p->cookie != -99){
    printf("error update_meanvar; cookie check failed!! exiting\n");
    exit(1);
  }
  //  printf("update mit %f\n",xi);
  p->i++;
  dx= (xi - p->mean)/(double)p->i;
  p->mean+=dx;
  p->var+= p->i*(p->i-1.)*dx*dx;
  if(p->max < xi) p->max=xi;
  if(p->min > xi) p->min=xi;

 

}

void printresult_meanvar(MEANVAR *p)
{
  if(p==NULL || p->cookie != -99){
    printf("error update_meanvar; cookie check failed!! exiting\n");
    exit(1);
  }
  printf("mean: %g stddev: %g\n", p->mean, sqrt(p->var/p->i));
  printf("min: %g max: %g\n", p->min, p->max);
}

void getresult_meanvar(
           MEANVAR *p,double *mean, double *var, double *min, double *max, int *n)
/*************************
 * calculates the results and puts them into the variables pointed to;
 * NULL pointer can be handed to indicate that that result is not needed
 */
{
  if(p==NULL || p->cookie != -99){
    printf("error getresult_meanvar; cookie check failed!! exiting\n");
    exit(1);
  }
  if(mean) *mean = p->mean;
  if(var) *var = p->var/p->i;
  if(max) *max = p->max;
  if(min) *min = p->min;
  if(n) *n = p->i;
}

void getresult_meanvar_ci(
           MEANVAR *p,double *mean, double *var, double *min, double *max, 
	   int *n, double *ci/*[2]*/)
/*************************
 * calculates the results and puts them into the variables pointed to;
 * NULL pointer can be handed to indicate that that result is not needed
 * can also calculate a 95 percent CI for the mean
 */
{ double emax;
  if(p==NULL || p->cookie != -99){
    printf("error getresult_meanvar; cookie check failed!! exiting\n");
    exit(1);
  }
  if(mean) *mean = p->mean;
  if(var) *var = p->var/p->i;
  if(max) *max = p->max;
  if(min) *min = p->min;
  if(n) *n = p->i;
  if(ci){ 
    emax = 1.96*sqrt((p->var/p->i)/p->i);
    ci[0] = p->mean-emax;
    ci[1] = p->mean+emax;
    //printf("emax %f lb %f ub %f\n",emax,ci[0],ci[1]);
  } 
}



/********************************************************/

