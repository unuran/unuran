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
#include <stdlib.h>
#include <unistd.h>

#include <unuran.h>
#include <unur_struct.h>
#include <unur_source.h>
#include <unuran_tests.h>

#include <src/utils/matrix_source.h>
#include <experimental/hitrou.h>
#include <experimental/gibbs.h>
#include "meanvarcor.c"

#define DISTRIBUTION_NORMAL 0 
#define DISTRIBUTION_STUDENT 1

#define COVAR_CONSTANT 0
#define COVAR_NEIGHBOURS 1
#define COVAR_POWER 2

#define METHOD_HITROU 0
#define METHOD_VMT 1
#define METHOD_GIBBS 2
#define METHOD_GIBBS_RANDOM 3
#define METHOD_HITROU_BOX_COORDINATE 4
#define METHOD_HITROU_BOX 5
#define METHOD_HITROU_STRIP_ADAPTIVE 6
#define METHOD_HITROU_BALL 7

#define MAXDIM 100

#define VERBOSE  0     /* output switch for test and moments  */

/* default values */
int DIM = 2;
long SAMPLESIZE = 10000;
long EXPERIMENTS = 100;
int DISTRIBUTION = 0;
int NU = 1;
int COVAR = 2;
double RHO = 0.0;
double SIGMA = 1.0;
int METHOD = 0; 
long SKIP = 0;
double BALL_RADIUS = 0.1;

/*---------------------------------------------------------------------------*/

void printarray(char *info, double *a) {
#define im(a,b)  ((a)*5+(b))    /* index for moments */
  int d, m;
  double average;
  
  printf("%s", info);

  for (m=1; m<=4; m++) {
    average=0;
    for (d=0; d<DIM; d++) {
      average += a[im(d,m)];
    }
    average /= DIM;
    printf("  % .3e", average);
    if (isinf(average)) printf("     ");
  }
  
  printf("\n");

#undef im
}


void math2(struct unur_gen *gen) {

  int i,j;
  double x[100];
  FILE *fn;
  
  fn = fopen("vertex.txt","w");
  
    fprintf(fn,"vertex={");
    for (i=1; i<=SAMPLESIZE; i++) {
      unur_sample_vec(gen, x);  
      
      fprintf(fn,"{");
      for (j=0; j<DIM; j++) {
        fprintf(fn,"%4.2f", x[j]);
	if (j<DIM-1) fprintf(fn,",");
      }
      fprintf(fn,"}");
      
      //printf("{%4.2f,%4.2f}", x[0], x[1]);
      if (i<SAMPLESIZE ) fprintf(fn,",");
    }
    fprintf(fn,"};\n");

    fclose(fn);
}




int main(int argc, char *argv[])
{
#define ic(a,b)  ((a)*DIM+(b))  /* index for covariance matrix */
#define im(a,b)  ((a)*5+(b))    /* index for moments */

  UNUR_DISTR *distr=NULL;
  UNUR_PAR *par=NULL;
  UNUR_GEN *gen=NULL;

  double *x;
  double *mean;
  double *covar;
  
  double *moments;           /* calculated moments */
  double *moments_expected;  /* expected moments   */
  
  double q, diff;
  char info[20];
  
  double time_setup, time_sample;
  
  int i, j, d, m;
  long loop;
  double s, si, sj;
    
  double a_mean[MAXDIM*5];
  double a_stddev[MAXDIM*5];
  double a_var[MAXDIM*5];
  double a_mse[MAXDIM*5];
  double a_min[MAXDIM*5];
  double a_max[MAXDIM*5];
  double a_q1[MAXDIM*5];
  double a_q2[MAXDIM*5];
  double a_q3[MAXDIM*5];
  
  MEANVAR *mv[MAXDIM*5];
  QUANTILE *quant[MAXDIM*5];  

  char c;

  /* read options */
  while ((c = getopt(argc, argv, "d:s:n:r:w:t:m:c:e:f:b:h")) != -1) {
    switch (c) {
    case 't':     /* type of distribution  */
      DISTRIBUTION=atol(optarg);
      break;
    case 'm':     /* sampling method  */
      METHOD=atol(optarg);
      break;
    case 'd':     /* dimension */
      DIM=atoi(optarg);
      break;
    case 'f':     /* degrees of freedom for student distribution */
      NU=atoi(optarg);
      break;
    case 'n':     /* number of samples */
      SAMPLESIZE=(long) atof(optarg);
      break;
    case 'e':     /* number of experiments */
      EXPERIMENTS=(long) atof(optarg);
      break;
    case 's':     /* skip size */
      SKIP=(long) atof(optarg);
      break;
    case 'c':     /* degrees of freedom for student distribution */
      COVAR=atoi(optarg);
      break;
    case 'r':     /* covariance rho */
      RHO=atof(optarg);
      break;
    case 'w':     /* covariance rho */
      SIGMA=atof(optarg);
      break;
    case 'b':     /* ball radius */
      BALL_RADIUS=atof(optarg);
      break;
    
    case 'h':     /* help */
      printf("options\n" );
      printf(" -d dim          : dimension (%d) \n", DIM );
      printf(" -t type         : 0=normal, 1=student (%d) \n", DISTRIBUTION );
      printf(" -f nu           : degrees of freedom for student (%d) \n", NU );
      printf(" -m method       : 0=H&R+RD+STRIP, 1=VMT, 2=GIBBS 3=GIBBS+RD \n" );
      printf("                 : 4=H&R+COORD+BOX, 5=H&R+RD+BOX 6=H&R+RD+ADAPTIVE STRIP\n" );
      printf("                 : 7=H&R+BALL (%d)\n", METHOD);
      printf(" -b ball_radius  : ball radius for ball sampler (%f)\n", BALL_RADIUS);
      printf(" -s skip         : skip parameter for HITROU (%ld) \n", SKIP );
      printf(" -c covar_matrix : 0=constant, 1=neighbours, 2=power (%d)\n", COVAR);
      printf(" -r rho          : covariance parameter (%f)\n", RHO);
      printf(" -w sigma        : variance parameter (%f)\n", SIGMA);
      printf(" -n samplesize   : number of samples (%ld) \n", SAMPLESIZE );
      printf(" -e experiments  : number of repetitions (%ld) \n", EXPERIMENTS );
      
      return 0;
      break;
    default:
      break;
    }
  }
  
  unur_set_default_debug(UNUR_DEBUG_OFF);
   
  for(d=0;d<DIM;d++){
    for (m=1; m<=4; m++) {
       mv[im(d,m)]=init_meanvar();
       if (EXPERIMENTS>1) quant[im(d,m)] = init_quantile(EXPERIMENTS);
    }
  }    
  
  x = _unur_vector_new(DIM);
  mean = _unur_vector_new(DIM);
  covar = _unur_vector_new(DIM*DIM);

  moments = _unur_vector_new(DIM*5); 
  moments_expected = _unur_vector_new(DIM*5);
      
  for (i=0; i<DIM; i++) {
  for (j=0; j<=i ; j++) {
    if (COVAR==COVAR_CONSTANT) covar[ic(i,j)] = (i==j) ? 1. : RHO;
    if (COVAR==COVAR_NEIGHBOURS) covar[ic(i,j)] = (i==j) ? 1. :  ((i-j)==1) ? RHO: 0;
    if (COVAR==COVAR_POWER) covar[ic(i,j)] = pow(RHO,(i-j));
    covar[ic(j,i)] = covar[ic(i,j)];
  }}
 
  for (i=0; i<DIM; i++) {
  for (j=0; j<DIM; j++) {      
    si=(i<DIM/2.) ? 1: sqrt(SIGMA);
    sj=(j<DIM/2.) ? 1: sqrt(SIGMA);
    covar[ic(i,j)] *= si * sj;
  }}
  
  
#if 0 
  _unur_matrix_print_matrix ( DIM, covar, "Covariance Matrix",
			      stdout, "", "---" );
#endif
			        
  if (DISTRIBUTION==DISTRIBUTION_NORMAL) {
    printf("DISTRIBUTION='MULTINORMAL'\n");
    distr = unur_distr_multinormal(DIM,mean,covar);
    for (d=0; d<DIM; d++) {
      moments_expected[im(d,1)] = 0.;
      moments_expected[im(d,2)] = 1.;
      moments_expected[im(d,3)] = 0.;
      moments_expected[im(d,4)] = 3.;
    }
  }
  
  if (DISTRIBUTION==DISTRIBUTION_STUDENT) {
    printf("DISTRIBUTION='MULTISTUDENT'");
    if (NU==1) printf(" (MULTICAUCHY)"); printf("\n");
    distr = unur_distr_multistudent(DIM,NU,mean,covar);
    for (d=0; d<DIM; d++) {
      moments_expected[im(d,1)] = INFINITY;
      moments_expected[im(d,2)] = INFINITY;
      moments_expected[im(d,3)] = INFINITY;
      moments_expected[im(d,4)] = INFINITY;
    
      if (NU>=2) moments_expected[im(d,1)]=0;
      if (NU>=3) moments_expected[im(d,2)]=NU/(NU-2.);
      if (NU>=4) moments_expected[im(d,3)]=0;
      if (NU>=5) moments_expected[im(d,4)]=3. * NU*NU / ((NU-2.)*(NU-4.));
    }
    printf("NU=%d\n", NU);
  }

  printf("DIM=%d\n", DIM);
  if (COVAR==COVAR_CONSTANT)   printf("COVAR=CONSTANT\n");
  if (COVAR==COVAR_NEIGHBOURS) printf("COVAR=NEIGHBOURS\n");
  if (COVAR==COVAR_POWER)      printf("COVAR=POWER\n");
  printf("RHO=%f\n", RHO);
  printf("SIGMA=%f\n", SIGMA);
  printf("SAMPLESIZE=%ld\n", SAMPLESIZE);
  printf("EXPERIMENTS=%ld\n", EXPERIMENTS);
  
  if (METHOD==METHOD_HITROU) {
    printf("METHOD=HITROU (STRIP)\n");
    par = unur_hitrou_new(distr);
    unur_hitrou_set_variant_random_direction(par);
    unur_hitrou_set_skip(par,SKIP);
  }
    
  if (METHOD==METHOD_VMT) {
    printf("METHOD=VMT\n");  
    par = unur_vmt_new(distr);
  }
  
  if (METHOD==METHOD_GIBBS) {
    printf("METHOD=GIBBS \n");  
    par = unur_gibbs_new(distr);
    unur_gibbs_set_skip(par,SKIP);
  }
  
  if (METHOD==METHOD_GIBBS_RANDOM) {
    printf("METHOD=GIBBS (RANDOM DIRECTIONS)\n");  
    par = unur_gibbs_new(distr);
    unur_gibbs_set_skip(par,SKIP);
    unur_gibbs_set_variant_random_direction(par);
  }

  if (METHOD==METHOD_HITROU_BOX) {
    printf("METHOD=HITROU (BOX)\n");
    par = unur_hitrou_new(distr);
    unur_hitrou_use_bounding_rectangle(par,1);
    unur_hitrou_set_variant_random_direction(par);
    unur_hitrou_set_skip(par,SKIP);
  }
  
  if (METHOD==METHOD_HITROU_BOX_COORDINATE) {
    printf("METHOD=HITROU (BOX + COORDINATE SAMPLER)\n");
    par = unur_hitrou_new(distr);
    unur_hitrou_use_bounding_rectangle(par,1);
    unur_hitrou_set_variant_coordinate(par);
    unur_hitrou_set_skip(par,SKIP);
  }
   
  if (METHOD==METHOD_HITROU_STRIP_ADAPTIVE) {
    printf("METHOD=HITROU (ADAPTIVE STRIP)\n");
    par = unur_hitrou_new(distr);
    unur_hitrou_set_variant_random_direction(par);
    unur_hitrou_set_adaptive_strip(par, 1);    
    //unur_hitrou_set_v(par, 0.1);
    unur_hitrou_set_skip(par,SKIP);
  }  

  if (METHOD==METHOD_HITROU_BALL) {
    printf("METHOD=HITROU (BALL)\n");
    printf("RADIUS=%f\n", BALL_RADIUS);
    par = unur_hitrou_new(distr);
    unur_hitrou_set_variant_ball(par);
    unur_hitrou_set_ball_radius(par, BALL_RADIUS);
    unur_hitrou_use_bounding_rectangle(par, 1); /* TEST : to compute vmax/2 for initial point */
    unur_hitrou_set_skip(par,SKIP);
  }  
  
       
  printf("SKIP=%ld\n", SKIP);

#if 0
  gen = unur_init(par);
  math2(gen);
  exit(0);
#else    
  gen = unur_test_timing(par, 3, &time_setup, &time_sample, TRUE, stdout);
#endif
    
  /* main loop */    
  for (loop=1; loop<=EXPERIMENTS; loop++) {
    
    unur_test_moments(gen, moments, 4, SAMPLESIZE, VERBOSE, stdout);
 
    for(d=0;d<DIM;d++){
      s=(d<DIM/2.) ? 1: sqrt(SIGMA);
      
      for (m=1; m<=4; m++) {
        update_meanvar(mv[im(d,m)],moments[im(d,m)]/pow(s,m));
        if (EXPERIMENTS>1) update_quantile(quant[im(d,m)],moments[im(d,m)]/pow(s,m));
      }
    }    
    
  }

  /* output of results */   
  printf("------------------------------------------------------------------\n");
  printf("                #1          #2          #3          #4\n");
  printf("expect   : ");
  for (m=1; m<=4; m++) {
    printf(" % .3e ", moments_expected[m]);
  }      
  printf("\n");
  
  for (m=1; m<=4; m++) {
    for(d=0;d<DIM;d++) {
      a_mean[im(d,m)] = mv[im(d,m)]->mean;  
      a_stddev[im(d,m)] = sqrt(mv[im(d,m)]->var/mv[im(d,m)]->i); 
      a_var[im(d,m)] = mv[im(d,m)]->var/mv[im(d,m)]->i;       
      diff=mv[im(d,m)]->mean-moments_expected[im(d,m)];
      a_mse[im(d,m)] = mv[im(d,m)]->var/mv[im(d,m)]->i + diff*diff;      
      a_min[im(d,m)] = mv[im(d,m)]->min;
      a_max[im(d,m)] = mv[im(d,m)]->max;
      if (EXPERIMENTS>1) {
        q=2.5 ; a_q1[im(d,m)] = get_quantile( quant[im(d,m)], q/100); 
        q=50.0; a_q2[im(d,m)] = get_quantile( quant[im(d,m)], q/100); 
        q=97.5; a_q3[im(d,m)] = get_quantile( quant[im(d,m)], q/100); 
      }
  }}
  
  printarray("mean     :", a_mean);
  printarray("stddev   :", a_stddev);
  printarray("var      :", a_var);
  printarray("mse      :", a_mse);
  printf("\n");
  printarray("min      :", a_min);
  if (EXPERIMENTS>1) {
    q= 2.5; sprintf(info, "q.%03d    :", (int)(q*10));  
    printarray(info, a_q1);
    q=50.0; sprintf(info, "q.%03d    :", (int)(q*10));  
    printarray(info, a_q2);
    q=97.5; sprintf(info, "q.%03d    :", (int)(q*10));  
    printarray(info, a_q3);
  }
  printarray("max      :", a_max);
  
  for(d=0;d<DIM;d++){
    for (m=1; m<=4; m++) {
      free_meanvar(mv[im(d,m)]);
      if (EXPERIMENTS>1) free_quantile(quant[im(d,m)]);
    }
  }
    
  unur_distr_free(distr);
  unur_free(gen);
  
  free(x); free(mean); free(covar);
  free(moments); free(moments_expected); 
    
  return 0;

#undef ic
#undef im

}

/*---------------------------------------------------------------------------*/

