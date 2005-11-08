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

#include <distr/distr_source.h>
#include <src/utils/matrix_source.h>
#include <experimental/hitrou.h>
#include <experimental/gibbs.h>
#include <experimental/ball.h>
#include <experimental/walk.h>
#include <methods/x_gen_source.h>
#include <tests/testdistributions/testdistributions.h>
#include "meanvarcor.c"
/*#include <src/utils/fft.c>*/

#define DISTRIBUTION_NORMAL 0 
#define DISTRIBUTION_STUDENT 1
#define DISTRIBUTION_CAUCHY_BALL 2

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
#define METHOD_BALL_ROU 7
#define METHOD_BALL_ROU_ADAPTIVE 8
#define METHOD_BALL_PDF 9
#define METHOD_BALL_PDF_ADAPTIVE 10
#define METHOD_WALK 11
#define METHOD_HITROU_UNIDIRECTIONAL 12

#define MAXDIM 100

#define VERBOSE  0     /* output switch for test and moments  */

#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

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
  double uv[100];
  FILE *fx;
  FILE *fuv;
 
  uv[0]=0; /* init ... not necessary */ 

  fx = fopen("x.txt","w");
  fuv = fopen("uv.txt","w");
  
    fprintf(fx,"x={");
    fprintf(fuv,"uv={");
    for (i=1; i<=SAMPLESIZE; i++) {
      unur_sample_vec(gen, x);  
#if 0            
      _unur_hitrou_get_point_current( gen, uv);
#endif      
      /* current point in x-y coordinates */
      fprintf(fx,"{");
      for (j=0; j<DIM; j++) {
        fprintf(fx,"%4.3f", x[j]);
	if (j<DIM-1) fprintf(fx,",");
      }
      fprintf(fx,"}");
      
#if 0
      /* current point in u-v coordinates */
      fprintf(fuv,"{");
      for (j=0; j<=DIM; j++) {
        fprintf(fuv,"%4.3f", uv[j]);
	if (j<DIM) fprintf(fuv,",");
      }
      fprintf(fuv,"}");
#endif      
      
      
      if (i<SAMPLESIZE ) fprintf(fx,",");
      if (i<SAMPLESIZE ) fprintf(fuv,",");
    }
    fprintf(fx,"};\n");
    fprintf(fuv,"};\n");

    fclose(fx);
    fclose(fuv);
}



int main(int argc, char *argv[])
{
#define ic(a,b)  ((a)*DIM+(b))  /* index for covariance matrix */
#define im(a,b)  ((a)*5+(b))    /* index for moments */

  UNUR_DISTR *distr=NULL;
  UNUR_PAR *par=NULL;
  UNUR_GEN *gen=NULL;

  UNUR_PAR *par_clone=NULL;  
  
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

/* discrete fft tests ... 

  int n=100;
  int ifac[n];
  float wsave[2*n+15];
  float r[n];
  _fdrffti( n,  wsave, ifac);
  
  for (i=2;i<=ifac[1]+1; i++) printf("ifac[% d]=%d\n", i, ifac[i]); 
  for (i=0; i<100; i++) r[i]=i % 17;
  for (i=0; i<100; i++) printf("r[% d]=%f\n", i, r[i]);
  
  _fdrfftf(n, r, wsave, ifac);
  for (i=0; i<100; i++) printf("r[% d]=%f\n", i, r[i]);

*/  
    
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
    
    case 'h':     /*[GWa92] help */
      printf("options\n" );
      printf(" -d dim          : dimension (%d) \n", DIM );
      printf(" -t type         : 0=normal, 1=student 2=cauchy_ball (%d) \n", DISTRIBUTION );
      printf(" -f nu           : degrees of freedom for student (%d) \n", NU );
      printf(" -m method       : 0=H&R+RD+STRIP, 1=VMT, 2=GIBBS 3=GIBBS+RD \n" );
      printf("                 : 4=H&R+COORD+BOX, 5=H&R+RD+BOX 6=H&R+RD+ADAPTIVE STRIP\n" );
      printf("                 : 7=BALL+RoU  8=BALL+RoU+ADAPTIVE RADIUS \n");
      printf("                 : 9=BALL+PDF 10=BALL+PDF+ADAPTIVE RADIUS \n");
      printf("                 : 11=WALK  12=H&R+RD+STRIP+UNI (%d)\n", METHOD );
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
//  unur_set_default_debug(UNUR_DEBUG_ALL);
     
  for(d=0;d<DIM;d++){
    for (m=1; m<=4; m++) {
       mv[im(d,m)]=init_meanvar();
       if (EXPERIMENTS>1) quant[im(d,m)] = init_quantile(EXPERIMENTS);
    }
  }    
  
  x = _unur_vector_new(DIM+1); /* we need the extra '+1' for the pdf ball sampler */
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

  if (DISTRIBUTION==DISTRIBUTION_CAUCHY_BALL) {
    printf("DISTRIBUTION='CAUCHY_BALL'");
    distr = unur_distr_multicauchy_RoU_ball(DIM);
    for (d=0; d<DIM; d++) {
      moments_expected[im(d,1)] = INFINITY;
      moments_expected[im(d,2)] = INFINITY;
      moments_expected[im(d,3)] = INFINITY;
      moments_expected[im(d,4)] = INFINITY;    
    }
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
    unur_hitrou_set_skip(par,SKIP);
  }  

  if (METHOD==METHOD_HITROU_UNIDIRECTIONAL) {
    printf("METHOD=HITROU (UNIDIRECTIONAL)\n");
    par = unur_hitrou_new(distr);
    unur_hitrou_set_variant_random_direction(par);
    unur_hitrou_set_adaptive_strip(par, 0);    
    unur_hitrou_use_bounding_rectangle(par,0);
    unur_hitrou_set_unidirectional(par, 1);    
    unur_hitrou_set_skip(par,SKIP);
  }  
  
  
  if (METHOD==METHOD_BALL_ROU) {
    printf("METHOD=BALL (ROU)\n");
    printf("RADIUS=%f\n", BALL_RADIUS);
    par = unur_ball_new(distr);
    unur_ball_set_variant_rou(par);
    unur_ball_set_ball_radius(par, BALL_RADIUS);
    unur_ball_set_skip(par,SKIP);
  }  
  
  if (METHOD==METHOD_BALL_ROU_ADAPTIVE) {
    printf("METHOD=BALL (ROU + ADAPTIVE RADIUS)\n");
    printf("INITIAL_RADIUS=%f\n", BALL_RADIUS);
    par = unur_ball_new(distr);
    unur_ball_set_variant_rou(par);
    unur_ball_set_ball_radius(par, BALL_RADIUS);
    unur_ball_set_adaptive_ball(par, 1);    
    unur_ball_set_skip(par,SKIP);
  }  

  if (METHOD==METHOD_BALL_PDF) {
    printf("METHOD=BALL (PDF)\n");
    printf("RADIUS=%f\n", BALL_RADIUS);
    par = unur_ball_new(distr);
    unur_ball_set_variant_pdf(par);
    unur_ball_set_ball_radius(par, BALL_RADIUS);
    unur_ball_set_skip(par,SKIP);
  }  
  
  if (METHOD==METHOD_BALL_PDF_ADAPTIVE) {
    printf("METHOD=BALL (PDF + ADAPTIVE RADIUS)\n");
    printf("INITIAL_RADIUS=%f\n", BALL_RADIUS);
    par = unur_ball_new(distr);
    unur_ball_set_variant_pdf(par);
    unur_ball_set_ball_radius(par, BALL_RADIUS);
    unur_ball_set_adaptive_ball(par, 1);    
    unur_ball_set_skip(par,SKIP);
  }  

  if (METHOD==METHOD_WALK) {
    printf("METHOD=WALK \n");
    printf("INITIAL_RADIUS=%f\n", BALL_RADIUS);
    par = unur_walk_new(distr);
    unur_walk_set_ball_radius(par, BALL_RADIUS);
    unur_walk_set_skip(par,SKIP);
  }  
  
  printf("SKIP=%ld\n", SKIP);

  par_clone = _unur_par_clone(par);
 
#if 0
  gen = unur_init(par);
  math2(gen);
  exit(0);
#else  
  gen = unur_test_timing(par, 3, &time_setup, &time_sample, TRUE, stdout);
#endif
    
  double *uv;
  double *eigenvalues;
  double *eigenvectors;
  
  eigenvalues = _unur_vector_new(DIM);
  eigenvectors = _unur_vector_new(DIM*DIM);
  _unur_matrix_eigensystem(DIM, covar, eigenvalues, eigenvectors);
  uv = _unur_vector_new(DIM +1);
  
  /* setting x to be (normalized) eigenvector corresponding to largest eigenvalue */
  for (d=0; d<DIM; d++) {
    x[d]=eigenvectors[DIM*(DIM-1)-1+d];
  }  
  double fx;  /* pdf value at the point x */
  fx=PDF(x);  
  x[DIM] = fx * 0.9; /* used by the pdf ball sampler */
  /* transforming x[] into uv[] coordinates (at RoU boundary */
  /* valid for the case, that the RoU r-parameter is 1 */
  double V;
  V = pow(PDF(x), 1./(DIM + 1.));
  uv[DIM] = V;
  for (d=0; d<DIM; d++) {
        uv[d] = x[d] * V;
  } 
  
  /* scaling uv[]-coordinates */
  for (d=0; d<=DIM; d++) {
    uv[d] *= 0.9;
  }
    
  par = _unur_par_clone(par_clone);
  gen = unur_init(par);

  /* main loop */    
  for (loop=1; loop<=EXPERIMENTS; loop++) {
      
    if (METHOD==METHOD_HITROU 
    || METHOD==METHOD_HITROU_BOX 
    || METHOD==METHOD_HITROU_BOX_COORDINATE
    || METHOD==METHOD_HITROU_UNIDIRECTIONAL
    || METHOD==METHOD_HITROU_STRIP_ADAPTIVE) {
      _unur_hitrou_set_point_current( gen, uv );
    }
    
    if (METHOD==METHOD_BALL_ROU
    || METHOD==METHOD_BALL_ROU_ADAPTIVE) { 
      _unur_ball_set_point_current( gen, uv );
    }
    
    if (METHOD==METHOD_BALL_PDF
    || METHOD==METHOD_BALL_PDF_ADAPTIVE) { 
      _unur_ball_set_point_current( gen, x );
    }
    
    if (METHOD==METHOD_GIBBS
    || METHOD==METHOD_GIBBS_RANDOM) {
      _unur_gibbs_set_point_current( gen, x);    
    }
    
    if (METHOD==METHOD_WALK) {    
      _unur_walk_set_point_current( gen, x);    
    }
    
    unur_test_moments(gen, moments, 4, SAMPLESIZE, VERBOSE, stdout);
    
    for(d=0;d<DIM;d++){
      s=(d<DIM/2.) ? 1: sqrt(SIGMA);
      
      for (m=1; m<=4; m++) {
        update_meanvar(mv[im(d,m)],moments[im(d,m)]/pow(s,m));
        if (EXPERIMENTS>1) update_quantile(quant[im(d,m)],moments[im(d,m)]/pow(s,m));
      }
    }    
    
  
  } /* next experiment */

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
    
  unur_test_count_pdf(gen, SAMPLESIZE, 2, stdout);
    
  for(d=0;d<DIM;d++){
    for (m=1; m<=4; m++) {
      if (mv[im(d,m)]) free_meanvar(mv[im(d,m)]);
      if (EXPERIMENTS>1) 
        if (quant[im(d,m)]) free_quantile(quant[im(d,m)]);
    }
  }
    
  if (par_clone) unur_par_free(par_clone);
  
  if (distr) unur_distr_free(distr);
  if (gen)   unur_free(gen);
  
  if (x) free(x); 
  if (mean) free(mean); 
  if (covar) free(covar);
  if (moments) free(moments); 
  if (moments_expected) free(moments_expected); 
    
  if (uv) free(uv); 
  if (eigenvalues) free(eigenvalues); 
  if (eigenvectors) free(eigenvectors);
  
  return 0;

#undef ic
#undef im

}

/*---------------------------------------------------------------------------*/

