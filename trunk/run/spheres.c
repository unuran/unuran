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
#include <methods/x_gen_source.h>

#define METHOD_COORDINATE 1
#define METHOD_RANDOM_DIRECTION 2
#define METHOD_REJECTION 3

#define MAXDIM 100

#define VERBOSE  0     /* output switch for test and moments  */
#define MATHEMATICA 0

#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/* default values */
int DIM = 2;
long SAMPLESIZE = 100;
long EXPERIMENTS = 1;
int METHOD = 1; 
int COORDINATE = 0; /* current coordinate for coordinate sampler */

FILE *fmath; /* mathematica output */


UNUR_GEN *UNIFORM=NULL;
#define NORMAL    UNIFORM->gen_aux        /* pointer to normal variate generator */

/*---------------------------------------------------------------------------*/

void
get_random_direction( double *direction)
     /*----------------------------------------------------------------------*/
     /* generte a random direction vector  (unit vector)                     */
     /*----------------------------------------------------------------------*/
{
  int d;

  for (d=0; d<DIM; d++) {
    direction[d] = unur_sample_cont(NORMAL);
  }  

  _unur_vector_normalize(DIM, direction);
}

/*---------------------------------------------------------------------------*/

void get_lambdas(double *x, double *direction, double *lambda1, double *lambda2) {
  double A, B, C, D;
  int i;
  
  /* solving quadratic equation for lambda1 and lambda2 */
  A=0; B=0; C=0;
  for (i=0; i<DIM; i++) {
    A += direction[i]*direction[i];
    B += 2.*direction[i]*x[i];
    C += x[i]*x[i];
  }
  C -= 1.;
  
  D = sqrt(B*B-4*A*C);
  
  *lambda1 = (-B-D)/(2*A);
  *lambda2 = (-B+D)/(2*A);

}


/*---------------------------------------------------------------------------*/


void math2_init() 
{
  fmath = fopen("spheres.txt","w");  
  fprintf(fmath,"x={");
}  
 
void math2_point(double *x) 
{ 
  int j;
    
  /* current point */
  fprintf(fmath,"{");
  for (j=0; j<DIM; j++) {
    fprintf(fmath,"%4.3f", x[j]);
    if (j<DIM-1) fprintf(fmath,",");
  }
  fprintf(fmath,"}");
      
  fprintf(fmath,",");
}

void math2_tail() 
{    
  int j;
  /* origin ... */
  fprintf(fmath,"{");
  for (j=0; j<DIM; j++) {
    fprintf(fmath,"%4.3f", 0.);
    if (j<DIM-1) fprintf(fmath,",");
  }  
  fprintf(fmath,"}");
    
  fprintf(fmath,"};\n");
  fclose(fmath);
}



/*---------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  double *x; /* current point */
  double *direction; /* current direction */
  double lambda, lambda1, lambda2; /* direction parameters */
  double U; /* U(0,1) */
  double R; /* radius of inner sphere */
  long inside;
  double diff, mse, r, sum_r, diff_r, mse_r;
  double expected_inside, expected_r;
  
  int i;
  long loop, sample;
    
  char c;

//  UNIFORM = unur_str2gen("uniform(0,1) & urng=mt19937(1173)");
  UNIFORM = unur_str2gen("uniform(0,1) & urng=mt19937(13731)");
  
  /* we need an generator for standard normal distributons */
  struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
  struct unur_par   *normalpar = unur_arou_new( normaldistr );
  unur_arou_set_usedars( normalpar, TRUE );
  NORMAL = unur_init( normalpar );
  _unur_distr_free( normaldistr );
  NORMAL->urng = UNIFORM->urng;
   
  /* read options */
  while ((c = getopt(argc, argv, "d:n:m:e:h")) != -1) {
    switch (c) {
    case 'm':     /* sampling method  */
      METHOD=atol(optarg);
      break;
    case 'd':     /* dimension */
      DIM=atoi(optarg);
      break;
    case 'n':     /* number of samples */
      SAMPLESIZE=(long) atof(optarg);
      break;
    case 'e':     /* number of experiments */
      EXPERIMENTS=(long) atof(optarg);
      break;
    
    case 'h':     /* help */
      printf("options\n" );
      printf(" -d dim          : dimension (%d) \n", DIM );
      printf(" -m method       : 1=COORD, 2=Random Direction 3=Rejection (%d)\n", METHOD );
      printf(" -n samplesize   : number of samples (%ld) \n", SAMPLESIZE );
      printf(" -e experiments  : number of repetitions (%ld) \n", EXPERIMENTS );
      
      return 0;
      break;
    default:
      break;
    }
  }
  
  // unur_set_default_debug(UNUR_DEBUG_OFF);
   
  x = _unur_vector_new(DIM); /* current point at origin */
  direction  = _unur_vector_new(DIM); /* sampling direction from current point  */
  
  //if (METHOD==METHOD_COORDINATE)        printf("METHOD=COORDINATE\n");
  //if (METHOD==METHOD_RANDOM_DIRECTION)  printf("METHOD=RANDOM DIRECTION\n");
  
  R = exp(-log(2)/DIM); /* radius of interior sphere (volume = 1/2 of unit sphere volume)*/
  //printf("R=%e\n", R); 
  
  expected_inside = SAMPLESIZE / 2.;
  expected_r = DIM/(DIM+1.);
  
#if MATHEMATICA 
  math2_init();
#endif
   
  /* main loop */    
  mse=0.; mse_r=0.;
  for (loop=1; loop<=EXPERIMENTS; loop++) {

    /* setting initial point at the origin */
    for (i=0; i<DIM; i++) x[i] = 0.;
    inside=0; 
    sum_r=0.;
    
    for (sample=1; sample<=SAMPLESIZE; sample++) {
    
      if (METHOD==METHOD_COORDINATE) {    
        for (i=0; i<DIM; i++) direction[i] = 0.;
	direction[COORDINATE]=1;
    
        /* cyclical coordinate directions */
        if (++COORDINATE >= DIM) COORDINATE=0; 
      }
      
      if (METHOD==METHOD_RANDOM_DIRECTION) {
        get_random_direction(direction);
      }

      get_lambdas(x, direction, &lambda1, &lambda2);
      
      U = unur_sample_cont(UNIFORM);
      lambda = lambda1 + U * (lambda2 - lambda1);
          
      /* set new point */
      for (i=0; i<DIM; i++) x[i]=x[i]+lambda*direction[i];

#if 1      
      if (METHOD==METHOD_REJECTION) {
        do {
	  for (i=0; i<DIM; i++) {
	    U = unur_sample_cont(UNIFORM);  
	    x[i] = 2*U-1.;
	  }
	} while (_unur_vector_norm(DIM, x)>=1);
      }
#endif


#if MATHEMATICA 
      math2_point(x);
#endif
      
      r = _unur_vector_norm(DIM, x);        
      sum_r += r;
      if (r<R) {
        inside++;
      }
    } /* next sample */
    
    
    diff = (inside - expected_inside);
    mse += diff*diff; 
  
    diff_r = (sum_r / SAMPLESIZE - expected_r) ;
    mse_r += diff_r*diff_r;
  
  } /* next experiment */
  mse /= EXPERIMENTS;
  mse_r /= EXPERIMENTS;
    
#if MATHEMATICA 
  math2_tail();
#endif
  
  /* output of results */   
  printf("dim=%2d rmse/e=%e rmse_r/e=%e\n", DIM, sqrt(mse)/expected_inside, sqrt(mse_r)/expected_r); 
      
  if (UNIFORM) unur_free(UNIFORM);
  
  free(x); free(direction);
      
  return 0;


}

/*---------------------------------------------------------------------------*/

