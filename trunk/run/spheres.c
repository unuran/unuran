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

#define MAXDIM 100

#define VERBOSE  0     /* output switch for test and moments  */


#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/* default values */
int DIM = 2;
long SAMPLESIZE = 100;
long EXPERIMENTS = 1;
int METHOD = 1; 
int COORDINATE = 0; /* current coordinate for coordinate sampler */

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

int main(int argc, char *argv[])
{
  double *x; /* current point */
  double *direction; /* current direction */
  double lambda, lambda1, lambda2; /* direction parameters */
  double U; /* U(0,1) */
  double R; /* radius of inner sphere */
  long inside;
  double diff, var;
  
  int i;
  long loop, sample;
    
  char c;

  UNIFORM = unur_str2gen("uniform(0,1) & urng=mt19937(1173)");
  
  /* we need an generator for standard normal distributons */
  struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
  struct unur_par   *normalpar = unur_arou_new( normaldistr );
  unur_arou_set_usedars( normalpar, TRUE );
  NORMAL = unur_init( normalpar );
  _unur_distr_free( normaldistr );
   
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
      printf(" -m method       : 1=COORD, 2=Random Direction (%d)\n", METHOD );
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
  
  if (METHOD==METHOD_COORDINATE)        printf("METHOD=COORDINATE\n");
  if (METHOD==METHOD_RANDOM_DIRECTION)  printf("METHOD=RANDOM DIRECTION\n");
  
  R = exp(-log(2)/DIM); /* radius of interior sphere (volume = 1/2 of unit sphere volume)*/
  //printf("R=%e\n", R); 
  
  /* main loop */    
  var=0.; 
  for (loop=1; loop<=EXPERIMENTS; loop++) {

    /* setting initial point at the origin */
    for (i=0; i<DIM; i++) x[i] = 0.;
    inside=0; 
    
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
  
      if (_unur_vector_norm(DIM, x)<R) {
        inside++;
      }
    }
    diff = (inside -SAMPLESIZE/2.);
    var += diff*diff; 
  } /* next experiment */
  var /= EXPERIMENTS;
  
  /* output of results */   
  printf("dim=%d 'var'=%e\n", DIM, var); 
  printf("------------------------------------------------------------------\n");
      
  if (UNIFORM) unur_free(UNIFORM);
  
  free(x); free(direction);
      
  return 0;


}

/*---------------------------------------------------------------------------*/

