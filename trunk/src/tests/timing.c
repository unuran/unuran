/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      timing.c                                                     *
 *                                                                           *
 *   estimate setup and (marginal) generation time                           *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/
/* define timer */

#if defined(HAVE_GETTIMEOFDAY) && defined(HAVE_SYS_TIME_H)
/* use gettimeofday() command. Not in ANSI C! */
#include <sys/time.h>
static struct timeval tv;
#define _unur_get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )
#else
/* use clock() command. ANSI C but less accurate */
#include <time.h>
#define _unur_get_time() ( (1.e6 * clock()) / CLOCKS_PER_SEC )
#endif

/*---------------------------------------------------------------------------*/
static char test_name[] = "Timing";
/*---------------------------------------------------------------------------*/

#define TIME_BOUND   1000000
/*---------------------------------------------------------------------------*/
/*  approximate time for running timing test unur_test_timing_total_run()    */
/*---------------------------------------------------------------------------*/

double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int repeat );
/*---------------------------------------------------------------------------*/
/*  estimate average time (in micro seconds) for sampling 1 random variate   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* compare two doubles (needed for sorting) */
inline static int
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

/*---------------------------------------------------------------------------*/

struct unur_gen*
unur_test_timing( struct unur_par *par, 
		  int log_samplesize, 
		  double *time_setup,
		  double *time_sample,
		  int verbosity,
		  FILE *out )
     /*----------------------------------------------------------------------*/
     /*  init generator and estimate setup and generation time.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par            ... pointer to parameters for generator object      */
     /*   log_samplesize ... common log of maximal sample size               */
     /*   time_setup     ... time for setup                                  */
     /*   time_sample    ... marginal generation time (i.e. for one r.n.)    */
     /*   verbosity      ... verbosity level, 0 = no output, 1 = output      */
     /*   out            ... output stream                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object.                                       */
     /*   setup time and marginal generation time are stored in              */
     /*   setup_time and marginal_time, respectively.                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  int k;
  double x;
  double *vec = NULL;
  double time_uniform, time_exponential;
  double time_start, *time_gen;
  long samples, samplesize, log_samples;

  /* check parameter */
  _unur_check_NULL(test_name,par,NULL);
  if (log_samplesize < 2) log_samplesize = 2;

  /* need an array to store timings */
  time_gen = _unur_malloc((log_samplesize+1) * sizeof(double));

  /* marginal generation time for one unifrom random number */
  time_uniform = unur_test_timing_uniform( par,log_samplesize );
  /* marginal generation time for one exponential random variate */
  time_exponential = unur_test_timing_exponential( par,log_samplesize );

  /* we need an array for the vector */
  if (_unur_gen_is_vec(par))
    vec = _unur_malloc( par->distr->dim * sizeof(double) );

  /* initialize generator (and estimate setup time) */
  time_start = _unur_get_time();
  gen = _unur_init(par);
  *time_setup = _unur_get_time();

  /* init successful ? */
  if (!gen) {
    free (time_gen);
    return NULL;
  }

  /* evaluate generation time */
  samplesize = 10;
  samples = 0;
  for( log_samples=1; log_samples<=log_samplesize; log_samples++ ) {

    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( ; samples < samplesize; samples++ )
	k = unur_sample_discr(gen);
      break;
    case UNUR_METH_CONT:
      for( ; samples < samplesize; samples++ )
	x = unur_sample_cont(gen);
      break;
    case UNUR_METH_VEC:
      for( ; samples < samplesize; samples++ )
	unur_sample_vec(gen,vec);
      break;
    default: /* unknown ! */
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }

    time_gen[log_samples] = _unur_get_time();
    samplesize *= 10;
  }

  /* compute generation times */

  /* marginal generation time */
  *time_sample = (time_gen[log_samplesize] - time_gen[log_samplesize-1]) / (0.09 * samplesize);
  /* mean time per random number including setup */
  samplesize = 1;
  for( log_samples=1; log_samples<=log_samplesize; log_samples++ ) {
    samplesize *= 10;
    time_gen[log_samples] = (time_gen[log_samples] - time_start) / samplesize;
  }
  /* setup time */
  *time_setup -= time_start;
  
  /* now print times */
  if (verbosity) {
    fprintf(out,"\nTIMING:\t\t    usec \t relative to \t relative to\n");
    fprintf(out,"\t\t\t\t uniform\t exponential\n\n");
    /* setup time */
    fprintf(out,"   setup time:\t    %#g \t %#g \t %#g\n",
	    (*time_setup),
	    (*time_setup)/time_uniform,
	    (*time_setup)/time_exponential);
    /* marginal generation time */
    fprintf(out,"   generation time: %#g \t %#g \t %#g\n",
	    (*time_sample),
	    (*time_sample)/time_uniform,
	    (*time_sample)/time_exponential);
    /* generation times */
    fprintf(out,"\n   average generation time for samplesize:\n");
    for( log_samples=1; log_samples<=log_samplesize; log_samples++ )
      fprintf(out,"\t10^%ld:\t    %#g \t %#g \t %#g\n",log_samples,
	      time_gen[log_samples],
	      time_gen[log_samples]/time_uniform,
	      time_gen[log_samples]/time_exponential);
  }

  /* free memory */
  free(time_gen);
  if (vec) free(vec);

  /* return generator object */
  return gen;

} /* end of unur_test_timing() */

/*---------------------------------------------------------------------------*/

double 
unur_test_timing_total( const UNUR_PAR *par, int samplesize, double max_duration )
     /*----------------------------------------------------------------------*/
     /*  estimate average time (in micro seconds) for generating a sample    */
     /*  of size `samplesize' (including setup) are estimated.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameters for generator object        */
     /*   samplesize   ... sample size                                       */
     /*   max_duration ... upper bound for total time (in seconds) for       */
     /*                    running test                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   median of several runs for generating sample (in micro seconds)    */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1.                                                         */
     /*----------------------------------------------------------------------*/
{
#define REPEAT_0_PILOT 11         /* size of pilot study */

  double time_pilot;
  int repeat;

  /* check parameter */
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0) return -1.;
  if (max_duration < 1.e-3) max_duration = 1.e-3; 

  /* pilot study */
  time_pilot = unur_test_timing_total_run(par,1,REPEAT_0_PILOT);
  if (time_pilot < 0)
    return -1.;

  /* now run timing test */
  repeat = (int) (TIME_BOUND / time_pilot);
  if (repeat <= 11) {
    /* there is no need to run this test again */
    return time_pilot;
  }
  else {
    if (repeat > 1000) 
      /* there is no need for more than 1000 repetitions */
      repeat = 1000;
    return unur_test_timing_total_run(par,1,repeat);
  }

  /* this should not happen */
  return -1.;

#undef SIZE_0_PILOT
} /* end of unur_test_timing_total() */

/*---------------------------------------------------------------------------*/

double unur_test_timing_total_run( const struct unur_par *par, int samplesize, int repeat )
     /*----------------------------------------------------------------------*/
     /*  estimate average time (in micro seconds) for sampling               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par        ... pointer to parameters for generator object          */
     /*   samplesize ... sample size                                         */
     /*   repeat     ... number of samples (repetitions of sampling)         */
     /*                                                                      */
     /* return:                                                              */
     /*   total time in micro seconds                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par_tmp;    /* temporary working copy of parameter opject */
  struct unur_gen *gen_tmp;    /* temporary generator object                 */
  double *time;
  int n, rep;
  double time_total;           /* total time for sampling                    */
  double time_start;
  int k;
  double x;
  double *vec = NULL;

  /* check parameter */
  _unur_check_NULL(test_name,par,-1.);
  if (samplesize < 0 || repeat < 1) 
    return -1.;

  /* we need an array for storing timing results */
  time = _unur_malloc( repeat * sizeof(double) );

  /* we need an array for an random vector */
  if (_unur_gen_is_vec(par))
    vec = _unur_malloc( par->distr->dim * sizeof(double) );

  /* make samples */
  for (rep = 0; rep < repeat; rep++) {

    /* make a working copy of parameter object */
    par_tmp = _unur_malloc(sizeof(struct unur_par));
    memcpy (par_tmp, par, sizeof(struct unur_par));

    /* start timer */
    time_start = _unur_get_time();

    /* make generator object (init) */
    gen_tmp = _unur_init(par_tmp);
    if (!gen_tmp) return -1.;

    /* run generator */
    switch (gen_tmp->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( n=1; n<samplesize; n++ )
	k = unur_sample_discr(gen_tmp);
      break;
    case UNUR_METH_CONT:
      for( n=1; n<samplesize; n++ )
	x = unur_sample_cont(gen_tmp);
      break;
    case UNUR_METH_VEC:
      for( n=1; n<samplesize; n++ )
	unur_sample_vec(gen_tmp,vec);
      break;
    default: /* unknown ! */
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }

    /* stop timer */
    time[rep]= _unur_get_time() - time_start;

    /* destroy parameter object */
    unur_free(gen_tmp);
  }

  /* compute median */
  qsort( time, repeat, sizeof(double), compare_doubles);
  time_total = time[repeat/2];

  /* free memory */
  if (vec) free(vec);
  free(time);

  return time_total;
} /* end of unur_test_timing_total_run() */

/*---------------------------------------------------------------------------*/

double
unur_test_timing_uniform( struct unur_par *par, int log_samplesize )
     /*----------------------------------------------------------------------*/
     /*  estimate generation time for URNG using UNURAN wrapper.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   mean generation time                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define TIMING_REPETITIONS 11

  struct unur_gen *gen_urng;
  static double uniform_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;

  if (uniform_time <= 0.) {  
    /* unknown. have to estimate time first */
    
    /* sample size */
    int samplesize = 1;
    for( j=0; j<log_samplesize; j++ )
      samplesize *= 10;

    /* make generator object for uniform generator */
    gen_urng = unur_init( unur_unif_new(NULL) );
    _unur_check_NULL( test_name,gen_urng,-1. );
    unur_chg_urng(gen_urng,par->urng);

    /* evaluate marginal generation times */
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(gen_urng);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }

    /* compute median */
    qsort( time, TIMING_REPETITIONS, sizeof(double), compare_doubles);

    /* store marginal generation time for uniforms */
    uniform_time = time[TIMING_REPETITIONS/2];

    /* free generator object for uniform random number generator */
    _unur_free(gen_urng);

  }

  return uniform_time;

#undef TIMING_REPETITIONS

} /* end of unur_test_timing_uniform() */

/*---------------------------------------------------------------------------*/

double
unur_test_timing_exponential( struct unur_par *par, int log_samplesize )
     /*----------------------------------------------------------------------*/
     /*  estimate generation time for URNG using UNURAN wrapper.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   mean generation time                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
#define TIMING_REPETITIONS 11

  struct unur_distr *unit_distr;
  struct unur_par   *unit_par;
  struct unur_gen   *unit_gen;
  static double exponential_time = -1.;
  double time[TIMING_REPETITIONS];
  double x;
  int j,n;

  if (exponential_time <= 0.) {  
    /* unknown. have to estimate time first */
    
    /* sample size */
    int samplesize = 1;
    for( j=0; j<log_samplesize; j++ )
      samplesize *= 10;

    /* make generator object for uniform generator */
    unit_distr = unur_distr_exponential(NULL,0);
    unit_par = unur_cstd_new(unit_distr);
    unur_cstd_set_variant(unit_par,UNUR_STDGEN_INVERSION);
    unit_gen = unur_init(unit_par); 
    _unur_check_NULL( test_name,unit_gen,-1. );
    unur_chg_urng(unit_gen,par->urng);

    /* evaluate marginal generation times */
    for( n=0; n<TIMING_REPETITIONS; n++ ) {
      time[n] = _unur_get_time();
      for( j=0; j<samplesize; j++ )
	x = unur_sample_cont(unit_gen);
      time[n] = (_unur_get_time() - time[n])/samplesize;
    }

    /* compute median */
    qsort( time, TIMING_REPETITIONS, sizeof(double), compare_doubles);

    /* store marginal generation time for uniforms */
    exponential_time = time[TIMING_REPETITIONS/2];

    /* free generator object for uniform random number generator */
    unur_distr_free(unit_distr);
    unur_free(unit_gen);

  }

  return exponential_time;

#undef TIMING_REPETITIONS

} /* end of unur_test_timing_exponential() */

/*---------------------------------------------------------------------------*/
