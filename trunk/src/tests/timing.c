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

#include <time.h>       /** TODO: only one of these two is necessary **/
#include <sys/time.h>

#include <source_unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

/* define timer */
#ifdef HAVE_GETTIMEOFDAY
/* use gettimeofday() command. Not in ANSI C! */
static struct timeval tv;
#define _unur_get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )
#else
/* use clock() command. ANSI C but less accurate */
#define _unur_get_time() ( (1.e6 * clock()) / CLOCKS_PER_SEC )
#endif

/*---------------------------------------------------------------------------*/
static char test_name[] = "Timing";
/*---------------------------------------------------------------------------*/

static double _unur_test_timing_uniform( struct unur_par *par, int log_samplesize );

/*---------------------------------------------------------------------------*/

struct unur_gen*
unur_test_timing( struct unur_par *par, int log_samplesize )
     /*----------------------------------------------------------------------*/
     /*  init generator and estimate setup and generation time.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par            ... pointer to parameters for generator object      */
     /*   log_samplesize ... common log of maximal sample size               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  double *vec = NULL;
  int dim = 1;
  double fast;  
  double time_start, time_setup, time_marg, *time_gen;
  long samples, samplesize, log_samples;

  /* check parameter */
  _unur_check_NULL(test_name,par,NULL);
  if (log_samplesize < 2) log_samplesize = 2;

  /* need an array to store timings */
  time_gen = _unur_malloc((log_samplesize+1) * sizeof(double));

  /* fastest possible generation time */
  fast = _unur_test_timing_uniform( par,log_samplesize );

  /* initialize generator (and estimate setup time) */
  time_start = _unur_get_time();
  gen = _unur_init(par);
  time_setup = _unur_get_time();

  /* init successful ? */
  if (!gen) return NULL;

  /* we need an array for the vector */
  if (unur_is_vec(par)) {
    dim = unur_get_dimension(gen);
    vec = _unur_malloc( dim * sizeof(double) );
  }

  /* evaluate generation time */
  samplesize = 10;
  samples = 0;
  for( log_samples=1; log_samples<=log_samplesize; log_samples++ ) {

    switch (gen->method & UNUR_MASK_TYPE) {
    case UNUR_METH_DISCR:
      for( ; samples < samplesize; samples++ )
	_unur_sample_discr(gen);
      break;
    case UNUR_METH_CONT:
      for( ; samples < samplesize; samples++ )
	_unur_sample_cont(gen);
      break;
    case UNUR_METH_VEC:
      for( ; samples < samplesize; samples++ )
	_unur_sample_vec(gen,vec);
      break;
    default: /* unknown ! */
      _unur_error(test_name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }

    time_gen[log_samples] = _unur_get_time();
    samplesize *= 10;
  }

  /* compute generation times */

  /* marginal generation time */
  time_marg = (time_gen[log_samplesize] - time_gen[log_samplesize-1]) / (0.09 * samplesize);
  /* mean time per random number including setup */
  samplesize = 1;
  for( log_samples=1; log_samples<=log_samplesize; log_samples++ ) {
    samplesize *= 10;
    time_gen[log_samples] = (time_gen[log_samples] - time_start) / samplesize;
  }
  /* setup time */
  time_setup -= time_start;
  
  /* now print times */
  printf("\nTIMING:\t\t    usec \t relative to \t relative to\n");
  printf("\t\t\t\t uniform\t marginal\n\n");
  /* setup time */
  printf("   setup time:\t    %#g \t %#g \t %#g\n",time_setup,time_setup/fast,time_setup/time_marg);
  /* marginal generation time */
  printf("   generation time: %#g \t %#g \t %#g\n",time_marg,time_marg/fast,1.);
  /* generation times */
  printf("\n   average generation time for samplesize:\n");
  for( log_samples=1; log_samples<=log_samplesize; log_samples++ )
    printf("\t10^%ld:\t    %#g \t %#g \t %#g\n",log_samples,
	   time_gen[log_samples],time_gen[log_samples]/fast,time_gen[log_samples]/time_marg);

  /* free memory */
  free(time_gen);
  if (vec) free(vec);

  /* return generator object */
  return gen;

} /* end of unur_test_timing() */

/*---------------------------------------------------------------------------*/

static double
_unur_test_timing_uniform( struct unur_par *par, int log_samplesize )
     /*----------------------------------------------------------------------*/
     /*  estimate fastest possible generation time, i.e., for a routine      */
     /*  that just calls the uniform random number generator  and returns    */
     /*  this value.                                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   mean generation time                                               */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen_urng;
  static double fastest_time = 0.;
  int j;

  if (fastest_time <= 0.) {  
    /* unknown. have to estimate time first */
    
    /* sample size */
    int samplesize = 1;
    for( j=0; j<log_samplesize; j++ )
      samplesize *= 10;

    /* make generator object for uniform generator */
    gen_urng = _unur_unif_init( unur_unif_new(0,1) );
    _unur_check_NULL( test_name,gen_urng,-1. );
    unur_chg_urng(gen_urng,par->urng);

    /* evaluate marginal generation time */
    fastest_time = _unur_get_time();
    for( j=0; j<samplesize; j++ )
      _unur_sample_cont(gen_urng);

    fastest_time = (_unur_get_time() - fastest_time)/samplesize;

    /* free generator object for uniform random number generator */
    _unur_free(gen_urng);

  }

  return fastest_time;
} /* end of _unur_test_timing_uniform() */

/*---------------------------------------------------------------------------*/
