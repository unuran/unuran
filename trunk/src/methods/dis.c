/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dis.c                                                        *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    indexed search (guide table)                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given N discrete events with different probabilities P[k]            *
 *      produce a value k consistent with its probability.                   *
 *                                                                           *
 *   REQUIRED:  pointer to probability vector                                *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *      int factor  ... relative size of guide table (relative to N)         *
 *                      (default: 1)                                         *
 *      int variant ... indicates variant for setup procedure                *
 *                      (default: 0)                                         *
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Chen, H. C. and Asau, Y. (1974): On generating random variates      *
 *       from an empirical distribution, AIIE Trans. 6, pp. 163-166          *
 *                                                                           *
 *   SUMMARY:                                                                *
 *   [2] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This method a varient of the inversion method, i.e.                     *
 *                                                                           *
 *   (1) Generate a random number U ~ U(0,1).                                *
 *   (2) Find largest integer I such that F(I) = P(X<=I) <= U.               *
 *                                                                           *
 *   Step (2) is the crucial step. Using sequential search requires O(N)     *
 *   comparisons. Indexed search however uses a guide table to jump to some  *
 *   I' <= I near I. In a preprossing step (0,1) is partitioned into N       *
 *   equal intervals and P(X<=I) <= (k-1)/N is solved for all k=1,...,N,     *
 *   and the solutions are stored in a table (guide table). Setup this       *
 *   table can be done in O(N). [2] has shown that the expected number of    *
 *   of comparisons is at most 2.                                            *
 *                                                                           *
 *   In the current implementation we use a variant that allows a different  *
 *   size for the guide table. Bigger guide tables reduce the expected       *
 *   number of comparisions but needs more memory and need more setup time.  *
 *   On the other hand, the guide table can be made arbitrarily small to     *
 *   save memory and setup time. Indeed, for size = 1 we have sequential     *
 *   search again that requires no preprocessing.                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   VARIANTS:                                                               *
 *                                                                           *
 *   We have three variants for the setup procedure:                         *
 *                                                                           *
 *   1 ... compute (k-1)/N for each k=1,...,N                                *
 *   2 ... compute (k-1)/N by summing up 1/N                                 *
 *   0 ... use 1 for large N and 2 for small N (default)                     *
 *                                                                           *
 *   variant 2 is faster but is more sensitive to roundoff errors when       *
 *   N is large.                                                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Wed Sep  8 16:15:41 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "DIS"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dis_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dis_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dis_debug_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for guide table.                                               */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.dis
#define GEN     gen->data.dis
#define SAMPLE  gen->sample.discr

/*---------------------------------------------------------------------------*/
/* Special debugging flags (do not use the first 3 bits) */
#define DIS_DEBUG_TABLE   0x010UL

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dis_new( double *prob, int len )
     /*---------------------------------------------------------------------------*/
     /* get default parameters                                                    */
     /*                                                                           */
     /* parameters:                                                               */
     /*   prob ... pointer to probability vector                                  */
     /*   len  ... length of probability vector                                   */
     /*                                                                           */
     /* return:                                                                   */
     /*   default parameters (pointer to structure)                               */
     /*                                                                           */
     /* error:                                                                    */
     /*   return NULL                                                             */
     /*---------------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_DIS_PAR);

  /* copy input */
  PAR.prob         = prob;
  PAR.len          = len;

  /* set default values */
  PAR.guide_factor = 1.;             /* use same size for guide table        */

  par->method      = UNUR_METH_DIS;  /* method and default variant           */
  par->set         = 0UL;            /* inidicate default parameters         */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  _unur_set_debugflag_default(par);  /* set default debugging flags          */
  _unur_set_genid(par,GENTYPE);      /* set generator identifier             */

  /* routine for starting generator */
  par->init = unur_dis_init;

  return par;

} /* end of unur_dis_new() */

/*****************************************************************************/

struct unur_gen *
unur_dis_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  double probh, gstep;
  int i,j;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_DIS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dis_create(par);
  if (!gen) { free(par); return NULL; }

  /* computation of cumulated probabilities */
  for( i=0, probh=0.; i<PAR.len; i++ ) {
    GEN.cumprob[i] = ( probh += PAR.prob[i] );
#if CHECKARGS
    /* ... and check probability vector */
    if (PAR.prob[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"probability < 0 not possible!");
      unur_dis_free(gen); free(par); 
      return NULL;
    }
#endif
  }
  GEN.sum = GEN.cumprob[PAR.len-1];

  /* computation of guide-table */
  
  if ((gen->method & UNUR_MASK_VARIANT) == 1) {
    GEN.guide_table[0] = 0;
    for( j=1, i=0; j<GEN.guide_size ;j++ ) {
      while( GEN.cumprob[i]/GEN.sum < ((double)j)/GEN.guide_size ) 
	i++;
      if (i >= PAR.len) {
	_unur_warning(gen->genid,UNUR_ERR_INIT,"roundoff error while making guide table!");
	break;
      }
      GEN.guide_table[j]=i;
    }
  }

  else { /* variant 2 */
    gstep = GEN.sum / GEN.guide_size;
    probh = 0.;
    for( j=0, i=0; j<GEN.guide_size ;j++ ) {
      while (GEN.cumprob[i] < probh) 
	i++;
      if (i >= PAR.len) {
	_unur_warning(gen->genid,UNUR_ERR_INIT,"roundoff error while making guide table!");
	break;
      }
      GEN.guide_table[j] = i;
      probh += gstep;
    }
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide_table[j] = PAR.len - 1;

  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_dis_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;
} /* end of unur_dis_init() */

/*****************************************************************************/

int
unur_dis_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (sample from random variate)                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{ 
  int j;
  double u;

  /* check arguments */
  CHECK_NULL(gen,0);
  COOKIE_CHECK(gen,CK_DIS_GEN,0);

  /* sample from U(0,1) */
  u = _unur_call_urng(gen);

  /* look up in guide table ... */
  j = GEN.guide_table[(int)(u * GEN.guide_size)];
  /* ... and search */
  u *= GEN.sum;
  while (GEN.cumprob[j] < u) j++;

  return j;

} /* end of unur_dis_sample() */

/*****************************************************************************/

void
unur_dis_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if (!gen) /* nothing to do */
    return;

  /* magic cookies */
  COOKIE_CHECK(gen,CK_DIS_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.guide_table);
  free(GEN.cumprob);
  if (GEN.prob) free(GEN.prob);
  free(gen);

} /* end of unur_dis_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_dis_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  unsigned long variant;
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_DIS_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DIS_GEN);

  /* copy some parameters into generator object */
  GEN.len = PAR.len;                /* length of probability vector          */
  GEN.prob = NULL;                  /* copy probability vector on demand     */
  _unur_copy_urng_pointer(par,gen); /* pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);    /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);        /* copy generator identifier             */

  /* routines for sampling and destroying generator */
  SAMPLE = unur_dis_sample;
  gen->destroy = unur_dis_free;

  /* set all pointers to NULL */
  GEN.cumprob = NULL;
  GEN.guide_table = NULL;

  /* which variant? */
  variant = par->method & UNUR_MASK_VARIANT;
  if (variant > 2) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"invalid variant! use default");
    variant = 0;
  }
  if (!variant)   /* default variant */
    variant = (PAR.len > 1000) ? 1 : 2;
  /* store method in generator structure */
  gen->method = (par->method & (~UNUR_MASK_VARIANT)) | variant;

  /* allocation for cummulated probabilities */
  GEN.cumprob = _unur_malloc( PAR.len * sizeof(double) );

  /* size of guide table */
  GEN.guide_size = (int)(PAR.len * PAR.guide_factor);
  if (GEN.guide_size <= 0)
    /* do not use a guide table whenever params->guide_factor is 0 or less */
    GEN.guide_size = 1;

  /* allocate memory for the guide table */
  GEN.guide_table = _unur_malloc( GEN.guide_size * sizeof(int) );

  /* copy probability vector on demand */
  if (par->method & UNUR_MASK_COPYALL) {
    /* need an arry for the probabilty vector */
    GEN.prob = _unur_malloc( GEN.len * sizeof(double) );
    /* copy */
    for (i=0; i<GEN.len; i++)
      GEN.prob[i] = PAR.prob[i];
  }

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dis_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

/*---------------------------------------------------------------------------*/

static void
_unur_dis_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = indexed search (guide table)\n",gen->genid);
  fprintf(log,"%s: variant = %ld ",gen->genid,gen->method & UNUR_MASK_VARIANT);
  _unur_print_if_default(par,UNUR_SET_VARIANT);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_dis_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: length of probability vector = %d\n",gen->genid,PAR.len);
  fprintf(log,"%s: length of guide table = %d   (rel. = %g%%",
	  gen->genid,GEN.guide_size,100.*PAR.guide_factor);
  _unur_print_if_default(par,UNUR_SET_FACTOR);
  if (GEN.guide_size == 1) 
    fprintf(log,") \t (-->sequential search");
  fprintf(log,")\n%s:\n",gen->genid);

  if (gen->debug & DIS_DEBUG_TABLE) {
    _unur_dis_debug_table(gen);
    fprintf(log,"%s:\n",gen->genid);
  }

} /* end of _unur_dis_debug_init() */

/*---------------------------------------------------------------------------*/

static void
_unur_dis_debug_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write guide table into logfile                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{   
  FILE *log;
 
  int i,j,m;
  int n_asts;

  log = unur_get_stream();
  
  fprintf(log,"%s: guide table:\n", gen->genid); 
  fprintf(log,"%s:\n", gen->genid);
  n_asts = 0;
  for (i=0; i<GEN.guide_size; i++){
    fprintf(log,"%s: [%5d] -> %5d ", gen->genid, i, GEN.guide_table[i]);
    /* print row of asterisks */
    if (i == GEN.guide_size-1)
      j = GEN.guide_size - GEN.guide_table[i];
    else
      j = GEN.guide_table[i+1] - GEN.guide_table[i] + 1;
    for (m=0; m<j && m<10; m++ ) {
      fprintf(log," *");
      ++n_asts;
    }
    /* too many asterisks print */
    if (m<j){
      n_asts += j-m;
      fprintf(log," ... %d", j);
    }
    fprintf(log,"\n");
  }

  /* print expected number of comparisons */
  fprintf(log,"%s:\n", gen->genid);
  fprintf(log,"%s: expected number of comparisons = %g\n",gen->genid,
          ((double)n_asts)/GEN.guide_size);
  fprintf(log,"%s:\n", gen->genid);

} /*  end of _unur_dis_debug_table() */

/*---------------------------------------------------------------------------*/
#endif

/*****************************************************************************/
