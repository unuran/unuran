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
 *****************************************************************************
 *****************************************************************************
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DIS_VAR_DIV         0x01u     /* compute guide table by division n/k */
#define DIS_VAR_ADD         0x02u     /* compute guide table by adding       */

#define DIS_VAR_THRESHOLD   1000      /* above this value: use variant 1, else 2 */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DIS_DEBUG_PRINTVECTOR   0x00000100u
#define DIS_DEBUG_TABLE         0x00000200u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DIS_SET_GUIDEFACTOR    0x010u
#define DIS_SET_VARIANT        0x020u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DIS"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dis_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
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

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       par->data.dis         /* data for parameter object         */
#define GEN       gen->data.dis         /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dis_new( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.prob == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.v."); return NULL;
  }

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_DIS_PAR);

  /* copy input */
  par->distr       = distr;          /* pointer to distribution object       */

  /* set default values */
  PAR.guide_factor = 1.;             /* use same size for guide table        */

  par->method      = UNUR_METH_DIS;  /* method                               */
  par->variant     = 0u;             /* default variant                      */
  par->set         = 0u;             /* inidicate default parameters         */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  par->genid    = _unur_set_genid(GENTYPE);/* set generator id               */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = unur_dis_init;

  return par;

} /* end of unur_dis_new() */

/*---------------------------------------------------------------------------*/

int
unur_dis_set_variant( struct unur_par *par, unsigned variant )
     /*----------------------------------------------------------------------*/
     /* set variant of method                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   variant ... indicator for variant                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( DIS );

  /* check new parameter for generator */
  if (variant != DIS_VAR_ADD && variant != DIS_VAR_DIV) {
    _unur_warning(par->genid,UNUR_ERR_PAR_VARIANT,"");
    return 0;
  }

  /* changelog */
  par->set |= DIS_SET_VARIANT;

  par->variant = variant;

  return 1;
} /* end of unur_dis_set_variant() */

/*---------------------------------------------------------------------------*/

int
unur_dis_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( DIS );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(par->genid,UNUR_ERR_PAR_SET,"relative table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= DIS_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_dis_set_guidefactor() */

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
  struct unur_gen *gen;         /* pointer to generator object */
  double *prob;                 /* pointer to probability vector */
  int n_prob;                   /* length of probability vector */
  double probh;                 /* aux variable for computing cummalate sums */
  double gstep;                 /* step size when computing guide table */
  int i,j;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_DIS ) {
    _unur_error(par->genid,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DIS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dis_create(par);
  if (!gen) { free(par); return NULL; }

  /* probability vector */
  prob = DISTR.prob;
  n_prob = DISTR.n_prob;

  /* computation of cumulated probabilities */
  for( i=0, probh=0.; i<n_prob; i++ ) {
    GEN.cumprob[i] = ( probh += prob[i] );
    /* ... and check probability vector */
    if (prob[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      unur_dis_free(gen); free(par); 
      return NULL;
    }
  }
  GEN.sum = GEN.cumprob[n_prob-1];

  /* computation of guide-table */
  
  if (gen->variant == DIS_VAR_DIV) {
    GEN.guide_table[0] = 0;
    for( j=1, i=0; j<GEN.guide_size ;j++ ) {
      while( GEN.cumprob[i]/GEN.sum < ((double)j)/GEN.guide_size ) 
	i++;
      if (i >= n_prob) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN.guide_table[j]=i;
    }
  }

  else { /* gen->variant == DIS_VAR_ADD */
    gstep = GEN.sum / GEN.guide_size;
    probh = 0.;
    for( j=0, i=0; j<GEN.guide_size ;j++ ) {
      while (GEN.cumprob[i] < probh) 
	i++;
      if (i >= n_prob) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN.guide_table[j] = i;
      probh += gstep;
    }
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide_table[j] = n_prob - 1;

  /* write info into log file */
#ifdef UNUR_ENABLE_LOGGING
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
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DIS_GEN,0);

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

  /* check input */
  if ( gen->method != UNUR_METH_DIS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DIS_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.guide_table);
  free(GEN.cumprob);
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
  struct unur_gen *gen;       /* pointer to generator object */
  int n_prob;                 /* length of probability vector */

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DIS_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DIS_GEN);

  /* copy generator identifier */
  gen->genid = par->genid;

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* routines for sampling and destroying generator */
  SAMPLE = unur_dis_sample;
  gen->destroy = unur_dis_free;

  /* set all pointers to NULL */
  GEN.cumprob = NULL;
  GEN.guide_table = NULL;

  /* copy some parameters into generator object */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  /* length of probability vector */
  n_prob = DISTR.n_prob;

  /* store method in generator structure */
  gen->method = par->method;

  /* which variant? */
  if (par->variant == 0)   /* default variant */
    par->variant = (n_prob > DIS_VAR_THRESHOLD) ? DIS_VAR_DIV : DIS_VAR_ADD;
  /* store variant in generator structure */
  gen->variant = par->variant;

  /* allocation for cummulated probabilities */
  GEN.cumprob = _unur_malloc( n_prob * sizeof(double) );

  /* size of guide table */
  GEN.guide_size = (int)( n_prob * PAR.guide_factor);
  if (GEN.guide_size <= 0)
    /* do not use a guide table whenever params->guide_factor is 0 or less */
    GEN.guide_size = 1;

  /* allocate memory for the guide table */
  GEN.guide_table = _unur_malloc( GEN.guide_size * sizeof(int) );

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dis_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
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

  /* check arguments */
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_DIS_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DIS_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = indexed search (guide table)\n",gen->genid);

  fprintf(log,"%s: variant = %d ",gen->genid,gen->variant);
  _unur_print_if_default(par,DIS_SET_VARIANT);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_discr_debug( &(gen->distr),gen->genid,(gen->debug & DIS_DEBUG_PRINTVECTOR));

  fprintf(log,"%s: sampling routine = unur_dis_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: length of probability vector = %d\n",gen->genid,DISTR.n_prob);
  fprintf(log,"%s: length of guide table = %d   (rel. = %g%%",
	  gen->genid,GEN.guide_size,100.*PAR.guide_factor);
  _unur_print_if_default(par,DIS_SET_GUIDEFACTOR);
  if (GEN.guide_size == 1) 
    fprintf(log,") \t (-->sequential search");
  fprintf(log,")\n%s:\n",gen->genid);

  if (gen->debug & DIS_DEBUG_TABLE)
    _unur_dis_debug_table(gen);

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

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DIS_GEN,/*void*/);

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

  fprintf(log,"%s:\n",gen->genid);

} /*  end of _unur_dis_debug_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
