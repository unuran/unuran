/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dgt.c                                                        *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    guide table (indexed search)                                 *
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

#define DGT_VARFLAG_DIV     0x01u     /* compute guide table by division n/k */
#define DGT_VARFLAG_ADD     0x02u     /* compute guide table by adding       */

#define DGT_VAR_THRESHOLD   1000      /* above this value: use variant 1, else 2 */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DGT_DEBUG_PRINTVECTOR   0x00000100u
#define DGT_DEBUG_TABLE         0x00000200u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DGT_SET_GUIDEFACTOR    0x010u
#define DGT_SET_VARIANT        0x020u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DGT"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dgt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dgt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_debug_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for guide table.                                               */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr      /* data for distribution object      */

#define PAR       par->data.dgt         /* data for parameter object         */
#define GEN       gen->data.dgt         /* data for generator object         */
#define DISTR     gen->distr.data.discr  /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dgt_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pv == NULL) {
    /* There is no PV try to compute it.                         */
    if ( DISTR_IN.pmf
	 && ( ((DISTR_IN.domain[1] - DISTR_IN.domain[0]) < UNUR_MAX_AUTO_PV)
	      || ( (distr->set & UNUR_DISTR_SET_PMFSUM) && DISTR_IN.domain[0] > INT_MIN ) ) ) {
      /* However this requires a PMF and either a bounded domain   */
      /* or the sum over the PMF.                                  */
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV. Try to compute it.");
    }
    else {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); return NULL;
    }
  }

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_DGT_PAR);

  /* copy input */
  par->distr       = distr;          /* pointer to distribution object       */

  /* set default values */
  PAR.guide_factor = 1.;             /* use same size for guide table        */

  par->method      = UNUR_METH_DGT;  /* method                               */
  par->variant     = 0u;             /* default variant                      */
  par->set         = 0u;             /* inidicate default parameters         */    
  par->urng        = unur_get_default_urng(); /* use default urng            */
  par->urng_aux    = NULL;                    /* no auxilliary URNG required */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dgt_init;

  return par;

} /* end of unur_dgt_new() */

/*---------------------------------------------------------------------------*/

int
unur_dgt_set_variant( struct unur_par *par, unsigned variant )
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
  _unur_check_par_object( par,DGT );

  /* check new parameter for generator */
  if (variant != DGT_VARFLAG_ADD && variant != DGT_VARFLAG_DIV) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
    return 0;
  }

  /* changelog */
  par->set |= DGT_SET_VARIANT;

  par->variant = variant;

  return 1;
} /* end of unur_dgt_set_variant() */

/*---------------------------------------------------------------------------*/

int
unur_dgt_set_guidefactor( struct unur_par *par, double factor )
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
  _unur_check_par_object( par,DGT );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"relative table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= DGT_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_dgt_set_guidefactor() */

/*****************************************************************************/

struct unur_gen *
_unur_dgt_init( struct unur_par *par )
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
  double *pv;                   /* pointer to probability vector */
  int n_pv;                     /* length of probability vector */
  double pvh;                   /* aux variable for computing cummalate sums */
  double gstep;                 /* step size when computing guide table */
  int i,j;
  
  /* check arguments */
  CHECK_NULL(par,NULL);
  
  /* check input */
  if ( par->method != UNUR_METH_DGT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DGT_PAR,NULL);
  
  /* create a new empty generator object */
  gen = _unur_dgt_create(par);
  if (!gen) { free(par); return NULL; }
  
  /* probability vector */
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;

  /* computation of cumulated probabilities */
  for( i=0, pvh=0.; i<n_pv; i++ ) {
    GEN.cumpv[i] = ( pvh += pv[i] );
    /* ... and check probability vector */
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      _unur_dgt_free(gen); free(par); 
      return NULL;
    }
  }
  GEN.sum = GEN.cumpv[n_pv-1];

  /* computation of guide-table */
  
  if (gen->variant == DGT_VARFLAG_DIV) {
    GEN.guide_table[0] = 0;
    for( j=1, i=0; j<GEN.guide_size ;j++ ) {
      while( GEN.cumpv[i]/GEN.sum < ((double)j)/GEN.guide_size ) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN.guide_table[j]=i;
    }
  }

  else { /* gen->variant == DGT_VARFLAG_ADD */
    gstep = GEN.sum / GEN.guide_size;
    pvh = 0.;
    for( j=0, i=0; j<GEN.guide_size ;j++ ) {
      while (GEN.cumpv[i] < pvh) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN.guide_table[j] = i;
      pvh += gstep;
    }
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide_table[j] = n_pv - 1;

  /* write info into log file */
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dgt_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;
} /* end of _unur_dgt_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_dgt_create( struct unur_par *par )
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
  int n_pv;                   /* length of probability vector */

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DGT_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DGT_GEN);

  /* copy distribution object into generator object */
  _unur_distr_discr_copy( &(gen->distr), par->distr );

  /* we need a PV */
  if (DISTR.pv == NULL) {
    /* try to compute PV */
    if (unur_distr_discr_make_pv(&(gen->distr)) <= 0) {
      /* not successful */
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); 
      if (DISTR.pv) free(DISTR.pv);
      free(gen);
      return NULL;
    }
  }

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dgt_sample;
  gen->destroy = _unur_dgt_free;
  gen->clone = _unur_dgt_clone;

  /* set all pointers to NULL */
  GEN.cumpv = NULL;
  GEN.guide_table = NULL;

  /* copy some parameters into generator object */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */

  /* length of probability vector */
  n_pv = DISTR.n_pv;

  /* store method in generator structure */
  gen->method = par->method;
  gen->set = par->set;              /* indicates parameter settings          */

  /* which variant? */
  if (par->variant == 0)   /* default variant */
    par->variant = (n_pv > DGT_VAR_THRESHOLD) ? DGT_VARFLAG_DIV : DGT_VARFLAG_ADD;
  /* store variant in generator structure */
  gen->variant = par->variant;

  /* allocation for cummulated probabilities */
  GEN.cumpv = _unur_malloc( n_pv * sizeof(double) );

  /* size of guide table */
  GEN.guide_size = (int)( n_pv * PAR.guide_factor);
  if (GEN.guide_size <= 0)
    /* do not use a guide table whenever params->guide_factor is 0 or less */
    GEN.guide_size = 1;

  /* allocate memory for the guide table */
  GEN.guide_table = _unur_malloc( GEN.guide_size * sizeof(int) );

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dgt_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dgt_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE clone->data.dgt

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DGT_GEN,NULL);

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  _unur_distr_discr_copy( &(clone->distr), &(gen->distr) );

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = unur_gen_clone( gen->gen_aux );

  /* copy data for distribution */
  CLONE.cumpv = _unur_malloc( DISTR.n_pv * sizeof(double) );
  memcpy( CLONE.cumpv, GEN.cumpv, DISTR.n_pv * sizeof(double) );
  CLONE.guide_table = _unur_malloc( GEN.guide_size * sizeof(int) );
  memcpy( CLONE.guide_table, GEN.guide_table, GEN.guide_size * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_dgt_clone() */

/*****************************************************************************/

int
_unur_dgt_sample( struct unur_gen *gen )
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
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_DGT_GEN,0);

  /* sample from U(0,1) */
  u = _unur_call_urng(gen->urng);

  /* look up in guide table ... */
  j = GEN.guide_table[(int)(u * GEN.guide_size)];
  /* ... and search */
  u *= GEN.sum;
  while (GEN.cumpv[j] < u) j++;

  return (j + DISTR.domain[0]);

} /* end of _unur_dgt_sample() */

/*****************************************************************************/

void
_unur_dgt_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DGT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DGT_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free two auxiliary tables */
  if (GEN.guide_table) free(GEN.guide_table);
  if (GEN.cumpv)       free(GEN.cumpv);

  /* free memory */
  _unur_free_genid(gen);
  _unur_distr_discr_clear(gen);

  COOKIE_CLEAR(gen);
  free(gen);

} /* end of _unur_dgt_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dgt_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_DGT_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DGT_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = indexed search (guide table)\n",gen->genid);

  fprintf(log,"%s: variant = %d ",gen->genid,gen->variant);
  _unur_print_if_default(par,DGT_SET_VARIANT);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_discr_debug( &(gen->distr),gen->genid,(gen->debug & DGT_DEBUG_PRINTVECTOR));

  fprintf(log,"%s: sampling routine = _unur_dgt_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: length of probability vector = %d\n",gen->genid,DISTR.n_pv);
  fprintf(log,"%s: length of guide table = %d   (rel. = %g%%",
	  gen->genid,GEN.guide_size,100.*PAR.guide_factor);
  _unur_print_if_default(par,DGT_SET_GUIDEFACTOR);
  if (GEN.guide_size == 1) 
    fprintf(log,") \t (-->sequential search");
  fprintf(log,")\n%s:\n",gen->genid);

  if (gen->debug & DGT_DEBUG_TABLE)
    _unur_dgt_debug_table(gen);

} /* end of _unur_dgt_debug_init() */

/*---------------------------------------------------------------------------*/

static void
_unur_dgt_debug_table( struct unur_gen *gen )
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DGT_GEN,/*void*/);

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

} /*  end of _unur_dgt_debug_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
