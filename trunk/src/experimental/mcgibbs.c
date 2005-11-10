/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mcgibbs.c                                                    *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    Gibbs sampler using full conditional distributions.          *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF                                                            *
 *      Produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function and its derivatives                  *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      mode of the density                                                  *
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
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/condi.h>
#include <distr/cvec.h>
/* #include <distributions/unur_distributions.h> */
#include <uniform/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "tdrgw.h"
/* #include "auto.h" */

#include "mcgibbs.h"
#include "mcgibbs_struct.h"


/** TODO: **/
#define CK_MCGIBBS_PAR  CK_GIBBS_PAR 
#define CK_MCGIBBS_GEN  CK_GIBBS_GEN 

#define UNUR_METH_MCGIBBS  UNUR_METH_GIBBS 

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define MCGIBBS_VARMASK_VARIANT          0x000fu /* indicates variant        */
#define MCGIBBS_VARIANT_COORDINATE       0x0001u /* coordinate sampler       */
#define MCGIBBS_VARIANT_RANDOM_DIRECTION 0x0002u /* random direction sampler */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define MCGIBBS_DEBUG_CONDI   0x01000000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define MCGIBBS_SET_THINNING   0x001u  /* set thinning factor                */
#define MCGIBBS_SET_X0         0x002u  /* set starting point                 */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MCGIBBS"      /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcgibbs_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcgibbs_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void _unur_mcgibbs_coord_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_mcgibbs_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcgibbs_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mcgibbs_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_mcgibbs_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mcgibbs_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

/* generators for conditional distributions */
#define GEN_CONDI     gen->gen_aux_list     

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mcgibbs_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF");
    return NULL;
  }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dlogPDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mcgibbs_par) );
  COOKIE_SET(par,CK_MCGIBBS_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_MCGIBBS ;     /* method                          */
  par->variant  = MCGIBBS_VARIANT_COORDINATE;  /* default variant            */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  PAR->thinning = 1;                /* thinning factor                       */
  PAR->x0       = NULL;             /* starting point of chain, default is 0 */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_mcgibbs_init;

  return par;

} /* end of unur_mcgibbs_new() */

/*****************************************************************************/

int
unur_gibbs_set_variant_coordinate( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Coordinate Sampler :                                                 */
     /* Sampling along the coordinate directions (cyclic).                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCGIBBS );

  /* we use a bit in variant */
  par->variant = (par->variant & ~MCGIBBS_VARMASK_VARIANT) | MCGIBBS_VARIANT_COORDINATE;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_mcgibbs_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int
unur_gibbs_set_variant_random_direction( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Random Direction Sampler :                                           */
     /* Sampling along the random directions.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCGIBBS );
            
  /* we use a bit in variant */
  par->variant = (par->variant & ~MCGIBBS_VARMASK_VARIANT) | MCGIBBS_VARIANT_RANDOM_DIRECTION;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_mcgibbs_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int 
unur_mcgibbs_set_startingpoint( struct unur_par *par, const double *x0)
     /*----------------------------------------------------------------------*/
     /* set thinning factor for chain                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   x0       ... starting point of chain                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCGIBBS );

  /* store data */
  PAR->x0 = x0;

  /* changelog */
  par->set |= MCGIBBS_SET_X0;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_mcgibbs_set_startingpoint() */

/*---------------------------------------------------------------------------*/

int
unur_mcgibbs_set_thinning( struct unur_par *par, int thinning )
     /*----------------------------------------------------------------------*/
     /* set thinning factor for chain                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   thinning ... thinning factor                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCGIBBS );

  /* check new parameter for generator */
  if (thinning < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->thinning = thinning;

  /* changelog */
  par->set |= MCGIBBS_SET_THINNING;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_mcgibbs_set_thinning() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_mcgibbs_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  struct unur_par *par_condi;
  struct unur_gen *gen_condi;
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_MCGIBBS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MCGIBBS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mcgibbs_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_mcgibbs_debug_init(gen);
#endif

  /* make generators for marginal distributions */
  switch (gen->variant & MCGIBBS_VARMASK_VARIANT) {
  case MCGIBBS_VARIANT_COORDINATE:
    GEN->distr_condi = unur_distr_condi_new( gen->distr, GEN->state, NULL, 0);

    par_condi = unur_tdrgw_new(GEN->distr_condi);
    unur_set_use_distr_privatecopy( par_condi, FALSE );
    unur_set_debug( par_condi, (gen->debug&MCGIBBS_DEBUG_CONDI)?gen->debug:1u);

    GEN_CONDI[0] = gen_condi = unur_init(par_condi);
    /** TODO: error handling!!!! **/

    for (i=1; i<GEN->dim; i++) {
      GEN_CONDI[i] = unur_gen_clone(gen_condi);
    }

    break;

  case MCGIBBS_VARIANT_RANDOM_DIRECTION:
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_par_free(par); _unur_mcgibbs_free(gen);
    return NULL;
  }

  /* free parameters */
  _unur_par_free(par);

  /* o.k. */
  return gen;

} /* end of _unur_mcgibbs_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_mcgibbs_create( struct unur_par *par )
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
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MCGIBBS_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mcgibbs_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MCGIBBS_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim;

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mcgibbs_coord_sample_cvec;
  gen->destroy = _unur_mcgibbs_free;
  gen->clone = _unur_mcgibbs_clone;

  /* variant of sampling method */
  gen->variant = par->variant;        
  
  /* copy parameters into generator object */
  GEN->thinning = PAR->thinning;                /* thinning factor           */
  
  /* allocate memory for state */
  GEN->state = _unur_xmalloc( GEN->dim * sizeof(double));
  if (PAR->x0) {
    memcpy( GEN->state, PAR->x0, GEN->dim * sizeof(double));
  }
  else {
    for (i=0; i<GEN->dim; i++) 
      GEN->state[i] = 0.;
  }

  /* generator(s) for conditional distributions */
  GEN->distr_condi = NULL;
  GEN_CONDI = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen *) );
  for (i=0; i<GEN->dim; i++) GEN_CONDI[i] = NULL;

  /* allocate memory for random direction */
  GEN->direction = _unur_xmalloc( GEN->dim * sizeof(double));

  /* defaults */
  GEN->coord = (GEN->dim)-1;      /* current coordinate of GIBBS chain.
				     we want to start with coordinate 0. */

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_mcgibbs_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mcgibbs_clone( const struct unur_gen *gen )
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
#define CLONE         ((struct unur_mcgibbs_gen*)clone->datap)
#define CLONE_CONDI   clone->gen_aux_list     

  int i;
  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MCGIBBS_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy state */
  CLONE->state = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy( CLONE->state, GEN->state, GEN->dim * sizeof(double));
  
  /* copy generators for conditional distributions */
  if (GEN->distr_condi) CLONE->distr_condi = _unur_distr_clone( GEN->distr_condi );

  /* GEN_CONDI is cloned by _unur_generic_clone */
  /* however, these use a pointer to GEN->distr which must be updated, too */
  if (CLONE_CONDI) {
    for (i=0; i<GEN->dim; i++)
      if (CLONE_CONDI[i])
	CLONE_CONDI[i]->distr = CLONE->distr_condi;
  }

  /* allocate memory for random direction */
  CLONE->direction = _unur_xmalloc( GEN->dim * sizeof(double));

  return clone;

#undef CLONE
#undef CLONE_CONDI
} /* end of _unur_mcgibbs_clone() */

/*****************************************************************************/

void
_unur_mcgibbs_coord_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  double X;
  int thinning;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_MCGIBBS_GEN,RETURN_VOID);

  for (thinning = GEN->thinning; thinning > 0; --thinning) {

    /* update coordinate direction */
    GEN->coord = (GEN->coord + 1) % GEN->dim;
    
    /* update conditional distribution */
    unur_distr_condi_set_condition( GEN->distr_condi, GEN->state, NULL, GEN->coord);
    unur_tdrgw_reinit(GEN_CONDI[GEN->coord]);
    /** TODO: error handline **/
    
    /* sample from distribution */
    X = unur_sample_cont(GEN_CONDI[GEN->coord]);
    
    /* update state */
    GEN->state[GEN->coord] = X;
  }
  
  /* copy current state into given vector */
  memcpy(vec, GEN->state, GEN->dim * sizeof(double)); 

  return;

} /* end of _unur_mcgibbs_coord_sample_cvec() */

/*****************************************************************************/

void
_unur_mcgibbs_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_MCGIBBS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MCGIBBS_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free vectors */
  if (GEN->state) free (GEN->state);
  if (GEN->direction) free (GEN->direction);

  /* free conditional distribution object */
  if (GEN->distr_condi) _unur_distr_free (GEN->distr_condi);

  /* GEN_CONDI is freed by _unur_generic_free */

  _unur_generic_free(gen);

} /* end of _unur_mcgibbs_free() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_mcgibbs_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MCGIBBS_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = MCGIBBS (Markov Chain - GIBBS sampler)\n",gen->genid);
  fprintf(log,"%s: variant = ",gen->genid);
  switch (gen->variant & MCGIBBS_VARMASK_VARIANT) {
  case MCGIBBS_VARIANT_COORDINATE:
    fprintf(log,"coordinate sampling (original Gibbs sampler)  [default]\n"); break;
  case MCGIBBS_VARIANT_RANDOM_DIRECTION:
    fprintf(log,"random directions\n"); break;
  }
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  switch (gen->variant & MCGIBBS_VARMASK_VARIANT) {
  case MCGIBBS_VARIANT_COORDINATE:
    fprintf(log,"%s: sampling routine = _unur_mcgibbs_coord_sample()\n",gen->genid);
    break;
  case MCGIBBS_VARIANT_RANDOM_DIRECTION:
    fprintf(log,"%s: sampling routine = _unur_mcgibbs_randomdir_sample()\n",gen->genid);
    break;
  }

  _unur_matrix_print_vector( GEN->dim, GEN->state, "starting point = ", log, gen->genid, "\t   ");

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: INIT completed **********************\n",gen->genid);

} /* end of _unur_mcgibbs_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
