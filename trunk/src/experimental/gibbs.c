/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gibbs.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    Gibbs sampler using the full-conditionals.                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF                                                            *
 *      Produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *   OPTIONAL:                                                               *
 *      mode of the density                                                  *
 *                                                                           *
 *****************************************************************************
 *   (c) 2000 Wolfgang Hoermann and Josef Leydold                            *
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
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distr/cont.h>
#include <distr/condi.h>
#include <distributions/unur_distributions.h>
#include <utils/matrix_source.h>
#include <uniform/urng.h>
#include <methods/tdr.h>
#include <methods/arou.h>
#include <methods/x_gen_source.h>
#include <methods/x_gen.h>
#include "unur_methods_source.h"

#include "gibbs.h"
#include "gibbs_struct.h"

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define GIBBS_SET_SKIP    0x004u     /* set skip-parameter                   */

/*---------------------------------------------------------------------------*/

#define GENTYPE "GIBBS"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void  _unur_gibbs_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_clone( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_gibbs_debug_init ( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

static void _unur_gibbs_random_unit_vector( struct unur_gen *gen,
                                          int dim, double *direction);
/*---------------------------------------------------------------------------*/
/* generte a random direction vector                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_gibbs_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_gibbs_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

#define NORMAL    gen->gen_aux        /* pointer to normal variate generator */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */
#define GIBBS_VARIANT_COORDINATE       0x0001u /* coordinate sampler(default)*/
#define GIBBS_VARIANT_RANDOM_DIRECTION 0x0002u /* random direction sampler   */

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_gibbs_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_gibbs_par) );
  COOKIE_SET(par,CK_GIBBS_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR->dim = distr->dim;

  /* set default values */
  PAR->skip      = 0;         /* number of skipped points in chain           */
  par->method   = UNUR_METH_GIBBS;    /* method and default variant          */
  par->variant  = GIBBS_VARIANT_COORDINATE;  /* default variant              */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_gibbs_init;

  return par;

} /* end of unur_gibbs_new() */

/*****************************************************************************/

int
unur_gibbs_set_skip( struct unur_par *par, long skip )
     /*----------------------------------------------------------------------*/
     /* Set the skip-parameter for the Gibbs sampler.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   skip  ... skip-parameter                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* check new parameter for generator */
  if (skip < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"skip<0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->skip = skip;

  /* changelog */
  par->set |= GIBBS_SET_SKIP;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_gibbs_set_skip() */

/*****************************************************************************/

int unur_gibbs_set_variant_coordinate( UNUR_PAR *par ) {
     /*----------------------------------------------------------------------*/
     /* Coordinate Sampler :                                                 */
     /* Sampling along the coordinate directions (cyclic).                   */
     /* This is the default.                                                 */
     /*----------------------------------------------------------------------*/

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
            
  par->variant  = GIBBS_VARIANT_COORDINATE;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_gibbs_set_variant_coordinate() */

/*****************************************************************************/

int unur_gibbs_set_variant_random_direction( UNUR_PAR *par ) {
     /*----------------------------------------------------------------------*/
     /* Random Direction Sampler :                                           */
     /* Sampling along the random directions.                                */
     /*----------------------------------------------------------------------*/
  
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
  
  par->variant  = GIBBS_VARIANT_RANDOM_DIRECTION;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_gibbs_set_variant_random_direction() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_gibbs_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_GIBBS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_gibbs_create(par);
  if (!gen) {
    _unur_par_free(par);
    return NULL;
  }

  /* we need a generator for standard normal distributons */
  if (NORMAL==NULL) {
    struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
    struct unur_par   *normalpar = unur_arou_new( normaldistr );
    unur_arou_set_usedars( normalpar, TRUE );
    NORMAL = unur_init( normalpar );
    _unur_distr_free( normaldistr );
    if (NORMAL == NULL) {
       _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
        "Cannot create aux Gaussian generator");
       _unur_free(gen); free (par);
       return NULL;
    }
    /* uniform random number generator and debugging flags */
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
  }  
  
  /* set initial point */
  /* all coordinates are already set to 0 by _unur_gibbs_create() */
  /* any coordinate differing from 0 can be set here ... */

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_gibbs_debug_init(gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_gibbs_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_gibbs_create( struct unur_par *par )
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
  int d;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_gibbs_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_GIBBS_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_gibbs_sample_cvec;

  gen->destroy = _unur_gibbs_free;
  gen->clone = _unur_gibbs_clone;

  /* variant of sampling method */
  gen->variant = par->variant;        
  
  /* allocate memory for current point and random direction */
  GEN->point_current = _unur_xmalloc( (PAR->dim) * sizeof(double));
  GEN->direction = _unur_xmalloc( (PAR->dim) * sizeof(double));

  /* copy parameters into generator object */
  GEN->dim   = PAR->dim;              /* dimension */
  GEN->skip  = PAR->skip;             /* number of skipped poins in chain */
  
  /* initialize parameters */
  GEN->pdfcount = 0;
  GEN->coordinate = 0;
  for (d=0; d<GEN->dim; d++) {
    GEN->point_current[d]=0.;
  }

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_gibbs_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_gibbs_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_gibbs_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_GIBBS_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  CLONE->point_current = _unur_xmalloc( (GEN->dim) * sizeof(double));
  CLONE->direction = _unur_xmalloc( (GEN->dim) * sizeof(double));

  /* copy parameters into clone object */
  CLONE->skip = GEN->skip;

  memcpy(CLONE->point_current, GEN->point_current, (GEN->dim) * sizeof(double));
  memcpy(CLONE->direction, GEN->direction, (GEN->dim) * sizeof(double));
  
  return clone;

#undef CLONE
} /* end of _unur_gibbs_clone() */

/*****************************************************************************/

void
_unur_gibbs_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  int dim; 
  long skip;
  
  /* possibly move this to the gibbs generator structure */
  UNUR_PAR *par_conditional = NULL;
  UNUR_DISTR *distr_conditional = NULL;
  UNUR_GEN *gen_conditional = NULL;
    
  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  dim = GEN->dim;
 
  for (skip=0; skip<=GEN->skip; skip++) {

      if ( gen->variant == GIBBS_VARIANT_COORDINATE ) {
          distr_conditional = unur_distr_condi_new(gen->distr, GEN->point_current, NULL, GEN->coordinate);
      }
      
      if ( gen->variant == GIBBS_VARIANT_RANDOM_DIRECTION ) {
          _unur_gibbs_random_unit_vector(gen, dim, GEN->direction);
          distr_conditional = unur_distr_condi_new(gen->distr, GEN->point_current, GEN->direction, GEN->coordinate);
      }
      
      par_conditional = unur_tdr_new(distr_conditional);
      gen_conditional = unur_init(par_conditional);
      
      if (gen_conditional==NULL) {
         _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,
          "Cannot create aux conditional generator");
      }
      else {
        GEN->point_current[GEN->coordinate] = unur_sample_cont(gen_conditional);
      }
    
      /* free allocated memory */
      if (distr_conditional) unur_distr_free(distr_conditional);
      if (gen_conditional)   unur_free(gen_conditional);

      /* prepare next coordinate direction */
      GEN->coordinate = GEN->coordinate + 1; 
      if (GEN->coordinate >= dim) GEN->coordinate = 0;
      
  }

  /* copy current point coordinates */
  memcpy(vec, GEN->point_current, GEN->dim*sizeof(double)); 

  return;
} /* end of _unur_gibbs_sample() */


/*****************************************************************************/

void
_unur_gibbs_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_GIBBS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->point_current) free(GEN->point_current);
  if (GEN->direction) free(GEN->direction);
  _unur_generic_free(gen);

} /* end of _unur_gibbs_free() */

/*****************************************************************************/
/**  Auxilliary routines                                                    **/
/*****************************************************************************/




/*****************************************************************************/
/**  Additional routines used for testing                                   **/
/*****************************************************************************/

void
_unur_gibbs_random_unit_vector( struct unur_gen *gen,
                               int dim, double *direction)
     /*--------------------------------------------*/
     /* generte a random unit direction vector     */
     /*--------------------------------------------*/
{
  int d;

  for (d=0; d<dim; d++) {
    do {
      direction[d] = unur_sample_cont(NORMAL);
    } while (direction[d]==0.); /* extremely seldom case */
  }

  /* normalize direction vector */
  
}

/*---------------------------------------------------------------------------*/



/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_gibbs_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int dim; 

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  log = unur_get_stream();
  dim = GEN->dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = gibbs ",gen->genid);
  if (gen->variant==GIBBS_VARIANT_COORDINATE) 
    fprintf(log," (coordinate sampler)\n");
  if (gen->variant==GIBBS_VARIANT_RANDOM_DIRECTION) 
    fprintf(log," (random direction sampler)\n");
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_gibbs_sample",gen->genid);

  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */

  /* skip */
  fprintf(log,"%s: skip = %ld",gen->genid, GEN->skip);
  _unur_print_if_default(gen,GIBBS_SET_SKIP);
  fprintf(log,"\n%s:\n",gen->genid);
  
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_gibbs_debug_init() */

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
