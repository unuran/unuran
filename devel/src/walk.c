/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      walk.c                                                       *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    random walk sampler                                          *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and (optionally) a bounding rectangle for the acceptance   *
 *      region.                                                              *
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
 *****************************************************************************
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <urng/urng.h>
#include <parser/parser.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "x_gen.h"
#include "arou.h"
#include "walk.h"
#include "walk_struct.h"

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define WALK_SET_THINNING  0x004u     /* set thinning-parameter              */
#define WALK_SET_RADIUS    0x010u     /* set ball radius                     */

/*---------------------------------------------------------------------------*/

#define GENTYPE "WALK"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_walk_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_walk_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_walk_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_walk_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_walk_clone( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following function print debugging information on output stream,      */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_walk_debug_init ( const struct unur_gen *gen );

#endif

static void _unur_walk_random_direction( struct unur_gen *gen,
                                           int dim, double *direction);
/*---------------------------------------------------------------------------*/
/* generte a random direction vector                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_walk_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_walk_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

#define NORMAL    gen->gen_aux        /* pointer to normal variate generator */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */


/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_walk_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_walk_par) );
  COOKIE_SET(par,CK_WALK_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR->dim = distr->dim;

  /* set default values */
  PAR->thinning    = 1;  /* thinning parameter of chain */
  PAR->ball_radius = 1.; /* ball radius (if not set)  */
  par->method   = UNUR_METH_WALK;   /* method and default variant          */
  par->variant  = 0;                  /* default variant        */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */
  
  /* routine for starting generator */
  par->init = _unur_walk_init;

  return par;

} /* end of unur_walk_new() */

/*****************************************************************************/


int
unur_walk_set_ball_radius( UNUR_PAR *par, double ball_radius )
     /*----------------------------------------------------------------------*/
     /* Sets radius of ball used for the random walk sampler.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ball_radius                                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, WALK );

  /* check new parameter for generator */
  if (ball_radius <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ball_radius <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store values */
  PAR->ball_radius = ball_radius;

  /* changelog */
  par->set |= WALK_SET_RADIUS;

  return UNUR_SUCCESS;

} /* end of unur_walk_set_ball_radius() */

/*****************************************************************************/

int
unur_walk_set_thinning( UNUR_PAR *par, long thinning )
     /*----------------------------------------------------------------------*/
     /* Set the thinning-parameter                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   thinning     ... thinning parameter                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, WALK );

  /* check new parameter for generator */
  if (thinning < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning<1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->thinning = thinning;

  /* changelog */
  par->set |= WALK_SET_THINNING;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_walk_set_thinning() */

/*****************************************************************************/



/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_walk_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_WALK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_WALK_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_walk_create(par);
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

  
#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_walk_debug_init(gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_walk_init() */

/*---------------------------------------------------------------------------*/


struct unur_gen *
_unur_walk_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_WALK_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_walk_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_WALK_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_walk_sample_cvec;

  gen->destroy = _unur_walk_free;
  gen->clone = _unur_walk_clone;

  /* variant of sampling method */
  gen->variant = par->variant;        
    

  /* allocate memory for random direction vector */
  GEN->direction = _unur_xmalloc( (PAR->dim) * sizeof(double));

  /* allocate memory for current interior and random candidate point */
  GEN->point_current = _unur_xmalloc( (PAR->dim) * sizeof(double));
  GEN->point_random  = _unur_xmalloc( (PAR->dim) * sizeof(double));

  /* allocate memory for working point (in the (x,y)-coordinate system */
  GEN->x  = _unur_xmalloc( (PAR->dim) * sizeof(double));

  /* copy parameters into generator object */
  GEN->dim   = PAR->dim;               /* dimension */
  GEN->thinning  = PAR->thinning;      /* thinning parameter of chain */
  GEN->ball_radius = PAR->ball_radius; /* ball radius */
  
  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);

  /* initialize parameters */
  for (d=0; d<GEN->dim; d++) {
    GEN->point_current[d]=0.;
    GEN->point_random[d]=0.;
  }

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_walk_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_walk_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_walk_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_WALK_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  CLONE->point_current = _unur_xmalloc( (GEN->dim) * sizeof(double));
  CLONE->point_random  = _unur_xmalloc( (GEN->dim) * sizeof(double));
  CLONE->direction     = _unur_xmalloc( (GEN->dim) * sizeof(double));
  CLONE->x             = _unur_xmalloc( (GEN->dim) * sizeof(double));

  memcpy(CLONE->point_current, GEN->point_current, (GEN->dim) * sizeof(double));
  memcpy(CLONE->point_random , GEN->point_random , (GEN->dim) * sizeof(double));
  memcpy(CLONE->direction,     GEN->direction,     (GEN->dim) * sizeof(double));
  memcpy(CLONE->x,             GEN->x,             (GEN->dim) * sizeof(double));

  /* copy data */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  return clone;

#undef CLONE
} /* end of _unur_walk_clone() */

/*****************************************************************************/

int
_unur_walk_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double lambda; 
  long thinning;
  double u;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_WALK_GEN,UNUR_ERR_COOKIE);

  dim = GEN->dim;

  for (thinning=1; thinning<=GEN->thinning; thinning++) {
 
    /* initialization */
    lambda=0;
   
    /* until we find an inside point */
    while (1) {

      u = _unur_call_urng(gen->urng);
      
      _unur_walk_random_direction(gen, dim, GEN->direction);
      lambda = GEN->ball_radius * pow(u, 1./dim);
      
      /* calculate the "candidate" point along the given random direction */
      for (d=0; d<dim; d++)
        GEN->point_random[d] = GEN->point_current[d] + lambda * GEN->direction[d];
            
      /* check if random point is to be accepted */
      u = _unur_call_urng(gen->urng);
      if ( u * PDF(GEN->point_current) <= PDF(GEN->point_random) ) {
        /* update current point */
        for (d=0; d<=dim; d++)
          GEN->point_current[d] = GEN->point_random[d] ;
        
	break; /* jump out of the while() loop */
      }
      else {
        /* we do not accept the random point  */
        /* no change of current point : returning the current point */
        break;
      }  

    } /* while() loop */
    
  } /* next thinning */

  
    /* returning the current sample point */
    for (d=0; d<dim; d++) 
      vec[d] = GEN->point_current[d];
  
      
#if 0
  printf("GEN->ball_radius = %f\n", GEN->ball_radius);
#endif
  
  return UNUR_SUCCESS;
} /* end of _unur_walk_sample() */


/*****************************************************************************/

void
_unur_walk_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_WALK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_WALK_GEN,RETURN_VOID);


#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) {} /* write additional info into log file */;
#endif

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->direction) free(GEN->direction);
  if (GEN->point_current) free(GEN->point_current);
  if (GEN->point_random)  free(GEN->point_random);
  if (GEN->x)  free(GEN->x);
  _unur_generic_free(gen);

} /* end of _unur_walk_free() */


/*****************************************************************************/
/**  Auxilliary routines                                                    **/
/*****************************************************************************/


/*---------------------------------------------------------------------------*/

void
_unur_walk_random_direction( struct unur_gen *gen,
                               int dim, double *direction)
     /*----------------------------------------------------------------------*/
     /* generte a random direction vector  (unit vector)                     */
     /*----------------------------------------------------------------------*/
{
  int d;

  for (d=0; d<dim; d++) {
    do {
      direction[d] = unur_sample_cont(NORMAL);
    } while (direction[d]==0.); /* extremely seldom case */
  }

  _unur_vector_normalize(dim, direction);
}

/*---------------------------------------------------------------------------*/

  

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Additional routines used for testing                                   **/
/*****************************************************************************/


/*---------------------------------------------------------------------------*/

void _unur_walk_set_point_current( UNUR_GEN *gen, double *x)
     /* set the current point (dimension=dim) */
{
  int d;

  for (d=0; d<GEN->dim; d++) {
    GEN->point_current[d]=x[d];
  }
}

/*---------------------------------------------------------------------------*/

void _unur_walk_get_point_current( UNUR_GEN *gen, double *x)
     /* get the current point (dimension=dim) */
{
  int d;

  for (d=0; d<GEN->dim; d++) {
    x[d]=GEN->point_current[d];
  }
}

/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_walk_debug_init( const struct unur_gen *gen )
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
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_WALK_GEN,RETURN_VOID);

  log = unur_get_stream();
  dim = GEN->dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = random walk sampler\n",gen->genid);

  /* ball radius */
  fprintf(log,"%s: ball radius = %g",gen->genid, GEN->ball_radius);
  _unur_print_if_default(gen,WALK_SET_RADIUS);
  fprintf(log,"\n%s:\n",gen->genid);
  
  /* thinning */
  fprintf(log,"%s: thinning = %ld",gen->genid, GEN->thinning);
  _unur_print_if_default(gen,WALK_SET_THINNING);
  fprintf(log,"\n%s:\n",gen->genid);
  
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_walk_sample",gen->genid);
  fprintf(log,"()\n%s:\n",gen->genid);

  /* print center */
  _unur_matrix_print_vector( GEN->dim, GEN->center, "center =", log, gen->genid, "\t   ");
  
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_walk_debug_init() */


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
