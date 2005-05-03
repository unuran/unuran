/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      hitrou.c                                                     *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    hit and run ratio-of-uniforms method                         *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and (optionally) a bounding rectangle for the acceptance   *
 *      region.                                                              *
 *      Produce a value x consistent with its density                        *
 *      The bounding rectangle is computed numerically if it is not given.   *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *   OPTIONAL:                                                               *
 *      mode of the density                                                  *
 *      bounding rectangle of acceptance region                              *
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
 *   [1] Wakefield J.C., Gelfand A.E., Smith A.F.M.                          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method.                                                             *
 *       Statistics and Computing (1991) 1, pp (129-133)                     *
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
#include <utils/fmax_source.h>
#include <utils/hooke_source.h>
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include <uniform/urng.h>
#include <parser/parser.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "x_gen.h"
#include "arou.h"
#include "hitrou.h"
#include "hitrou_struct.h"

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define HITROU_SET_U       0x001u     /* set u values of bounding rectangle   */
#define HITROU_SET_V       0x002u     /* set v values of bounding rectangle   */
#define HITROU_SET_SKIP    0x004u     /* set skip-parameter                   */
#define HITROU_SET_R       0x008u     /* set r-parameter                      */
#define HITROU_SET_RADIUS  0x010u     /* set ball radius                      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "HITROU"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_hitrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void  _unur_hitrou_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_hitrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitrou_clone( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following function print debugging information on output stream,      */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_hitrou_debug_init ( const struct unur_gen *gen );

#endif

static void _unur_hitrou_random_direction( struct unur_gen *gen,
                                           int dim, double *direction);
/*---------------------------------------------------------------------------*/
/* generte a random direction vector                                         */
/*---------------------------------------------------------------------------*/

static int _unur_hitrou_inside_shape( UNUR_GEN *gen, double *uv );
/*---------------------------------------------------------------------------*/
/* check if point uv[] is inside shape                                       */
/*---------------------------------------------------------------------------*/

static void _unur_hitrou_x_to_uv( UNUR_GEN *gen, double *x, double *uv );
static void _unur_hitrou_uv_to_x( UNUR_GEN *gen, double *uv, double *x );
/*---------------------------------------------------------------------------*/
/* transforming between x[] and uv[] coordinates                             */
/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_hitrou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_hitrou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

#define NORMAL    gen->gen_aux        /* pointer to normal variate generator */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */
#define HITROU_VARIANT_COORDINATE       0x0001u /* coordinate sampler        */
#define HITROU_VARIANT_RANDOM_DIRECTION 0x0002u /* random direction sampler  */
#define HITROU_VARIANT_BALL             0x0003u /* ball sampler              */

#define HITROU_INITIAL_POINTS 20 /* random directions used for initial point */

/* increase strip positions (vmax) by this factor */
#define HITROU_ADAPTIVE_MULTIPLIER 1.1 

/* do we need these limits ? */
#define HITROU_BALL_RADIUS_MIN 1e-5
#define HITROU_BALL_RADIUS_MAX 1e5

#define HITROU_BALL_RADIUS_FACTOR 1.1

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_hitrou_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_hitrou_par) );
  COOKIE_SET(par,CK_HITROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR->dim = distr->dim;

  /* set default values */
  PAR->r   = 1.;         /* r-parameter of the generalized method       */
  PAR->skip      = 0;    /* number of skipped points in chain           */
  PAR->vmax      = 0.;   /* v-boundary of bounding rectangle (unknown)  */
  PAR->umin  = NULL;     /* u-boundary of bounding rectangle (unknown)  */
  PAR->umax  = NULL;     /* u-boundary of bounding rectangle (unknown)  */
  PAR->bounding_rectangle = 0; /* should we calculate the bounding rect */
  PAR->adaptive_points = 1; /* reusing outside points as new line-segment ends */
  PAR->adaptive_strip = 0;  /* usage of adaptive algorithm for the strip position */
  PAR->adaptive_ball = 0;   /* usage of adaptive algorithm for the ball radius */
  PAR->ball_radius = 1.; /* initial ball radius for ball sampler (if not set)  */
  par->method   = UNUR_METH_HITROU;   /* method and default variant          */
  par->variant  = HITROU_VARIANT_RANDOM_DIRECTION; /* default variant        */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */
  
  /* routine for starting generator */
  par->init = _unur_hitrou_init;

  return par;

} /* end of unur_hitrou_new() */

/*****************************************************************************/

int
unur_hitrou_set_u( struct unur_par *par, double *umin, double *umax )
     /*----------------------------------------------------------------------*/
     /* Sets left and right u-boundary of bounding rectangle.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   umin ... left boundary of rectangle                                */
     /*   umax ... right boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );

  /* check new parameter for generator */
  dim = PAR->dim; /* making source code more readable */
  for (d=0; d<dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }

  /* set values */
  PAR->umin = umin;
  PAR->umax = umax;

  /* changelog */
  par->set |= HITROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_u() */

/*---------------------------------------------------------------------------*/

int
unur_hitrou_set_v( struct unur_par *par, double vmax )
     /*----------------------------------------------------------------------*/
     /* Sets upper v-boundary of bounding rectangle.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   vmax ... upper boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store values */
  PAR->vmax = vmax;

  /* changelog */
  par->set |= HITROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_v() */

/*---------------------------------------------------------------------------*/

int
unur_hitrou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* Set the r-parameter for the generalized ratio-of-uniforms method.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   r  ... r-parameter                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->r = r;

  /* changelog */
  par->set |= HITROU_SET_R;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_hitrou_set_ball_radius( struct unur_par *par, double ball_radius )
     /*----------------------------------------------------------------------*/
     /* Sets radius of ball used for the ball-sampler.                       */
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
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (ball_radius <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ball_radius <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store values */
  PAR->ball_radius = ball_radius;

  /* changelog */
  par->set |= HITROU_SET_RADIUS;

  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_ball_radius() */

/*---------------------------------------------------------------------------*/


int
unur_hitrou_set_skip( struct unur_par *par, long skip )
     /*----------------------------------------------------------------------*/
     /* Set the skip-parameter for the generalized ratio-of-uniforms method. */
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
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (skip < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"skip<0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->skip = skip;

  /* changelog */
  par->set |= HITROU_SET_SKIP;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_skip() */


/*****************************************************************************/

int
unur_hitrou_use_bounding_rectangle( struct unur_par *par, int flag )
     /*----------------------------------------------------------------------*/
     /* Set the flag for the calculating and usage of the bounding rectangle */
     /*                                                                      */
     /* parameters:                                                          */
     /*   flag ... 1:use bounding rectangle, 0:use infinite strip            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (flag != 0 && flag != 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid flag for bounding rect");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->bounding_rectangle = flag;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_use_bounding_rectangle() */

/*****************************************************************************/

int 
unur_hitrou_set_adaptive_points( struct unur_par *par, int adaptive_flag )
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (adaptive_flag!=0 && adaptive_flag!=1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"adaptive_flag must be 0 or 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->adaptive_points = adaptive_flag;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_adaptive_points() */

/*****************************************************************************/

int 
unur_hitrou_set_adaptive_strip( struct unur_par *par, int adaptive_flag )
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (adaptive_flag!=0 && adaptive_flag!=1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"adaptive_flag must be 0 or 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->adaptive_strip = adaptive_flag;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_adaptive_strip() */

/*****************************************************************************/

int 
unur_hitrou_set_adaptive_ball( struct unur_par *par, int adaptive_flag )
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* check new parameter for generator */
  if (adaptive_flag!=0 && adaptive_flag!=1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"adaptive_flag must be 0 or 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->adaptive_ball = adaptive_flag;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_adaptive_ball() */

/*****************************************************************************/

int 
unur_hitrou_set_variant_coordinate( UNUR_PAR *par ) {
     /*----------------------------------------------------------------------*/
     /* Coordinate Sampler :                                                 */
     /* Sampling along the coordinate directions (cyclic).                   */
     /*----------------------------------------------------------------------*/

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );
            
  par->variant  = HITROU_VARIANT_COORDINATE;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitrou_set_variant_coordinate() */

/*****************************************************************************/

int 
unur_hitrou_set_variant_random_direction( UNUR_PAR *par ) {
     /*----------------------------------------------------------------------*/
     /* Random Direction Sampler :                                           */
     /* Sampling along the random directions.                                */
     /* This is the default.                                                 */
     /*----------------------------------------------------------------------*/
  
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );
  
  par->variant  = HITROU_VARIANT_RANDOM_DIRECTION;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitrou_set_variant_random_direction() */


/*****************************************************************************/

int 
unur_hitrou_set_variant_ball( UNUR_PAR *par ) {
     /*----------------------------------------------------------------------*/
     /* Ball Sampler :                                                       */
     /* Sampling along random directions uniformly distributed inside ball   */
     /*----------------------------------------------------------------------*/
  
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );
  
  par->variant  = HITROU_VARIANT_BALL;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitrou_set_variant_ball() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_hitrou_init( struct unur_par *par )
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
  int d, i;
  double pdf, pdf_max;
  
  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_HITROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HITROU_PAR,NULL);

  if ( PAR->bounding_rectangle==0 && par->variant == HITROU_VARIANT_COORDINATE) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"strip + coordinate samlper not possible");
    return NULL; 
  }

  /* create a new empty generator object */
  gen = _unur_hitrou_create(par);
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

  if (GEN->bounding_rectangle==0 && GEN->adaptive_strip==1) {
    /* we are using an adaptive strip */
    
    /* do we have an vmax set ? */
    if (GEN->vmax==0) GEN->vmax=1.; /* initial value for strip position */

    /* we look for an interior point along different random directions */
    for (pdf_max=0, i=1; i<=HITROU_INITIAL_POINTS; i++) {
    
      /* choose a random direction */
      _unur_hitrou_random_direction(gen, GEN->dim+1, GEN->direction);
      _unur_hitrou_uv_to_x( gen, GEN->direction, GEN->x );

      /* check if this point is better (having higher value of pdf) */
      pdf = PDF(GEN->x);
      if (pdf_max < pdf ) {
        /* update current max of the pdf */
	pdf_max = pdf;
	
	/* point_random[] on the RoU-surface along the given direction */
        _unur_hitrou_x_to_uv( gen, GEN->x, GEN->point_random );
                
        /* get the middle point */
        for (d=0; d<=GEN->dim; d++)
          GEN->point_current[d] = GEN->point_random[d] / 2;	  
      }
      
    } /* next random direction */    
        
  }
  else {

    /* compute bounding rectangle */
    if (_unur_hitrou_rectangle(gen)!=UNUR_SUCCESS) {
      _unur_par_free(par); _unur_hitrou_free(gen);
      return NULL;
    }
    /* set initial point inside RoU shape */
    GEN->point_current[GEN->dim]=GEN->vmax/2.; /* all other coordinates are already 0 */      
  }
  
#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_hitrou_debug_init(gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_hitrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_hitrou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  int d; /* index used in dimension loops (0 <= d < dim) */
  struct MROU_RECTANGLE *rr;

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_HITROU_GEN, UNUR_ERR_COOKIE );

  /* Boundary rectangle is already set */
  if ((gen->set & HITROU_SET_U) && (gen->set & HITROU_SET_V)) {
    return UNUR_SUCCESS;
  }

  /* Allocating and filling mrou_rectangle struct */
  rr = _unur_mrou_rectangle_new();

  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->umin;
  rr->umax   = GEN->umax;
  rr->r      = GEN->r;
  rr->center = GEN->center;
  rr->bounding_rectangle = GEN->bounding_rectangle;
  rr->genid  = gen->genid;
  
  /* setting initialization values */
  for (d=0; d<GEN->dim; d++) {
    rr->umin[d]=0;
    rr->umax[d]=0;
  } 
  
  /* calculate bounding rectangle */
  if (!(gen->set & HITROU_SET_U) || !(gen->set & HITROU_SET_V)) {
    _unur_mrou_rectangle_compute(rr);
  }

  if (!(gen->set & HITROU_SET_V)) {
     /* user has not provided any upper bound for v */
     GEN->vmax = rr->vmax;
  }

  if (!(gen->set & HITROU_SET_U)) {
    /* user has not provided any bounds for u */
    for (d=0; d<GEN->dim; d++) {
      GEN->umin[d] = rr->umin[d];
      GEN->umax[d] = rr->umax[d];
    }
  }

  free(rr);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_hitrou_rectangle() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hitrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HITROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_hitrou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_HITROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_hitrou_sample_cvec;

  gen->destroy = _unur_hitrou_free;
  gen->clone = _unur_hitrou_clone;

  /* variant of sampling method */
  gen->variant = par->variant;        
    
  /* allocate memory for u-boundary arrays */
  GEN->umin = _unur_xmalloc( PAR->dim * sizeof(double)); /* bounding rectangle */
  GEN->umax = _unur_xmalloc( PAR->dim * sizeof(double)); /* bounding rectangle */

  /* allocate memory for random direction vector in (U,V) space */
  GEN->direction = _unur_xmalloc( (PAR->dim+1) * sizeof(double));

  /* allocate memory for current interior and random candidate point */
  GEN->point_current = _unur_xmalloc( (PAR->dim+1) * sizeof(double));
  GEN->point_random  = _unur_xmalloc( (PAR->dim+1) * sizeof(double));

  /* allocate memory for working point (in the (x,y)-coordinate system */
  GEN->x  = _unur_xmalloc( (PAR->dim) * sizeof(double));

  /* copy parameters into generator object */
  GEN->dim   = PAR->dim;              /* dimension */
  GEN->r     = PAR->r;                /* r-parameter of the hitrou method */
  GEN->vmax  = PAR->vmax;             /* upper v-boundary of bounding rectangle */
  GEN->skip  = PAR->skip;             /* number of skipped poins in chain */
  GEN->bounding_rectangle = PAR->bounding_rectangle; /* using bounding rect flag */
  GEN->adaptive_points = PAR->adaptive_points;  /* reusing outside points for line-segment */
  GEN->adaptive_strip = PAR->adaptive_strip;    /* using adaptive strip flag */
  GEN->adaptive_ball = PAR->adaptive_ball;    /* using adaptive ball flag */
  GEN->ball_radius = PAR->ball_radius; /* ball radius of ball sampler */
  
  if (PAR->umin != NULL) memcpy(GEN->umin, PAR->umin, GEN->dim * sizeof(double));
  if (PAR->umax != NULL) memcpy(GEN->umax, PAR->umax, GEN->dim * sizeof(double));

  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);

  /* initialize parameters */
  GEN->pdfcount = 0;
  GEN->coordinate = 0;  
  for (d=0; d<GEN->dim+1; d++) {
    GEN->point_current[d]=0.;
    GEN->point_random[d]=0.;
  }

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_hitrou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hitrou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_hitrou_gen*)clone->datap)

  struct unur_gen *clone;


  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HITROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory for arrays */
  CLONE->umin = _unur_xmalloc( GEN->dim * sizeof(double));
  CLONE->umax = _unur_xmalloc( GEN->dim * sizeof(double));

  CLONE->point_current = _unur_xmalloc( (GEN->dim+1) * sizeof(double));
  CLONE->point_random  = _unur_xmalloc( (GEN->dim+1) * sizeof(double));
  CLONE->direction = _unur_xmalloc( (GEN->dim+1) * sizeof(double));

  /* copy parameters into clone object */
  CLONE->skip = GEN->skip;
  CLONE->r = GEN->r;
  CLONE->vmax = GEN->vmax;
  CLONE->bounding_rectangle = GEN->bounding_rectangle;
  CLONE->adaptive_points = GEN->adaptive_points;
  CLONE->adaptive_strip = GEN->adaptive_strip;
  CLONE->adaptive_ball = GEN->adaptive_ball;
  CLONE->ball_radius = GEN->ball_radius;
  
  memcpy(CLONE->umin, GEN->umin, GEN->dim * sizeof(double));
  memcpy(CLONE->umax, GEN->umax, GEN->dim * sizeof(double));
  memcpy(CLONE->point_current, GEN->point_current, (GEN->dim+1) * sizeof(double));
  memcpy(CLONE->point_random , GEN->point_random , (GEN->dim+1) * sizeof(double));
  memcpy(CLONE->direction, GEN->direction, (GEN->dim+1) * sizeof(double));
  memcpy(CLONE->x, GEN->x, GEN->dim * sizeof(double));

  /* copy data */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  return clone;

#undef CLONE
} /* end of _unur_hitrou_clone() */

/*****************************************************************************/

void
_unur_hitrou_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double lambda, lmin, lmax; /* lambda parameters of line */
  long skip;
  int inside;
  double u;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);
  COOKIE_CHECK(gen,CK_HITROU_GEN,RETURN_VOID);

  dim = GEN->dim;

  for (skip=0; skip<=GEN->skip; skip++) {
 
    /* initialization */
    lmin=0; lmax=0;
   
    if (gen->variant == HITROU_VARIANT_RANDOM_DIRECTION) {
      /* generate random direction vector in (U,V) space */
      _unur_hitrou_random_direction(gen, dim+1, GEN->direction);
    
      if (GEN->adaptive_strip==1) {
        /* adaptive step for the upper strip position vmax */
        inside=1; 
        while (inside==1) {
          /* calculate lambda parameters for the intersections with v=vmax */
          lambda = (GEN->vmax - GEN->point_current[dim]) / GEN->direction[dim];
          if (lambda>0) lmax = lambda;
          if (lambda<0) lmin = lambda;
      
          /* calculate the "candidate" point along the given random direction */
          for (d=0; d<=dim; d++)
            GEN->point_random[d] = GEN->point_current[d] + lambda * GEN->direction[d];
      
          inside = _unur_hitrou_inside_shape(gen, GEN->point_random);

          /* increase vmax if necessary */
          if (inside==1) GEN->vmax *= HITROU_ADAPTIVE_MULTIPLIER;    
        }
      }
      
      /* calculate lambda parameters for the intersections with v=vmax and v=0 */
      lambda = (GEN->vmax - GEN->point_current[dim]) / GEN->direction[dim];
      if (lambda>0) lmax = lambda;
      if (lambda<0) lmin = lambda;

      lambda = (0 - GEN->point_current[dim]) / GEN->direction[dim];
      if (lambda>0) lmax = lambda;
      if (lambda<0) lmin = lambda;
      
      if (GEN->bounding_rectangle==1) {
        /* calculate the intersections alont all other coordinate directions */
        for (d=0; d<dim; d++) {
          lambda = (GEN->umin[d] - GEN->point_current[d]) / GEN->direction[d];
          if (lambda>0 && lambda<lmax) lmax = lambda;
          if (lambda<0 && lambda>lmin) lmin = lambda;

          lambda = (GEN->umax[d] - GEN->point_current[d]) / GEN->direction[d];
          if (lambda>0 && lambda<lmax) lmax = lambda;
          if (lambda<0 && lambda>lmin) lmin = lambda;
        }
      }
    }
    
    if (gen->variant == HITROU_VARIANT_COORDINATE) {
        
	/* calculate the intersections alont the current coordinate direction */
	
	/* setting current direction vector */
	for (d=0; d<=dim; d++) 
	  GEN->direction[d]=0.;
	
	d = GEN->coordinate;
        GEN->direction[d]=1.;
	
	if (d==dim) {
	  /* calculate lambda parameters for the intersections with v=vmax and v=0 */
          lambda = (GEN->vmax - GEN->point_current[dim]) / GEN->direction[dim];
          if (lambda>0) lmax = lambda;
          if (lambda<0) lmin = lambda;

          lambda = (0 - GEN->point_current[dim]) / GEN->direction[dim];
          if (lambda>0) lmax = lambda;
          if (lambda<0) lmin = lambda;

	}
	
	else {
	  /* calculate lambda parameters for the intersections with u[d]-planes */
	  lambda = (GEN->umin[d] - GEN->point_current[d]) / GEN->direction[d];
          if (lambda>0) lmax = lambda;
          if (lambda<0) lmin = lambda;

          lambda = (GEN->umax[d] - GEN->point_current[d]) / GEN->direction[d];
          if (lambda>0) lmax = lambda;
          if (lambda<0) lmin = lambda;
        }
    
    }
        
    /* until we find an inside point */
    while (1) {

      u = _unur_call_urng(gen->urng);
      lambda = lmin + (lmax-lmin) * u;
      if (GEN->adaptive_points==1) {
        if (lambda>0) lmax=lambda;
        if (lambda<0) lmin=lambda;
      }
      
      if (gen->variant == HITROU_VARIANT_BALL) {
        _unur_hitrou_random_direction(gen, dim+1, GEN->direction);
        lambda = GEN->ball_radius * pow(u, 1./(dim+1));
      }
      
      /* calculate the "candidate" point along the given random direction */
      for (d=0; d<=dim; d++)
        GEN->point_random[d] = GEN->point_current[d] + lambda * GEN->direction[d];
	
      /* check if random point is inside domain */
      if (_unur_hitrou_inside_shape(gen, GEN->point_random)) {
        /* update current point */
        for (d=0; d<=dim; d++)
          GEN->point_current[d] = GEN->point_random[d] ;
        
	if (gen->variant == HITROU_VARIANT_BALL && 
	    GEN->adaptive_ball==1 && 
	    GEN->ball_radius*HITROU_BALL_RADIUS_FACTOR < HITROU_BALL_RADIUS_MAX) 
	      GEN->ball_radius *= HITROU_BALL_RADIUS_FACTOR; 
	       
	break; /* jump out of the while() loop */
      }

      else {
        /* we are outside shape */
        if (gen->variant == HITROU_VARIANT_BALL) {
          /* no change of current point : returning the current point */
  	  if (GEN->adaptive_ball==1 && GEN->ball_radius/HITROU_BALL_RADIUS_FACTOR > HITROU_BALL_RADIUS_MIN) 
	    GEN->ball_radius /= HITROU_BALL_RADIUS_FACTOR;  
	  break;
        }
      }
    
    }

    if (gen->variant == HITROU_VARIANT_COORDINATE) {
      /* prepare next coordinate direction */
      GEN->coordinate = GEN->coordinate + 1; 
      if (GEN->coordinate > dim) GEN->coordinate = 0;
    }
    
  }

  /* calculate the sample point in the X[]-coordinate system            */
  _unur_hitrou_uv_to_x( gen, GEN->point_current, vec );
  
#if 0
  _unur_matrix_print_vector ( dim+1, GEN->point_current, "uv :", stdout, "", "---" );
  _unur_matrix_print_vector ( dim, vec, "x :", stdout, "", "---" );
  printf("GEN->ball_radius = %f\n", GEN->ball_radius);
#endif
  
  return;
} /* end of _unur_hitrou_sample() */


/*****************************************************************************/

void
_unur_hitrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_HITROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HITROU_GEN,RETURN_VOID);


#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) {} /* write additional info into log file */;
#endif

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->umin) free(GEN->umin);
  if (GEN->umax) free(GEN->umax);
  if (GEN->direction) free(GEN->direction);
  if (GEN->point_current) free(GEN->point_current);
  if (GEN->point_random)  free(GEN->point_random);
  if (GEN->x)  free(GEN->x);
  _unur_generic_free(gen);

} /* end of _unur_hitrou_free() */

/*****************************************************************************/
/**  Auxilliary routines                                                    **/
/*****************************************************************************/


void 
_unur_hitrou_x_to_uv( UNUR_GEN *gen, double *x, double *uv ) {
    /* transform the point X[] into the UV[]-coordinate system */
    int d;
    double V;
    
    V = pow(PDF(x), 1./(GEN->r * GEN->dim + 1.));
    
    uv[GEN->dim] = V;
    for (d=0; d<GEN->dim; d++) {
      if (GEN->r==1)
        uv[d] = (x[d] - GEN->center[d]) * V;
      else
        uv[d] = (x[d] - GEN->center[d]) * pow(V,GEN->r) ;
    } 
}
/*---------------------------------------------------------------------------*/

void 
_unur_hitrou_uv_to_x( UNUR_GEN *gen, double *uv, double *x ) {
    /* transform the point UV[] into the X[]-coordinate system */
    int d;
    double U,V;
    
    V = uv[GEN->dim];
    
    if (V==0) {
      for (d=0; d<GEN->dim; d++) 
        x[d]=0;
    }
    
    else {
      for (d=0; d<GEN->dim; d++) {
        U = uv[d];
        if (GEN->r==1)
          x[d] = U/V + GEN->center[d];
        else
          x[d] = U/pow(V,GEN->r) + GEN->center[d];
      }
    }
}

/*---------------------------------------------------------------------------*/

void
_unur_hitrou_random_direction( struct unur_gen *gen,
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

int
_unur_hitrou_inside_shape( UNUR_GEN *gen, double *uv )
     /* check if uv[] is inside RoU shape */
{
  double V;
  int inside=0;

  _unur_hitrou_uv_to_x( gen, uv, GEN->x );

  /* point inside domain ? */
  V = uv[GEN->dim];
  if (V>0 && V < pow(PDF(GEN->x),1./(GEN->r * GEN->dim + 1.)))
    inside=1;
  else
    inside=0;

  GEN->pdfcount++;
  
  return inside;
}

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Additional routines used for testing                                   **/
/*****************************************************************************/

long
_unur_hitrou_get_pdfcount( UNUR_GEN *gen)
     /* Return the number of PDF calls */
{
  return GEN->pdfcount;
}

/*---------------------------------------------------------------------------*/

void
_unur_hitrou_reset_pdfcount( UNUR_GEN *gen)
     /* Reset the number of PDF calls to 0 */
{
  GEN->pdfcount = 0;
}


/*---------------------------------------------------------------------------*/

void _unur_hitrou_set_point( UNUR_GEN *gen, double *uv)
     /* set the current point (dimension=dim+1) inside the testrectangle */
{
  int d;

  for (d=0; d<=GEN->dim; d++) {
    GEN->point_current[d]=uv[d];
  }
}

/*---------------------------------------------------------------------------*/

void _unur_hitrou_get_point( UNUR_GEN *gen, double *uv)
     /* get the current point (dimension=dim+1) inside the testrectangle */
{
  int d;

  for (d=0; d<=GEN->dim; d++) {
    uv[d]=GEN->point_current[d];
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
_unur_hitrou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double vol;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITROU_GEN,RETURN_VOID);

  log = unur_get_stream();
  dim = GEN->dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = hitrou (hit and run ratio-of-uniforms)\n",gen->genid);
  if (gen->variant == HITROU_VARIANT_COORDINATE) 
    fprintf(log,"%s: variant = coordinate sampler \n",gen->genid);
  if (gen->variant == HITROU_VARIANT_RANDOM_DIRECTION) 
    fprintf(log,"%s: variant = random direction sampler \n",gen->genid);
  if (gen->variant == HITROU_VARIANT_BALL) 
    fprintf(log,"%s: variant = ball sampler \n",gen->genid);  
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_hitrou_sample",gen->genid);

  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(log,"%s: r-parameter = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,HITROU_SET_R);
  fprintf(log,"\n%s:\n",gen->genid);

  /* skip */
  fprintf(log,"%s: skip = %ld",gen->genid, GEN->skip);
  _unur_print_if_default(gen,HITROU_SET_SKIP);
  fprintf(log,"\n%s:\n",gen->genid);
  
  /* print center */
  _unur_matrix_print_vector( GEN->dim, GEN->center, "center =", log, gen->genid, "\t   ");

  /* print adaptive flag for points */
  fprintf(log,"%s: adaptive_points:", gen->genid);  
  if (GEN->adaptive_points==0)
    fprintf(log,"\tno (direction line segment is kept constant)");  
  else
    fprintf(log,"\tyes (direction line segment is adapted by each step)");
  fprintf(log,"\n");

  /* print adaptive flag for strip */
  fprintf(log,"%s: adaptive_strip:", gen->genid);  
  if (GEN->adaptive_strip==0)
    fprintf(log,"\tno (position is kept constant at vmax)");  
  else
    fprintf(log,"\tyes (position is adapted by each step)");
  fprintf(log,"\n");
  
  /* print bounding rectangle */
  fprintf(log,"%s: Rectangle:",gen->genid);
  if (!((gen->set & HITROU_SET_U) && (gen->set & HITROU_SET_V)))
    fprintf(log,"\t[computed]");
  else
    fprintf(log,"\t[input]");
  fprintf(log,"\n");

  vol = GEN->vmax;
  fprintf(log,"%s:\tvmax = %g\n",gen->genid, GEN->vmax);

  if (GEN->bounding_rectangle==1) {
    for (d=0; d<dim; d++) {
      vol *= (GEN->umax[d]-GEN->umin[d]);
      fprintf(log,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid,
        d, d, GEN->umin[d], GEN->umax[d]);
    }
    fprintf(log,"%s:\n",gen->genid);
    fprintf(log,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN->r*GEN->dim+1));
  }
  else {
    fprintf(log,"%s:\tumin[],umax[] are not used\n" ,gen->genid);
  }

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_hitrou_debug_init() */


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
