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

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define HITROU_VARFLAG_VERIFY   0x002u   /* run verify mode                   */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define HITROU_SET_U       0x001u     /* set u values of bounding rectangle   */
#define HITROU_SET_V       0x002u     /* set v values of bounding rectangle   */
#define HITROU_SET_SKIP    0x004u     /* set skip-parameter                   */
#define HITROU_SET_R       0x008u     /* set r-parameter                      */

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
static void  _unur_hitrou_sample_check( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_hitrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitrou_clone( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* clone of generator object.                                                */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_hitrou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

static void _unur_hitrou_random_direction( struct unur_gen *gen, 
                                           int dim, double *direction);
/*---------------------------------------------------------------------------*/
/* generte a random direction vector                                         */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       par->data.hitrou       /* data for parameter object        */
#define GEN       gen->data.hitrou       /* data for generator object        */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */     
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

#define NORMAL    gen->gen_aux        /* pointer to normal variate generator */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
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
  par = _unur_xmalloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_HITROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR.dim = distr->dim;


  /* set default values */
  PAR.r		= 1.; 	      /* r-parameter of the generalized method       */
  PAR.skip      = 10;         /* number of skipped points in chain           */
  PAR.vmax      = 0.;         /* v-boundary of bounding rectangle (unknown)  */
  PAR.umin 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  PAR.umax 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  par->method   = UNUR_METH_HITROU;   /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
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
  dim = PAR.dim; /* making source code more readable */ 
  for (d=0; d<dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }
  
  /* set values */
  PAR.umin = umin;
  PAR.umax = umax;
  
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
  PAR.vmax = vmax;

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
  PAR.r = r;

  /* changelog */
  par->set |= HITROU_SET_R;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_r() */
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
  PAR.skip = skip;

  /* changelog */
  par->set |= HITROU_SET_SKIP;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_hitrou_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | HITROU_VARFLAG_VERIFY) : (par->variant & (~HITROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_hitrou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, HITROU, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify bounding rectangle on */
    gen->variant |= HITROU_VARFLAG_VERIFY;
    SAMPLE = _unur_hitrou_sample_check;
  }
  else {
    /* turn verify bounding rectangle off */
    gen->variant &= ~HITROU_VARFLAG_VERIFY;
    SAMPLE = _unur_hitrou_sample_cvec;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitrou_chg_verify() */

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

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_HITROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HITROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_hitrou_create(par);
  if (!gen) { 
    free(par); 
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

  /* compute bounding rectangle */
  if (_unur_hitrou_rectangle(gen)!=UNUR_SUCCESS) {
    free(par); _unur_hitrou_free(gen);
    return NULL;
  }

  /* set initial point inside RoU shape */
  GEN.point[GEN.dim]=GEN.vmax/2.; /* all other coordinates are already 0 */


#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_hitrou_debug_init(gen);
#endif

  /* free parameters */
  free(par);

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
  rr->dim    = GEN.dim;
  rr->umin   = GEN.umin;
  rr->umax   = GEN.umax;
  rr->r      = GEN.r;
  rr->center = GEN.center; 
  rr->genid  = gen->genid;
  
  /* calculate bounding rectangle */
  _unur_mrou_rectangle_compute(rr);


  if (!(gen->set & HITROU_SET_V)) {
     /* user has not provided any upper bound for v */
     GEN.vmax = rr->vmax;
  }

  if (!(gen->set & HITROU_SET_U)) {
    /* user has not provided any bounds for u */
    for (d=0; d<GEN.dim; d++) {
      GEN.umin[d] = rr->umin[d];
      GEN.umax[d] = rr->umax[d];
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
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_HITROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & HITROU_VARFLAG_VERIFY) ? _unur_hitrou_sample_check : _unur_hitrou_sample_cvec;

  gen->destroy = _unur_hitrou_free;
  gen->clone = _unur_hitrou_clone;

  /* allocate memory for u-boundary arrays */
  GEN.umin = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */
  GEN.umax = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */

  /* allocate memory for random direction vector in (U,V) space */
  GEN.direction = _unur_xmalloc( (PAR.dim+1) * sizeof(double)); 
  
  /* allocate memory for current interior point */
  GEN.point = _unur_xmalloc( (PAR.dim+1) * sizeof(double)); 

  /* copy parameters into generator object */
  GEN.dim   = PAR.dim;              /* dimension */
  GEN.r     = PAR.r;                /* r-parameter of the hitrou method */  
  GEN.vmax  = PAR.vmax;             /* upper v-boundary of bounding rectangle */
  GEN.skip  = PAR.skip;             /* number of skipped poins in chain */
  
  if (PAR.umin != NULL) memcpy(GEN.umin, PAR.umin, GEN.dim * sizeof(double));
  if (PAR.umax != NULL) memcpy(GEN.umax, PAR.umax, GEN.dim * sizeof(double));

  /* get center of the distribution */
  GEN.center = unur_distr_cvec_get_center(gen->distr);

  /* initialize parameters */
  for (d=0; d<GEN.dim+1; d++) GEN.point[d]=0.;

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
#define CLONE clone->data.hitrou

  struct unur_gen *clone;


  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HITROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory for arrays */
  CLONE.umin = _unur_xmalloc( GEN.dim * sizeof(double));
  CLONE.umax = _unur_xmalloc( GEN.dim * sizeof(double));
 
  CLONE.point = _unur_xmalloc( (GEN.dim+1) * sizeof(double));
  CLONE.direction = _unur_xmalloc( (GEN.dim+1) * sizeof(double));

  /* copy parameters into clone object */
  CLONE.skip = GEN.skip;
  CLONE.r = GEN.r;
  CLONE.vmax = GEN.vmax;
  memcpy(CLONE.umin, GEN.umin, GEN.dim * sizeof(double));
  memcpy(CLONE.umax, GEN.umax, GEN.dim * sizeof(double));
  memcpy(CLONE.point, GEN.point, (GEN.dim+1) * sizeof(double));
  memcpy(CLONE.direction, GEN.direction, (GEN.dim+1) * sizeof(double));

  /* copy data */
  CLONE.center = unur_distr_cvec_get_center(clone->distr);

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
  double U,V;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double lambda, lmin, lmax; /* lambda parameters of line */
  long skip;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  
  COOKIE_CHECK(gen,CK_HITROU_GEN,RETURN_VOID); 

  dim = GEN.dim;

  for (skip=0; skip<=GEN.skip; skip++) {
 
    /* generate random direction vector in (U,V) space */
    _unur_hitrou_random_direction(gen, dim+1, GEN.direction);

    /* some random initialization in order to avoid compiler warnings */
    lmin=0; lmax=0;

    /* calculate lambda parameters for the intersections with v=vmax and v=0 */
    lambda = (GEN.vmax - GEN.point[dim]) / GEN.direction[dim];
    if (lambda>0) lmax = lambda;    
    if (lambda<0) lmin = lambda;    
    
    lambda = (0 - GEN.point[dim]) / GEN.direction[dim];
    if (lambda>0) lmax = lambda;    
    if (lambda<0) lmin = lambda;    

    /* calculate the intersections alont all other coordinate directions */
    for (d=0; d<dim; d++) {
      lambda = (GEN.umin[d] - GEN.point[d]) / GEN.direction[d];
      if (lambda>0 && lambda<lmax) lmax = lambda;    
      if (lambda<0 && lambda>lmin) lmin = lambda;    
     
      lambda = (GEN.umax[d] - GEN.point[d]) / GEN.direction[d];
      if (lambda>0 && lambda<lmax) lmax = lambda;    
      if (lambda<0 && lambda>lmin) lmin = lambda;    
    }

    while (1) {
  
      /* TODO: ntrials++ */
      lambda = lmin + (lmax-lmin) * _unur_call_urng(gen->urng);
      if (lambda>0) lmax=lambda;
      if (lambda<0) lmin=lambda;

      V = GEN.point[dim] + lambda * GEN.direction[dim];
      for (d=0; d<dim; d++) {
        U = GEN.point[d] + lambda * GEN.direction[d];
        vec[d] = U/pow(V,GEN.r) + GEN.center[d];
      }
    
      /* UV[] inside domain ? */
      if (V <= pow(PDF(vec),1./(GEN.r * dim + 1.))) {
        /* update current point */
        GEN.point[dim] = V;
        for (d=0; d<dim; d++) GEN.point[d] += lambda * GEN.direction[d] ; 
      
        break; /* jump out of the while() loop */
      }
    }

  }

  return;
} /* end of _unur_hitrou_sample() */

/*---------------------------------------------------------------------------*/

void
_unur_hitrou_sample_check( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random sample vector (return)                              */
     /*----------------------------------------------------------------------*/
{ 
  _unur_hitrou_sample_cvec(gen, vec);

} /* end of _unur_hitrou_sample_check() */

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

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN.umin) free(GEN.umin); 
  if (GEN.umax) free(GEN.umax);
  if (GEN.direction) free(GEN.direction);
  if (GEN.point) free(GEN.point);
  _unur_generic_free(gen);

} /* end of _unur_hitrou_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

void 
_unur_hitrou_random_direction( struct unur_gen *gen, 
                               int dim, double *direction)
     /*----------------------------------------------------------------------*/
     /* generte a random direction vector (not necessarily unit vector)      */
     /*----------------------------------------------------------------------*/
{
  int d;
  
  for (d=0; d<dim; d++) {
    direction[d] = unur_sample_cont(NORMAL);
  }

}

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
  dim = GEN.dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = hitrou (hit and run ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_hitrou_sample",gen->genid);
  if (gen->variant & HITROU_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(log,"%s: r-parameter = %g",gen->genid, GEN.r);
  _unur_print_if_default(gen,HITROU_SET_R);
  fprintf(log,"\n%s:\n",gen->genid);

  /* print center */
  _unur_matrix_print_vector( GEN.dim, GEN.center, "center =", log, gen->genid, "\t   ");

  /* print bounding rectangle */
  fprintf(log,"%s: Rectangle:",gen->genid);
  if (!((gen->set & HITROU_SET_U) && (gen->set & HITROU_SET_V)))
    fprintf(log,"\t[computed]");
  else 
    fprintf(log,"\t[input]");
  fprintf(log,"\n");

  vol = GEN.vmax;
  fprintf(log,"%s:\tvmax = %g\n",gen->genid, GEN.vmax);
  for (d=0; d<dim; d++) {
    vol *= (GEN.umax[d]-GEN.umin[d]);
    fprintf(log,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
	    d, d, GEN.umin[d], GEN.umax[d]);
  }
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN.r*GEN.dim+1));
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_hitrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
