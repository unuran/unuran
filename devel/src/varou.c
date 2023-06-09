/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      varou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    adaptive ratio-of-uniforms method                            *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ...                                                                  *
 *      Produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      PDF, dPDF                                                            *
 *   OPTIONAL:                                                               *
 *      ...                                                                  *
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
 *   [1] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h> 
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "varou.h"
#include "varou_struct.h"

/*---------------------------------------------------------------------------*/
/* Debugging flag for additional data dump                                   */

#define VAROU_DEBUG 0 

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define VAROU_VARFLAG_VERIFY   0x002u   /* run verify mode                   */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define VAROU_SET_SOMETHING    0x001u     /* set something                   */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VAROU"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

/* this parameter define the maximal number of cones that will be created.   */
#define VAROU_MAX_CONES_DEFAULT 1000

/*---------------------------------------------------------------------------*/

/* number of intervals to find an initial bracketing */
int N_ALPHA_INTERVALS=20; /* should we move this parameter into generator object ? */

/* accuracy parameter for the brent algorithm              */
/* this is used in the calculation of minimal cone volumes */
#define VAROU_BRENT_EPSILON 0.001

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_varou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_varou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_varou_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_varou_sample_check( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_varou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_varou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_varou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif


/*---------------------------------------------------------------------------*/
/*   Auxilliary Routines                                                     */
/*---------------------------------------------------------------------------*/

static struct unur_varou_cone *_unur_varou_cone_new(int dim);
/*---------------------------------------------------------------------------*/
/* allocate memory for cone structure and data                               */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cones_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute initial (2*dim+1) verteces on (half) unit sphere and 2^dim cones  */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cone_copy(int dim, struct unur_varou_cone *c_destination, 
                                           struct unur_varou_cone *c_source);
/*---------------------------------------------------------------------------*/
/* copies the structure and data block between two cones                     */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cones_split( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* split cones with volume=inf or volume >= mean_volume                      */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cones_prepare( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* calculation of cummulated volumes, eventually sorting, guiding table ...   */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cone_parameters( struct unur_gen *gen, 
            struct unur_varou_cone *c, double alpha);
/*---------------------------------------------------------------------------*/
/* adjust cone parameters of cone (norms, volume, etc)                       */
/*---------------------------------------------------------------------------*/

static void _unur_varou_cone_set(struct unur_gen *gen, struct unur_varou_cone *c, 
            double *spoint, double *normal); 
/*---------------------------------------------------------------------------*/
/* sets the spoint and normal vectors in cone c                              */
/*---------------------------------------------------------------------------*/

static void _unur_varou_minimize_volume(struct unur_gen *gen, struct unur_varou_cone *c);
/*---------------------------------------------------------------------------*/
/* find appropriate tangent plane by choosing smallest cone volume           */
/*---------------------------------------------------------------------------*/

static double _unur_varou_volume(double alpha, void *p );
/*---------------------------------------------------------------------------*/
/* auxiliary function used in the computation of cone volumes                */
/*---------------------------------------------------------------------------*/

static void _unur_varou_sample_cone( struct unur_gen *gen, 
                                     struct unur_varou_cone *c, double *UV );
/*---------------------------------------------------------------------------*/
/* generate random vector UV uniformly distributed in the cone c             */
/*---------------------------------------------------------------------------*/

static int _unur_varou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static double _unur_varou_f( struct unur_gen *gen, double *uv); 
/*---------------------------------------------------------------------------*/
/* return value of PDF(u_0/v+c_0, u_1/v+c_1, ..., u_{dim-1}/v+c_{dim-1})     */
/* where (c_0, c_1, ...) is the center, u_i=uv[i] for i<dim and v=uv[dim]    */
/*---------------------------------------------------------------------------*/

static void _unur_varou_dF( struct unur_gen *gen, double *uv, double *dF);
/*---------------------------------------------------------------------------*/
/* calculate the value the gradient of v^(dim+1)-PDF(...) at the point (u,v) */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_varou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_varou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */     
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_varou_new( const struct unur_distr *distr )
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

  /* check PDF */
  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }
  
  /* check dPDF */
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dPDF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_varou_par) );
  COOKIE_SET(par,CK_VAROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR->dim = distr->dim;

  /* setting default number of cones to be used in the setup */
  PAR->max_cones = VAROU_MAX_CONES_DEFAULT; 

  /* set other default values */
  par->method   = UNUR_METH_VAROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_varou_init;

  return par;

} /* end of unur_varou_new() */

/*---------------------------------------------------------------------------*/

int
unur_varou_set_cones( struct unur_par *par, long ncones  )
     /*----------------------------------------------------------------------*/
     /* sets the number of cones to be used in the setup                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   ncones ... maximal number of cones to use                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   UNUR_FAILURE ... when number of cones is too small to provide      */
     /*                    an initialization (ncones < 2^dim)                */
     /*                                                                      */
     /* comment:                                                             */
     /*   when not set, the default VAROU_MAX_CONES_DEFAULT will be used     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VAROU );

  /* check if the number of cones sufficies for an initial splitting */
  if (pow(2., (double)PAR->dim) > ncones) return UNUR_FAILURE;
  
  /* setting the number of cones */
  PAR->max_cones = ncones ;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_varou_set_cones() */

/*---------------------------------------------------------------------------*/


int
unur_varou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, VAROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | VAROU_VARFLAG_VERIFY) : 
                            (par->variant & (~VAROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_varou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_varou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, VAROU, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cvec_error) 
    return UNUR_FAILURE;

  if (verify) {
    /* turn verify on */
    gen->variant |= VAROU_VARFLAG_VERIFY;
    SAMPLE = _unur_varou_sample_check;
  }
  else {
    /* turn verify bounding rectangle off */
    gen->variant &= ~VAROU_VARFLAG_VERIFY;
    SAMPLE = _unur_varou_sample_cvec;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_varou_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_varou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_VAROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VAROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_varou_create(par);

  /* free parameters */
  _unur_par_free(par);

  if (!gen) { return NULL; }

  /* check if initial cones can be computed */ 
  if (pow(2., (double)GEN->dim) > GEN->max_cones) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"max_cones too small");
    _unur_varou_free(gen);
    return NULL;
  }
  
  /* compute bounding rectangle */
  if (_unur_varou_rectangle(gen)!=UNUR_SUCCESS) {
    _unur_varou_free(gen);
    return NULL;
  }

  /* compute initial triangulation ... etc  */
  _unur_varou_cones_init(gen);

  /* cone splitting ... */
  _unur_varou_cones_split(gen);

  /* calculation of cummulated volumes, eventually sorting ... */
  _unur_varou_cones_prepare(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_varou_debug_init(gen);
#endif

  /* do we have an infinite hut ? */
  if (_unur_isinf(GEN->cone_list[GEN->n_cone-1]->sum_volume) ) {
    _unur_varou_free(gen);
    return NULL;
  }

  return gen;

} /* end of _unur_varou_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_varou_create( struct unur_par *par )
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

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VAROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_varou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_VAROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & VAROU_VARFLAG_VERIFY) ? _unur_varou_sample_check : 
                                                   _unur_varou_sample_cvec;

  gen->destroy = _unur_varou_free;
  gen->clone = _unur_varou_clone;

  /* allocate memory for u-boundary arrays */
  GEN->umin = _unur_xmalloc( PAR->dim * sizeof(double)); /* bounding rectangle */
  GEN->umax = _unur_xmalloc( PAR->dim * sizeof(double)); /* bounding rectangle */

  /* copy parameters into generator object */
  GEN->dim   = PAR->dim;  /* dimension */
  
  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);

  /* maximal number of verteces and cones                                      */
  /* during the triangulation of the upper half-unit-sphere, the parameter     */
  /* max_verteces define the maximal number of verteces that will be created   */
  /* during the adaptive process.                                              */
  /* initial triangulation create v0=2*dim+1 verteces and c0=2^dim cones.      */
  /* each subsequently added vertex lead to a cone-splitting and thus to one   */
  /* additional cone structure.                                                */
  /* if we thus have v verteces (v>=v0) then the number of cones c is given by */
  /*   c = c0 + (v-v0)                                                         */
  /*     = 2^dim + (v - (2*dim+1))                                             */
  /* i.e.                                                                      */
  /*   v_max =  c_max - 2^dim + (2*dim+1)                                      */
  /*         <= c_max + 1                                                      */
  GEN->max_cones = PAR->max_cones; 
  GEN->max_verteces = GEN->max_cones + 1;
  
  /* allocate memory for the list of the verteces on the half-unit-sphere */
  GEN->vertex_list = (double **) _unur_xmalloc(GEN->max_verteces*sizeof(long));
  GEN->n_vertex = 0;
  
  /* allocate memory for the list of the verteces on the half-unit-sphere */
  GEN->cone_list = (struct unur_varou_cone **) 
                  _unur_xmalloc(GEN->max_cones*sizeof(long));
  GEN->n_cone = 0;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_varou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_varou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_varou_gen*)clone->datap)

  struct unur_gen *clone;
  long i;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VAROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory and copy data for u-arrays */
  CLONE->umin = _unur_xmalloc( GEN->dim * sizeof(double));
  CLONE->umax = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy(CLONE->umin, GEN->umin, GEN->dim * sizeof(double));
  memcpy(CLONE->umax, GEN->umax, GEN->dim * sizeof(double));

  /* allocate memory and copy data for the vertex list */
  CLONE->vertex_list = (double **) _unur_xmalloc(GEN->max_verteces*sizeof(long));
  for (i=0; i<GEN->n_vertex; i++) {
    CLONE->vertex_list[i] = _unur_xmalloc((GEN->dim+1)*sizeof(double));
    memcpy(CLONE->vertex_list[i], GEN->vertex_list[i], (GEN->dim+1)*sizeof(double));
  }
		    
  /* allocate memory and copy data for the cone list */
  CLONE->cone_list = (struct unur_varou_cone **) 
                    _unur_xmalloc(GEN->max_cones*sizeof(long));
  for (i=0; i<GEN->n_cone; i++) {
    CLONE->cone_list[i] = _unur_varou_cone_new(GEN->dim);
    _unur_varou_cone_copy(GEN->dim, CLONE->cone_list[i], GEN->cone_list[i]);
  }

  /* copy other data */
  CLONE->dim = GEN->dim;
  CLONE->vmax = GEN->vmax;
  CLONE->center = unur_distr_cvec_get_center(clone->distr);
  CLONE->n_vertex = GEN->n_vertex;
  CLONE->max_cones = GEN->max_cones;
  CLONE->max_verteces = GEN->max_verteces;
 
  return clone;

#undef CLONE
} /* end of _unur_varou_clone() */

/*---------------------------------------------------------------------------*/

int
_unur_varou_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{ 
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  long ic, ic_sample=0; /* cone index */
  double *UV; /* (dim+1) uniformly distributed in cone */
  double vol;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VAROU_GEN,UNUR_ERR_COOKIE); 

  dim = GEN->dim;

  /* theese are the coordinates in the u-v space : u[i]=UV[i] and v=UV[dim] */
  UV = _unur_xmalloc((dim+1)*sizeof(double));
  
  while(1) {
    /* vol is a uniformly distributed variable < sum of all cone volumes */
    vol = GEN->cone_list[GEN->n_cone-1]->sum_volume * _unur_call_urng(gen->urng);
    
    /* obtain cone index of cone in which we should sample */
    for (ic=0; ic<GEN->n_cone; ic++) {  
      if ( GEN->cone_list[ic]->sum_volume >= vol ) { 
        ic_sample=ic; 
        break; 
      }
    }
  
    /* sampling in cone #ic_sample */
    _unur_varou_sample_cone(gen, GEN->cone_list[ic_sample], UV);

    /* obtaining the (u/v)-coordinates */ 
    for (d=0; d<dim; d++) {  
      vec[d] = UV[d]/UV[dim] + GEN->center[d] ;
    }
  
    /* check if sample point is inside the potato volume */
    if ( pow(UV[dim], 1.+dim) <= PDF(vec) ) {
       free(UV);
       return UNUR_SUCCESS;
    }  
  }
  
   
} /* end of _unur_varou_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_varou_sample_check( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random sample vector (return)                              */
     /*----------------------------------------------------------------------*/
{ 
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  long ic, ic_sample=0; /* cone index */
  double *UV; /* (dim+1) uniformly distributed in cone */
  double vol;
  double norm_factor;
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VAROU_GEN,UNUR_ERR_COOKIE); 

  dim = GEN->dim;
  
  /* theese are the coordinates in the u-v space : u[i]=UV[i] and v=UV[dim] */
  UV = _unur_xmalloc((dim+1)*sizeof(double));
  
  while(1) {
    /* vol is a uniformly distributed variable < sum of all cone volumes */
    vol = GEN->cone_list[GEN->n_cone-1]->sum_volume * _unur_call_urng(gen->urng);
  
    /* obtain cone index of cone in which we should sample */
    for (ic=0; ic<GEN->n_cone; ic++) {  
      if ( GEN->cone_list[ic]->sum_volume >= vol ) {
        ic_sample=ic;
	break;
      }
    }
        
    /* sampling in cone #ic_sample */
    _unur_varou_sample_cone(gen, GEN->cone_list[ic_sample], UV);

    /* obtaining the (u/v)-coordinates */ 
    for (d=0; d<dim; d++) {  
      vec[d] = UV[d]/UV[dim] + GEN->center[d] ;
    }

    /* check if sample point is inside the potato volume */
    if ( pow(UV[dim], 1.+dim) <= PDF(vec) ) {

       /* check hat */
       norm_factor = _unur_vector_scalar_product( dim+1, 
                            GEN->cone_list[ic_sample]->normal, 
                            GEN->cone_list[ic_sample]->spoint )  
                       / _unur_vector_scalar_product( dim+1, 
                            GEN->cone_list[ic_sample]->normal, 
                            UV  ) ; 
	   
       if (norm_factor<1) { 
         _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
       }
       
       free(UV);
       return UNUR_SUCCESS;
      
    }  
  }
  
} /* end of _unur_varou_sample_check() */

/*---------------------------------------------------------------------------*/

void
_unur_varou_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  long i;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_VAROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->umin) free(GEN->umin); 
  if (GEN->umax) free(GEN->umax);

  /* verteces */
  if (GEN->vertex_list) {
    for (i=0; i<GEN->n_vertex; i++) {
      if (GEN->vertex_list[i]) free(GEN->vertex_list[i]);
    }
    free(GEN->vertex_list);
  }

  /* cones */
  if (GEN->cone_list) {
    for (i=0; i<GEN->n_cone; i++) {
      if (GEN->cone_list[i]) free(GEN->cone_list[i]);
    }
    free(GEN->cone_list);
  }

  _unur_generic_free(gen);

} /* end of _unur_varou_free() */

/*---------------------------------------------------------------------------*/



/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

struct unur_varou_cone *
_unur_varou_cone_new(int dim)
     /*----------------------------------------------------------------------*/
     /* allocate memory for cone structure and data                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_varou_cone *s;   
  char *b;

  /* allocate memory block for struct + data */
  s = (struct unur_varou_cone *) 
      _unur_xmalloc( sizeof(struct unur_varou_cone)
                   + sizeof(long)*(dim+1)      /* vertex index  */
                   + sizeof(double)*(dim+1)    /* lengths       */
                   + sizeof(double)*(dim+1)    /* spoint        */
                   + sizeof(double)*(dim+1) ); /* normal vector */

  /* we need a char * to get the following pointer arithmetic right */
  b = (char *) s; 
  s->index   = (long *)   (b + sizeof(struct unur_varou_cone));
  s->length  = (double *) (b + sizeof(struct unur_varou_cone)
                             + sizeof(long)*(dim+1) );
  s->spoint  = (double *) (b + sizeof(struct unur_varou_cone)
                             + sizeof(long)*(dim+1) 
			     + sizeof(double)*(dim+1) );
  s->normal  = (double *) (b + sizeof(struct unur_varou_cone)
                             + sizeof(long)*(dim+1) 
			     + sizeof(double)*(dim+1) 
			     + sizeof(double)*(dim+1) );

  return s;
} /* end of _unur_varou_cone_new() */

/*---------------------------------------------------------------------------*/

void 
_unur_varou_cone_copy(int dim, struct unur_varou_cone *c_destination, 
                               struct unur_varou_cone *c_source)
{
     /*----------------------------------------------------------------------*/
     /* copies the structure and data block between two cones                */
     /*----------------------------------------------------------------------*/

  long block_size;
  char *b;

  block_size= sizeof(struct unur_varou_cone)
            + sizeof(long)*(dim+1)      /* vertex index  */
            + sizeof(double)*(dim+1)    /* lengths       */
	    + sizeof(double)*(dim+1)    /* spoint        */
	    + sizeof(double)*(dim+1) ;  /* normal vector */

  memcpy(c_destination, c_source, block_size);

  /* adjusting pointers in detination cone */
  b = (char *) c_destination; 
  c_destination->index   = (long *)   (b + sizeof(struct unur_varou_cone));
  c_destination->length  = (double *) (b + sizeof(struct unur_varou_cone)
                                         + sizeof(long)*(dim+1) );
  c_destination->spoint  = (double *) (b + sizeof(struct unur_varou_cone)
                                         + sizeof(long)*(dim+1) 
			                 + sizeof(double)*(dim+1) );
  c_destination->normal  = (double *) (b + sizeof(struct unur_varou_cone)
                                         + sizeof(long)*(dim+1) 
			                 + sizeof(double)*(dim+1) 
			                 + sizeof(double)*(dim+1) );


  return;
}

/*---------------------------------------------------------------------------*/

void
_unur_varou_cones_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute initial (2*dim+1) verteces on the upper half-unit-sphere.    */
     /* for each u-coordinate dimension, d, two verteces are initialized     */
     /* having u_d=-1 and u_d=+1 (i.e. along the -/+ d-coordinate direction) */
     /* for the v-coordinate we have a vertex at the north-pole.             */
     /* calculate the 2^dim initial cones spanned by these verteces.         */
     /*----------------------------------------------------------------------*/
{
  int i, ic, iv, dim;

  dim = GEN->dim;

  /* The first vertex is at the north-pole, i.e. v=1, u_i=0 for 0<=i<dim */
  GEN->vertex_list[0] = _unur_vector_new(dim+1);
  GEN->vertex_list[0][dim]=1.;
     
  /* 2*dim verteces along the positive and negative coordinate directions */
  for (i=0; i<dim; i++) {
    iv=1+2*i; /* current running vertex number */
    GEN->vertex_list[iv] = _unur_vector_new(dim+1);
    GEN->vertex_list[iv][i] = -1.;
    GEN->vertex_list[iv+1] = _unur_vector_new(dim+1);
    GEN->vertex_list[iv+1][i] = +1.;
  }

  /* updating the current number of allocated unit verteces */
  GEN->n_vertex  = 2*dim + 1;

  /* updating (initial) number of cones */
  GEN->n_cone = pow(2., (double)dim);

  /* setting parameters for the 2^dim initial cones */
  for (ic=0; ic<GEN->n_cone; ic++) {
    GEN->cone_list[ic] = _unur_varou_cone_new(dim);
    GEN->cone_list[ic]->index[0]=0; /* initial cones contain the north-pole */
    GEN->cone_list[ic]->length[0]=GEN->vmax; 

    GEN->cone_list[ic]->unit_volume=1.; /* 1/(dim+1)! is calculated below */
    for (i=0; i<dim; i++) {
      /* obtain vertex indices of the ic'th cone */
      GEN->cone_list[ic]->index[i+1]=((GEN->n_cone-ic-1)>>i & 1)?(1+2*i):(2+2*i); 
      GEN->cone_list[ic]->unit_volume *= 1./(2.+i);
    }

    /* setup tangent plane and calculate cone volume */
    _unur_varou_minimize_volume(gen, GEN->cone_list[ic]);
  }

  return;
} /* end of _unur_varou_cones_init() */
 
/*---------------------------------------------------------------------------*/

void
_unur_varou_minimize_volume(struct unur_gen *gen, struct unur_varou_cone *c)
     /*----------------------------------------------------------------------*/
     /* find appropriate tangent plane by choosing smallest cone volume      */
     /*----------------------------------------------------------------------*/
{
  long i;
  double alpha, alpha_min;
  double alpha_start, alpha_end; /* interval in which a minimum is sought */
  int finite_volume=0; /* set to 1 as soon as a finite volume is found */
  struct unur_funct_generic fs;
  
  /* filling function structure to be used by the brent algorithm */
  fs.f = (UNUR_FUNCT_GENERIC *) _unur_varou_volume;
  fs.params = gen;

  /* storing the actual cone */
  GEN->current_cone = c;

  /* setting initial values */
  alpha_start=0.;
  alpha_end=0.;
  alpha=0.;

  /* initial bracketing of alpha between alpha_start and alpha_end             */
  /* alpha_start is the first value of alpha, where the volume becomes  finite */
  /* alpha_end   is tha last  value of alpha, where the volume is still finite */
  for (i=0; i<N_ALPHA_INTERVALS; i++) {
    alpha = i/((double)N_ALPHA_INTERVALS); /* alpha in [0,1) */
    _unur_varou_cone_parameters(gen, c, alpha); 
    if ( !_unur_isinf(c->volume) ) {
      if ( finite_volume==0 ) { 
        finite_volume=1; 
        alpha_start=alpha;
      }
      alpha_end=alpha;
    }
  }
  
  /* minimization step */
  if (alpha_start == alpha_end) { 
    alpha_min=alpha_start;
  }
  else {
    alpha=_unur_util_brent(fs, alpha_start, alpha_end, (alpha_start+alpha_end)/2., 
                           VAROU_BRENT_EPSILON);

    /* just to be sure, that we have obtained a good value for alpha */
    alpha_min=alpha_start;
    if (alpha_start <= alpha && alpha <= alpha_end)  alpha_min=alpha;
  }
  
  _unur_varou_cone_parameters(gen, c, alpha_min); 

  return;
} /* end of _unur_varou_minimize_volume() */

/*---------------------------------------------------------------------------*/

/* double  */
/* _unur_varou_calculate_unit_volume(struct unur_gen *gen, struct unur_varou_cone *c ) */
/*      /\* calculation of the unit volume using determinant calculation *\/ */
/*      /\* currently this is only implemented for dim=2                 *\/ */
/*      /\* and is only needed for testing purposes (may be removed)     *\/ */
/* { */
/*   int dim; */
/*   double x11,x12,x13, x21,x22,x23, x31,x32,x33; */
/*   double det; */
  
/*   dim = GEN->dim; */
  
/*   x11=GEN->vertex_list[ c->index[0] ][0]; */
/*   x12=GEN->vertex_list[ c->index[0] ][1]; */
/*   x13=GEN->vertex_list[ c->index[0] ][2]; */

/*   x21=GEN->vertex_list[ c->index[1] ][0]; */
/*   x22=GEN->vertex_list[ c->index[1] ][1]; */
/*   x23=GEN->vertex_list[ c->index[1] ][2]; */

/*   x31=GEN->vertex_list[ c->index[2] ][0]; */
/*   x32=GEN->vertex_list[ c->index[2] ][1]; */
/*   x33=GEN->vertex_list[ c->index[2] ][2]; */

/*   det = fabs( x11*x22*x33 + x12*x23*x31 + x13*x21*x32 */
/*         - x31*x22*x13 - x32*x23*x11 - x33*x21*x12 ) / 6.; */

/*   return det; */
  
/* } */

/*---------------------------------------------------------------------------*/

void
_unur_varou_cones_split( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* split cones with volume=inf or volume >= mean_volume                 */
     /*----------------------------------------------------------------------*/
{
  int dim;
  long ic; /* running cone index */
  long nc; /* last cone index */
  long nv; /* last vertex index */
  long iv1, iv2, it; /* vertex indices */
  double v, vt; /* v-coordinates of the verteces */
  double sum_volume;  /* volume sum of all bounded cones */
  double potato_volume;  /* volume of potato = 1/(dim+1)  */
  long   n_bounded;   /* number of bounded cones */
  double mean_volume; /* sum_volume / n_bounded */
  double norm; 
  long i;
  
  struct unur_varou_cone *c1; 
  struct unur_varou_cone *c2;
  
  dim = GEN->dim;

  /* volume enclosed by surface v^(dim+1)=pdf(u,v) */
  potato_volume=1./(dim+1.);

#if VAROU_DEBUG
  /* quick and dirty debugging */
  printf("n_cone; n_inf; vol; vol-vol0; rho\n"); 
#endif

  while (GEN->n_cone <= GEN->max_cones) {

    /* obtaining the number of bounded cones and their volume sum */
    n_bounded = 0;
    sum_volume = 0.;
    for (ic=0; ic<GEN->n_cone; ic++) {
      if (!_unur_isinf(GEN->cone_list[ic]->volume)) {
        n_bounded++;
        sum_volume += GEN->cone_list[ic]->volume;
      }
    }

    /* calculation of mean_volume */
    if (n_bounded==GEN->n_cone) {
      mean_volume = sum_volume / n_bounded ;
    }
    else {
      /* some volumes are (still) infinite */
      /* we set the mean volume to be higher than every bounded volume */
      mean_volume=sum_volume+1.; 
    }

#if VAROU_DEBUG
    /* quick and dirty debugging */
    printf("%8ld; %8ld; %e; %e; %e\n", 
        GEN->n_cone, GEN->n_cone - n_bounded, sum_volume, 
	sum_volume-potato_volume, sum_volume/potato_volume - 1.);
#endif

    /* splitting of cones with volume >= mean_volume */
    nc=GEN->n_cone; 
    for (ic=0; ic<nc; ic++) {
      if (GEN->cone_list[ic]->volume >= mean_volume ) {
	
        /* obtaining top vertex  */
        vt=0.; 
        it=0;
        for (i=0; i<=dim; i++) {
          v = GEN->vertex_list[ GEN->cone_list[ic]->index[i] ][dim];
          if (v>=vt) {vt=v; it=i;}
        }
    
        /* determining two distinct verteces of edge to be splitted    */
	/* we exclude the top vertex as long as the volume is infinite */
	iv1=0; iv2=0;
        while (iv1==iv2 || 
	      (_unur_isinf(GEN->cone_list[ic]->volume) && (iv1 == it || iv2 == it)) ) {
          iv1=(long) (dim+1)*_unur_call_urng(gen->urng);
          iv2=(long) (dim+1)*_unur_call_urng(gen->urng);
        }

        /*------------------------------------------------------------*/
        /* create new vertex */
        nv=GEN->n_vertex;
        if (nv>=GEN->max_verteces) goto done_splitting; 
	GEN->vertex_list[nv] = _unur_vector_new(dim+1);

        /* vertex coordinates are mean values of the iv1 and iv2 coordinates */
        for (i=0; i<=dim; i++) {
          GEN->vertex_list[nv][i] = 
	    (GEN->vertex_list[GEN->cone_list[ic]->index[iv1]][i]
	    +GEN->vertex_list[GEN->cone_list[ic]->index[iv2]][i])/2.;
        }

        /* scale this new vertex-vector to have norm=1 */
        norm = _unur_vector_norm(dim+1, GEN->vertex_list[nv]);
        for (i=0; i<=dim; i++)  GEN->vertex_list[nv][i] /= norm; 
      
        GEN->n_vertex++; /* adjust current number of verteces */
        
	/*------------------------------------------------------------*/
        
	if (GEN->n_cone>=GEN->max_cones) goto done_splitting; 
        
	/* create two new cones c1 and c2 */
        c1=_unur_varou_cone_new(dim);
        c2=_unur_varou_cone_new(dim);

        /* initially they are set to be exact copies of the original cone to 
	   be splitted */
        _unur_varou_cone_copy(dim, c1, GEN->cone_list[ic]);        
        _unur_varou_cone_copy(dim, c2, GEN->cone_list[ic]);        

        /* adjust first cone */
        c1->index[iv1] = nv;      
        c1->unit_volume /= 2.*norm; 

        /* adjust second cone */
        c2->index[iv2] = nv;      
        c2->unit_volume = c1->unit_volume; 
     
        /* searching for volume minimum of c1 and c2 */
        _unur_varou_minimize_volume(gen, c1);
        _unur_varou_minimize_volume(gen, c2);

	/* comparing the old volume with the sum of the two new volumes  */
        if (!_unur_isinf(GEN->cone_list[ic]->volume) &&
	    GEN->cone_list[ic]->volume < ( c1->volume + c2->volume ) ) {
           /* we have split a bounded cone -> obtaining something bigger  */
	  
	   /* re-use the hat of the original cone */ 
           _unur_varou_cone_set(gen, c1, GEN->cone_list[ic]->spoint,
                                         GEN->cone_list[ic]->normal);  
           _unur_varou_cone_set(gen, c2, GEN->cone_list[ic]->spoint,
                                         GEN->cone_list[ic]->normal);  
	}

        /* c1 ==> original cone, c2 ==> new cone */
        _unur_varou_cone_copy(dim, GEN->cone_list[ic], c1);        
        GEN->cone_list[GEN->n_cone]= c2;        

        free(c1); /* not c2 !!! */
        GEN->n_cone++; /* adjust current number of cones */

      }
    } 

  }

done_splitting:

  /* obtaining the number of bounded cones and their volume sum */
  /* TODO: move calculation of volume_sum into own subroutine ... */
  n_bounded = 0;
  sum_volume = 0.;
  for (ic=0; ic<GEN->n_cone; ic++) {
    if (!_unur_isinf(GEN->cone_list[ic]->volume)) {
      n_bounded++;
      sum_volume += GEN->cone_list[ic]->volume;
    }
  }

#if VAROU_DEBUG
  printf("%8ld; %8ld; %e; %e; %e\n", 
         GEN->n_cone, GEN->n_cone - n_bounded, sum_volume, 
         sum_volume-potato_volume, sum_volume/potato_volume - 1.);
#endif

  return;
} /* end of _unur_varou_cones_split() */

/*---------------------------------------------------------------------------*/

void
_unur_varou_cone_parameters( struct unur_gen *gen, struct unur_varou_cone *c, double alpha)
     /*----------------------------------------------------------------------*/
     /* adjust cone parameters of cone (norms, volume, etc)                  */
     /*----------------------------------------------------------------------*/
{
  int dim;
  double v,vt;
  double *p; /* position vector to surface */
  double *t; /* top vertex vector */
  double *ft; /* centre vector of face opposite to top vertex */
  double normp, normn;
  long i, it, itop, iv;

  dim = GEN->dim;
  
  /* obtaining maximal v coordinate of all cone verteces */
  vt=0.; 
  it=0;
  for (i=0; i<=dim; i++) {
    v = GEN->vertex_list[ c->index[i] ][dim];
    if (v>=vt) {vt=v; it=i;}
  }
    
    
  /*------------------------------------------------------------*/
      
  /* top vertex vector */
  t=_unur_vector_new(dim+1); 
  itop =  c->index[it];
  for (i=0; i<=dim; i++) {
    t[i]=GEN->vertex_list[itop][i];
  }

  /*------------------------------------------------------------*/
  
  /* calculating centre of face opposite of the top vector */
  ft=_unur_vector_new(dim+1); 
  for (iv=0; iv<=dim; iv++) {
    if (iv != it) {
      for (i=0; i<=dim; i++) {
        ft[i] += GEN->vertex_list[ c->index[iv] ][i] / dim; 
      }
    }
  }
   
  /*------------------------------------------------------------*/

  p=_unur_vector_new(dim+1); /* position vector to surface */

  /* setting p[] to be on the line connecting t[] and f[] */
  for (i=0; i<=dim; i++) {
    p[i] = t[i] + alpha * (ft[i]-t[i]);
  }

  /* scaling p[] to be unit vector */
  normp=_unur_vector_norm(dim+1, p);
  for (i=0; i<=dim; i++) {
     p[i] /= normp;
  }

  v = pow(_unur_varou_f(gen, p), 1./(1.+dim)) ;

  /* scaling p[] to touch surface */
  normp = v/p[dim]; 
    
  for (i=0; i<=dim; i++) {
    p[i] *= normp;
  }

  /*------------------------------------------------------------*/

  /* setting normal vector to tangent plane through p[]  */
  _unur_varou_dF(gen, p, c->normal);
  
  normn=_unur_vector_norm(dim+1, c->normal);    
  for (i=0; i<=dim; i++) c->normal[i] /= normn;

  /*------------------------------------------------------------*/
  
  _unur_varou_cone_set(gen, c, p, c->normal);

   free(t); free(p); free(ft); 

   return;
} /* end of _unur_varou_cone_parameters */

/*---------------------------------------------------------------------------*/

void
_unur_varou_cone_set(struct unur_gen *gen, struct unur_varou_cone *c, 
                     double *spoint, double *normal) 
     /*----------------------------------------------*/
     /* sets the spoint and normal vectors in cone c */
     /*----------------------------------------------*/
{
  int i, dim;
  double norm;

  dim = GEN->dim;

  /* setting spoint and normal vectors */
  for (i=0; i<=dim; i++) {
    c->spoint[i]=spoint[i];
    c->normal[i]=normal[i];
  }

  /* calculate lengths of all rays */
  for (i=0; i<=dim; i++) {
    norm = _unur_vector_scalar_product( dim+1, 
        c->normal, 
        c->spoint )  
         / _unur_vector_scalar_product( dim+1, 
        c->normal, 
        GEN->vertex_list[ c->index[i] ] ) ; 
      
    c->length[i] = (norm>0) ? norm: UNUR_INFINITY;
  }

  /* calculate total volume */
  c->volume = c->unit_volume;
  for (i=0; i<=dim; i++) {
    c->volume *= c->length[i];
  }

  return;

}

/*---------------------------------------------------------------------------*/

void 
_unur_varou_cones_prepare( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* calculation of cummulated volumes, eventually sorting, guiding table */
     /*----------------------------------------------------------------------*/
{
  long ic; /* running cone index */
  double sum_volume;

  /* calculate the cummulated volume sum */
  sum_volume = 0.;
  for (ic=0; ic<GEN->n_cone; ic++) {
    sum_volume += GEN->cone_list[ic]->volume;
    GEN->cone_list[ic]->sum_volume = sum_volume;    
  }

  return;
}


/*---------------------------------------------------------------------------*/

void
_unur_varou_sample_cone( struct unur_gen *gen, struct unur_varou_cone *c, double *UV )
     /*----------------------------------------------------------------------*/
     /* generate random vector UV uniformly distributed in the cone c        */
     /* having (dim+2) verteces (vert_0, ... , vert_dim, (0) )               */
     /*----------------------------------------------------------------------*/
{
  double  U; /* uniform variate */
  double *E; /* exponential variates */
  double *S; /* uniform spacings */
  double  E_sum; /* sum of all e[] */
  int i,j; /* index variables used in for loops */
  int dim;
  
  dim = GEN->dim; 
  
  S=_unur_xmalloc((dim+2)*sizeof(double));
  E=_unur_xmalloc((dim+2)*sizeof(double));
  
  /* calculating exponential variates */
  E_sum = 0; 
  for (i=0; i<dim+2; i++) {
    /* sample from U(0,1) */
    while ( (U = _unur_call_urng(gen->urng)) == 0.);
    /* sample from exponential distribution */
    E[i] = -log(1.-U);
    E_sum += E[i];	 
  }
  
  /* calculating uniform spacing */
  for (i=0; i<dim+2; i++) {
    S[i]=E[i]/E_sum;
  }

  /* calculating uniform vector in simplex */
  for (j=0; j<dim+1; j++) {
    UV[j] = 0;
    for (i=0; i<dim+1; i++) {
      /* the sum actually run over all (dim+2) verteces              */
      /* but, since all cones have the origin (0,0,...0) in common,  */
      /* a multiplication with its coordinates is not necessary      */
      /* and our sum run over the remaining (dim+1) verteces         */
      UV[j] += S[i] * GEN->vertex_list[c->index[i]][j] * c->length[i];
    }
  }

  free(S); free(E);

  return;
} /* end of _unur_varou_sample_cone() */

/*---------------------------------------------------------------------------*/

double
_unur_varou_volume(double alpha, void *p )
     /*----------------------------------------------------------------------*/
     /* auxiliary function used in the computation of cone volumes           */
     /*----------------------------------------------------------------------*/
{

  struct unur_gen *gen;
  gen = p; /* typecast from void* to unur_gen* */

  _unur_varou_cone_parameters(gen, GEN->current_cone, alpha);
 
  return -(GEN->current_cone->volume);

} /* end of _unur_varou_volume() */

/*---------------------------------------------------------------------------*/


double 
_unur_varou_f( struct unur_gen *gen, double *uv) 
     /*----------------------------------------------------------------------*/
     /* return calculated value of                                           */
     /* PDF(u_0/v+c_0, u_1/v+c_1, ..., u_{dim-1}/v+c_{dim-1})                */
     /* where (c_0, c_1, ...) is the center,                                 */
     /* u_i=uv[i] for i<dim and v=uv[dim].                                   */
     /*                                                                      */
     /* UNUR_INFINITY is returned when |v|<UNUR_EPSILON                      */
     /*----------------------------------------------------------------------*/
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double *x;  /* u/v+c vector   */
  double v;
  double f;   /* function value */

  dim = GEN->dim;
  v=uv[dim];
  if (fabs(v) <= UNUR_EPSILON) return UNUR_INFINITY;

  x = _unur_vector_new(dim);
  for (d=0; d<dim; d++) {
    x[d] = uv[d] / v + GEN->center[d];
  }

  f = PDF(x);

  free(x);
  
  return f;
}

/*---------------------------------------------------------------------------*/

void 
_unur_varou_dF( struct unur_gen *gen, double *uv, double *dF )
     /*----------------------------------------------------------------------*/
     /* calculate the gradient of v^(dim+1)-PDF(...) at the point (u,v)      */
     /* the gradient is written into the (dim+1) double array dF[]           */
     /*                                                                      */
     /* NULL is returned when |v|<UNUR_EPSILON                               */
     /*----------------------------------------------------------------------*/
{
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double v;
  double *x; /* u/v+c vector     */
  double *dPDF;

  dim = GEN->dim;
  v=uv[dim]; 
  if (fabs(v) <= UNUR_EPSILON) {dF= NULL; return;}
   
  /* reserving memory for arrays */ 
  x = _unur_vector_new(dim);
  dPDF = _unur_vector_new(dim);
  
  /* calculating the u/v ratios */
  for (d=0; d<dim; d++) {
    x[d] = uv[d] / v + GEN->center[d];
  }

  /* dPDF at the point x[] */
  _unur_cvec_dPDF(dPDF, x, gen->distr);

  /* gradient coordinates of F() */
  dF[dim] = (1.+dim) * pow(v, (double)dim);
  for (d=0; d<dim; d++) {
    dF[d] = - dPDF[d] / v;
    dF[dim] += uv[d]*dPDF[d]/(v*v);
  }

  /* freeing allocated memory */
  _unur_vector_free(x);
  _unur_vector_free(dPDF);
  
  return;
}

/*---------------------------------------------------------------------------*/

int
_unur_varou_rectangle( struct unur_gen *gen )
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
  COOKIE_CHECK( gen,CK_VAROU_GEN, UNUR_ERR_COOKIE );

  /* Allocating and filling rou_rectangle struct */
  rr = _unur_mrou_rectangle_new();

  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->umin;
  rr->umax   = GEN->umax;
  rr->r      = 1.;
  rr->center = GEN->center;
  rr->genid  = gen->genid;

  /* calculate bounding rectangle */
  _unur_mrou_rectangle_compute(rr);
		      
  GEN->vmax = rr->vmax;
  for (d=0; d<GEN->dim; d++) {
    GEN->umin[d] = rr->umin[d];
    GEN->umax[d] = rr->umax[d];
  }

  free(rr);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_varou_rectangle() */

/*---------------------------------------------------------------------------*/





/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_varou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  long ic, iv; /* indices for cones and verteces */
  double vol;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VAROU_GEN,RETURN_VOID);
  log = unur_get_stream();
  dim = GEN->dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = varou (adaptive ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  
  if (dim>1) _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_varou_sample",gen->genid);
  if (gen->variant & VAROU_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */

  /* print center */
  _unur_matrix_print_vector( GEN->dim, GEN->center, "center =", log, gen->genid, "\t   ");

  /* print bounding rectangle */
  fprintf(log,"%s: Rectangle:",gen->genid);
  fprintf(log,"\t[computed]\n");
 
  vol = GEN->vmax;
  fprintf(log,"%s:\tvmax = %g\n",gen->genid, GEN->vmax);
  for (d=0; d<dim; d++) {
    vol *= (GEN->umax[d]-GEN->umin[d]);
    fprintf(log,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
                d, d, GEN->umin[d], GEN->umax[d]);
  }
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN->dim+1));
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s:\tsum_of_cone_volumes = %g\n", gen->genid,GEN->cone_list[GEN->n_cone-1]->sum_volume);
  fprintf(log,"%s:\n",gen->genid);

  /* vertex list on upper unit-half-sphere */
  fprintf(log,"%s: Number of unit verteces = %ld\n",gen->genid, GEN->n_vertex);
  fprintf(log,"%s: Vertex list:\n",gen->genid);
     
  for (iv=0; iv<GEN->n_vertex; iv++) {
    fprintf(log,"%s:\t%6ld : (", gen->genid, iv);
    for (d=0; d<dim; d++) {
       fprintf(log,"% f,", GEN->vertex_list[iv][d]);
    }   
    fprintf(log,"% f)\n", GEN->vertex_list[iv][dim]);
  }
  fprintf(log,"%s:\n",gen->genid);
  
  /* cone list and parameters */
  fprintf(log,"%s: Number of cones = %ld\n",gen->genid, GEN->n_cone);
  fprintf(log,"%s: Cone list:\n",gen->genid);

  /* vertex indices of the cones */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : verteces = {", gen->genid, ic);
    for (d=0; d<dim; d++) {
       fprintf(log,"%ld, ", GEN->cone_list[ic]->index[d]);
    }   
    fprintf(log,"%ld}\n", GEN->cone_list[ic]->index[dim]);
  }
  fprintf(log,"%s:\n",gen->genid);

  /* vertex-vector lengths of the cones */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : norms = {", gen->genid, ic);
    for (d=0; d<dim; d++) {
       fprintf(log,"%g, ", GEN->cone_list[ic]->length[d]);
    }   
    fprintf(log,"%g}\n", GEN->cone_list[ic]->length[dim]);
  }
  fprintf(log,"%s:\n",gen->genid);
  
  /* cone volumes */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : volume = %g  \tsum = %g\n", gen->genid, ic, 
            GEN->cone_list[ic]->volume, GEN->cone_list[ic]->sum_volume);
  }
  fprintf(log,"%s:\n",gen->genid);
  /* cone unit volumes */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : unit_volume = %g\n", gen->genid, ic, 
                                                  GEN->cone_list[ic]->unit_volume);
  }
  fprintf(log,"%s:\n",gen->genid);

  /* touching point of surface and tangential plane */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : spoint = {", gen->genid, ic);
    for (d=0; d<dim; d++) {
       fprintf(log,"%g, ", GEN->cone_list[ic]->spoint[d]);
    }   
    fprintf(log,"%g}\n", GEN->cone_list[ic]->spoint[dim]);
  }
  fprintf(log,"%s:\n",gen->genid);

  /* normal vector to tangential plane */
  for (ic=0; ic<GEN->n_cone; ic++) {
    fprintf(log,"%s:\t%6ld : normal = {", gen->genid, ic);
    for (d=0; d<dim; d++) {
       fprintf(log,"%g, ", GEN->cone_list[ic]->normal[d]);
    }   
    fprintf(log,"%g}\n", GEN->cone_list[ic]->normal[dim]);
  }
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_varou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
