/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vnrou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    naive ratio-of-uniforms method                               *
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
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
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
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "vnrou.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define VNROU_VARFLAG_VERIFY   0x002u   /* run verify mode                   */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define VNROU_SET_U       0x001u     /* set u values of bounding rectangle   */
#define VNROU_SET_V       0x002u     /* set v values of bounding rectangle   */
#define VNROU_SET_R       0x008u     /* set r-parameter                      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VNROU"         /* type of generator                         */

/*---------------------------------------------------------------------------*/
/* convergence parameters for the hooke optimization algorithm */

#define VNROU_HOOKE_RHO 0.5 
#define VNROU_HOOKE_EPSILON 1e-7 
#define VNROU_HOOKE_MAXITER 10000

/* scaling factor of computed minimum boundary rectangle */
/* after we have computed the boundary rectangle (0, vmax)x(umin[d], umax[d])*/
/* we scale the obtained boundaries with this factor, i.e. :                 */
/* vmax = vmax * ( 1+ VNROU_RECT_SCALING)                                    */
/* umin[d] = umin[d] - (umax[d]-umin[d])*VNROU_RECT_SCALING/2.               */
/* umax[d] = umax[d] + (umax[d]-umin[d])*VNROU_RECT_SCALING/2.               */
#define VNROU_RECT_SCALING 1e-4

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vnrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vnrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vnrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static void  _unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec );
static void  _unur_vnrou_sample_check( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_vnrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_vnrou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       par->data.vnrou        /* data for parameter object        */
#define GEN       gen->data.vnrou        /* data for generator object        */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */     
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_vnrou_new( const struct unur_distr *distr )
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
  COOKIE_SET(par,CK_VNROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* copy number of dimensions from the distribution object */
  PAR.dim = distr->dim;


  /* set default values */
  PAR.r		= 1.; 	      /* r-parameter of the generalized method       */
  PAR.vmax      = 0.;         /* v-boundary of bounding rectangle (unknown)  */
  PAR.umin 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  PAR.umax 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  par->method   = UNUR_METH_VNROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_vnrou_init;

  return par;

} /* end of unur_vnrou_new() */

/*****************************************************************************/

int
unur_vnrou_set_u( struct unur_par *par, double *umin, double *umax )
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
  _unur_check_par_object( par, VNROU );
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
  par->set |= VNROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_u() */


/*---------------------------------------------------------------------------*/


int
unur_vnrou_set_v( struct unur_par *par, double vmax )
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
  _unur_check_par_object( par, VNROU );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR.vmax = vmax;

  /* changelog */
  par->set |= VNROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_v() */


/*---------------------------------------------------------------------------*/

int
unur_vnrou_set_r( struct unur_par *par, double r )
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
  _unur_check_par_object( par, VNROU );
  if (r <= 0.) {
    _unur_warning(GENTYPE, UNUR_ERR_PAR_SET, "r<=0; set to r=1");
    r = 1.;
  }
  /* store data */
  PAR.r = r;

  /* changelog */
  par->set |= VNROU_SET_R;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, VNROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | VNROU_VARFLAG_VERIFY) : (par->variant & (~VNROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify bounding rectangle on */
    gen->variant |= VNROU_VARFLAG_VERIFY;
    SAMPLE = _unur_vnrou_sample_check;
  }
  else {
    /* turn verify bounding rectangle off */
    gen->variant &= ~VNROU_VARFLAG_VERIFY;
    SAMPLE = _unur_vnrou_sample_cvec;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_chg_verify() */

/*****************************************************************************/

struct unur_gen *
_unur_vnrou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_VNROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_vnrou_create(par);
  if (!gen) { 
    free(par); 
    return NULL; 
  }

  /* compute bounding rectangle */
  if (_unur_vnrou_rectangle(gen)!=UNUR_SUCCESS) {
    free(par); _unur_vnrou_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_vnrou_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_vnrou_init() */

/*---------------------------------------------------------------------------*/

double 
_unur_vnrou_aux_vmax(double *x, void *p ) 
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  
  gen = p; /* typecast from void* to unur_gen* */
  return -pow( _unur_cvec_PDF((x),(gen->distr)) , 
                              1./(1.+ GEN.r * GEN.dim) ); 
}

/*---------------------------------------------------------------------------*/

double
_unur_vnrou_aux_umin(double *x, void *p) 
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	
{
  struct unur_gen *gen;
  
  gen = p; /* typecast from void* to unur_gen* */
  return (x[GEN.aux_dim] - GEN.center[GEN.aux_dim]) 
         * pow( _unur_cvec_PDF((x),(gen->distr)), 
                GEN.r / (1.+ GEN.r * GEN.dim) );
}

/*---------------------------------------------------------------------------*/

double
_unur_vnrou_aux_umax(double *x, void *p) 
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	
{
  return (- _unur_vnrou_aux_umin(x,p)) ;
}

/*---------------------------------------------------------------------------*/

int
_unur_vnrou_rectangle( struct unur_gen *gen )
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

  struct unur_funct_vgeneric faux; /* function to be minimized/maximized    */
  double *xstart, *xend, *xumin, *xumax; /* coordinate arrays used in maximum/minimum calculations */
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  int hooke_iters_vmax;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umin;  /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umax;  /* actual number of min/max iterations = return value of hooke()*/
  double scaled_epsilon; /* to be used in the hooke algorithm */

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_VNROU_GEN, UNUR_ERR_COOKIE );

  /* Boundary rectangle is already set */
  if ((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)) {
    return UNUR_SUCCESS;
  }

  /* dimension of the distribution */
  dim = GEN.dim;

  /* allocate memory for the coordinate vectors */
  xstart = _unur_xmalloc(dim * sizeof(double));
  xend   = _unur_xmalloc(dim * sizeof(double));
  xumin  = _unur_xmalloc(dim * sizeof(double));
  xumax  = _unur_xmalloc(dim * sizeof(double));
  
  /* calculation of vmax */
  if (!(gen->set & VNROU_SET_V)) {
    /* user has not provided any upper bound for v */
    /* calculating vmax as maximum of sqrt(f(x)) in the domain */
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_vnrou_aux_vmax;
      faux.params = gen;

      memcpy(xstart, GEN.center, dim * sizeof(double)); 
      hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend, 
				      VNROU_HOOKE_RHO, VNROU_HOOKE_EPSILON, VNROU_HOOKE_MAXITER);

      GEN.vmax = -faux.f(xend, faux.params);
      
      if (hooke_iters_vmax >= VNROU_HOOKE_MAXITER) {
	 scaled_epsilon = VNROU_HOOKE_EPSILON * GEN.vmax;
	 if (scaled_epsilon>VNROU_HOOKE_EPSILON) scaled_epsilon=VNROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         memcpy(xstart, xend, dim * sizeof(double)); 
         hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend, 
                              VNROU_HOOKE_RHO, scaled_epsilon , VNROU_HOOKE_MAXITER);
         GEN.vmax = -faux.f(xend, faux.params);
         if (hooke_iters_vmax >= VNROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (vmax)");  
         }
      }

  }

  /* calculation of umin and umax */
  if (!(gen->set & VNROU_SET_U)) {
    
    for (d=0; d<dim; d++) {

      /* setting coordinate dimension to be used by the auxiliary functions */
      GEN.aux_dim  = d;
  
      /*-----------------------------------------------------------------------------*/
      /* calculation for umin */
      
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_vnrou_aux_umin;
      faux.params = gen;

      hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend, 
                           VNROU_HOOKE_RHO, VNROU_HOOKE_EPSILON, VNROU_HOOKE_MAXITER);
      GEN.umin[d] = faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumin, xend, dim * sizeof(double)); 

      /*-----------------------------------------------------------------------------*/
      /* and now, an analogue calculation for umax */

      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_vnrou_aux_umax;
      faux.params = gen;

      hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend, 
                           VNROU_HOOKE_RHO, VNROU_HOOKE_EPSILON, VNROU_HOOKE_MAXITER);
      GEN.umax[d] = -faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumax, xend, dim * sizeof(double)); 

      /*-----------------------------------------------------------------------------*/
      /* checking if we need to recalculate umin */
      if (hooke_iters_umin >= VNROU_HOOKE_MAXITER) {
	 scaled_epsilon = VNROU_HOOKE_EPSILON * (GEN.umax[d]-GEN.umin[d]);
	 if (scaled_epsilon>VNROU_HOOKE_EPSILON) scaled_epsilon=VNROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         memcpy(xstart, xumin, dim * sizeof(double)); 
         hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend, 
                              VNROU_HOOKE_RHO, scaled_epsilon , VNROU_HOOKE_MAXITER);
         GEN.umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umin >= VNROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umin)");  
         }
      }

      /* checking if we need to recalculate umax */
      if (hooke_iters_umax >= VNROU_HOOKE_MAXITER) {
	 scaled_epsilon = VNROU_HOOKE_EPSILON * (GEN.umax[d]-GEN.umin[d]);
	 if (scaled_epsilon>VNROU_HOOKE_EPSILON) scaled_epsilon=VNROU_HOOKE_EPSILON;

         /* recalculating extremum with scaled_epsilon and new starting point */
         memcpy(xstart, xumax, dim * sizeof(double)); 
         hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend, 
                              VNROU_HOOKE_RHO, scaled_epsilon , VNROU_HOOKE_MAXITER);
         GEN.umin[d] = faux.f(xend, faux.params);
         if (hooke_iters_umax >= VNROU_HOOKE_MAXITER) {
           _unur_warning(gen->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umax)");  
         }
      }
    
      /*-----------------------------------------------------------------------------*/
      /* additional scaling of boundary rectangle */   
      GEN.vmax = GEN.vmax * ( 1+ VNROU_RECT_SCALING);
      GEN.umin[d] = GEN.umin[d] - (GEN.umax[d]-GEN.umin[d])*VNROU_RECT_SCALING/2.;
      GEN.umax[d] = GEN.umax[d] + (GEN.umax[d]-GEN.umin[d])*VNROU_RECT_SCALING/2.;        
    
    }
  }


  free(xstart); free(xend); free(xumin); free(xumax);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_vnrou_rectangle() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vnrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_VNROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & VNROU_VARFLAG_VERIFY) ? _unur_vnrou_sample_check : _unur_vnrou_sample_cvec;

  gen->destroy = _unur_vnrou_free;
  gen->clone = _unur_vnrou_clone;

  /* allocate memory for u-boundary arrays */
  GEN.umin = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */
  GEN.umax = _unur_xmalloc( PAR.dim * sizeof(double)); /* bounding rectangle */

  /* copy parameters into generator object */
  GEN.dim   = PAR.dim;              /* dimension */
  GEN.r     = PAR.r;                /* r-parameter of the vnrou method */  
  GEN.vmax  = PAR.vmax;             /* upper v-boundary of bounding rectangle */
  
  if (PAR.umin != NULL) memcpy(GEN.umin, PAR.umin, GEN.dim * sizeof(double));
  if (PAR.umax != NULL) memcpy(GEN.umax, PAR.umax, GEN.dim * sizeof(double));
  /* get center of the distribution */
  GEN.center = unur_distr_cvec_get_center(gen->distr);
 
  /* initialize parameters */

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_vnrou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vnrou_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.vnrou

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VNROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory for u-arrays and center */
  CLONE.umin = _unur_xmalloc( GEN.dim * sizeof(double));
  CLONE.umax = _unur_xmalloc( GEN.dim * sizeof(double));
  
  /* copy parameters into clone object */
  memcpy(CLONE.umin, GEN.umin, GEN.dim * sizeof(double));
  memcpy(CLONE.umax, GEN.umax, GEN.dim * sizeof(double));

  /* copy data */
  CLONE.center = unur_distr_cvec_get_center(clone->distr);

  return clone;

#undef CLONE
} /* end of _unur_vnrou_clone() */

/*****************************************************************************/

void
_unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{ 
  double U, V;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID); 

  dim = GEN.dim;
 
  while (1) {

    /* generate point uniformly on rectangle */
    while ( (V = _unur_call_urng(gen->urng)) == 0.);
    V *= GEN.vmax;
    for (d=0; d<dim; d++) {
      U = GEN.umin[d] + _unur_call_urng(gen->urng) * (GEN.umax[d] - GEN.umin[d]);
      vec[d] = U/pow(V,GEN.r) + GEN.center[d];
    }
    
    /* X[] inside domain ? */
    
    /* accept or reject */
    if (V <= pow(PDF(vec),1./(GEN.r * dim + 1.)))
      return;
  }

} /* end of _unur_vnrou_sample() */

/*---------------------------------------------------------------------------*/

void
_unur_vnrou_sample_check( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random sample vector (return)                              */
     /*----------------------------------------------------------------------*/
{ 
  double U, V;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  int hat_error;

  double fx,sfx,xfx;
  
  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID); 

  dim = GEN.dim;
 
  while (1) {
    /* generate point uniformly on rectangle */
    while ( (V = _unur_call_urng(gen->urng)) == 0.);
    V *= GEN.vmax;
    for (d=0; d<dim; d++) {
      U = GEN.umin[d] + _unur_call_urng(gen->urng) * (GEN.umax[d] - GEN.umin[d]);
      vec[d] = U/pow(V,GEN.r) + GEN.center[d];
    }
    
    /* X[] inside domain ? */

    /* evaluate PDF */
    fx = PDF(vec);
    
    /* a point on the boundary of the region of acceptance
       has the coordinates ( (vec[]-center[]) * (fx)^(r/r*dim+1)), fx^(1/r*dim+1) ). */
    sfx = pow( fx, 1./(GEN.r * dim+1.) );
    /* check hat */
    hat_error=0;
    if ( sfx > (1.+DBL_EPSILON) * GEN.vmax ) hat_error++;  
   
    sfx = pow( fx, GEN.r/(GEN.r * dim + 1.) );
    for (d=0; d<dim; d++) {
     xfx = (vec[d]-GEN.center[d]) * sfx;
     if ( (xfx < (1.+UNUR_EPSILON) * GEN.umin[d]) 
       || (xfx > (1.+UNUR_EPSILON) * GEN.umax[d]))
       hat_error++;
    }

    if (hat_error>0) _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
 
    /* accept or reject */
    if (V <= pow(PDF(vec),1./( GEN.r * dim + 1.)))
      return;
  }

} /* end of _unur_vnrou_sample_check() */

/*****************************************************************************/

void
_unur_vnrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VNROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN.umin) free(GEN.umin); 
  if (GEN.umax) free(GEN.umax);
  _unur_generic_free(gen);

} /* end of _unur_vnrou_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_vnrou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  
  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);

  log = unur_get_stream();
  dim = GEN.dim;

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = vnrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_vnrou_sample",gen->genid);
  if (gen->variant & VNROU_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(log,"%s: r-parameter = %g",gen->genid, GEN.r);
  _unur_print_if_default(gen,VNROU_SET_R);
  fprintf(log,"\n%s:\n",gen->genid);

  /* print center */
  _unur_matrix_print_vector( GEN.dim, GEN.center, "center =", log, gen->genid, "\t   ");

  /* print bounding rectangle */
  fprintf(log,"%s: Rectangle:",gen->genid);
  if (!((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)))
    fprintf(log,"\t[computed]");
  else 
    fprintf(log,"\t[input]");
  fprintf(log,"\n");

  fprintf(log,"%s:\tvmax = %g\n",gen->genid, GEN.vmax);
  for (d=0; d<dim; d++) {
    fprintf(log,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
	    d, d, GEN.umin[d], GEN.umax[d]);
  }
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_vnrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
