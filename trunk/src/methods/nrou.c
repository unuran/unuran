/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      nrou.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
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
 *   [1] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [1] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (U,V) is uniformly distributed in                                   *
 *       A = {(u,v): 0 < v <= sqrt(g(u/v)), x_0 < u/v < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (U,V) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <utils/fmax_source.h>
#include <utils/unur_fp_source.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "nrou.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define NROU_VARFLAG_VERIFY   0x002u   /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NROU_SET_U       0x001u     /* set u values of bounding rectangle    */
#define NROU_SET_V       0x002u     /* set v values of bounding rectangle    */
#define NROU_SET_CENTER  0x004u     /* set approximate mode of distribution  */
#define NROU_SET_R       0x008u     /* set r-parameter                       */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NROU"         /* type of generator                          */

/* scaling factor of computed minimum boundary rectangle                     */
/* after we have computed the boundary rectangle (0, vmax)x(umin[d], umax[d])*/
/* we scale the obtained boundaries with this factor, i.e. :                 */
/* vmax = vmax * ( 1+ NROU_RECT_SCALING)                                     */
/* umin = umin - (umax-umin)*NROU_RECT_SCALING/2.                            */
/* umax = umax + (umax-umin)*NROU_RECT_SCALING/2.                            */
#define NROU_RECT_SCALING 1e-4


/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_nrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_nrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_nrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_nrou_sample( struct unur_gen *gen );
static double _unur_nrou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_nrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_nrou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.nrou        /* data for parameter object         */
#define GEN       gen->data.nrou        /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define BD_MAX    (DBL_MAX/1000.)       /* constant used for bounding rectangles */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_nrou_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_xmalloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_NROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.umin      = 0.;          /* u-boundary of bounding rectangle (unknown) */
  PAR.umax      = 0.;          /* u-boundary of bounding rectangle (unknown) */
  PAR.vmax      = 0.;          /* v-boundary of bounding rectangle (unknown) */
  PAR.center    = 0.;          /* center of distribution (default: 0)        */
  PAR.r         = 1.;          /* r-parameter of the generalized method      */
  
  par->method   = UNUR_METH_NROU;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_nrou_init;

  return par;

} /* end of unur_nrou_new() */

/*****************************************************************************/

int
unur_nrou_set_u( struct unur_par *par, double umin, double umax )
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (!_unur_FP_greater(umax,umin)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR.umin = umin;
  PAR.umax = umax;

  /* changelog */
  par->set |= NROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_u() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_v( struct unur_par *par, double vmax )
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
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR.vmax = vmax;

  /* changelog */
  par->set |= NROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_v() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_center( struct unur_par *par, double center )
     /*----------------------------------------------------------------------*/
     /* Set the center (approximate mode) of the PDF.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   center ... center of distribution                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* store data */
  PAR.center = center;

  /* changelog */
  par->set |= NROU_SET_CENTER;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_set_center() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* Sets r-parameter of the generalized ratio-of-uniforms method         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   r ... r-parameter                                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0 set to r=1");
    r=1.;
  }
  
  /* store values */
  PAR.r = r;

  /* changelog */
  par->set |= NROU_SET_R;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, NROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | NROU_VARFLAG_VERIFY) : (par->variant & (~NROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, NROU, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= NROU_VARFLAG_VERIFY;
    SAMPLE = _unur_nrou_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~NROU_VARFLAG_VERIFY;
    SAMPLE = _unur_nrou_sample;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_chg_verify() */

/*****************************************************************************/

struct unur_gen *
_unur_nrou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_NROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_nrou_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute bounding rectangle */
  if (_unur_nrou_rectangle(gen)!=UNUR_SUCCESS) {
    free(par); _unur_nrou_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_nrou_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_nrou_init() */

/*---------------------------------------------------------------------------*/

double 
_unur_aux_bound_vmax(double x, void *p) {
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
  struct unur_gen *gen;
  gen=p; /* typecast from void* to unur_gen* */
   
  return pow( _unur_cont_PDF((x),(gen->distr)), 1./(1.+GEN.r) ); 
}

/*---------------------------------------------------------------------------*/

double
_unur_aux_bound_umax(double x, void *p) {
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	
  struct unur_gen *gen;
  gen=p; /* typecast from void* to unur_gen* */
  
  return (x-GEN.center) * pow( _unur_cont_PDF((x),(gen->distr)),
                               GEN.r / (1.+ GEN.r) );
}

/*---------------------------------------------------------------------------*/

double
_unur_aux_bound_umin(double x, void *p) {
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	

  return (- _unur_aux_bound_umax(x,p)) ;
}

/*---------------------------------------------------------------------------*/

int
_unur_nrou_rectangle( struct unur_gen *gen )
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
  struct unur_funct_generic faux; /* function to be minimized/maximized    */
  double x, cx, sx, bx; /* parameters to be used in min/max search */ 

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_NROU_GEN, UNUR_ERR_COOKIE );

  /* Boundary rectangle is already set */
  if ((gen->set & NROU_SET_U) && (gen->set & NROU_SET_V)) {
    return UNUR_SUCCESS;
  }

  /* If center has not been set, we'll set it to the (optional) mode          */
  /* (provided the user has not set any parameters of the bounding rectangle) */
  /* Otherwise the center defaults to the value set in unur_nrou_new() i.e. 0 */
  if (!(gen->set & NROU_SET_CENTER) &&
       (gen->distr->set & UNUR_DISTR_SET_MODE) && 
      !(gen->set & NROU_SET_U) && 
      !(gen->set & NROU_SET_V) ) {
    GEN.center = DISTR.mode; 
  }  
	  
  /* starting point in min/max algorithm */
  cx=GEN.center;  
 
  /* --------------------------------------------------------------------- */
 
  /* calculation of vmax */
  if (!(gen->set & NROU_SET_V)) {
    /* user has not provided any upper bound for v */

    /* we start to check if optional mode is present */
    if ( (gen->distr->set & UNUR_DISTR_SET_MODE) &&
         (DISTR.mode >= DISTR.BD_LEFT) &&
	 (DISTR.mode <= DISTR.BD_RIGHT) ) {
	 
      /* setting vmax to be (f(mode))^(1/(1+r)) */
      GEN.vmax = pow(PDF(DISTR.mode), 1./(1.+GEN.r));
    }
    else {
      /* calculating vmax as maximum of (f(x))^(1/(1+r)) in the domain */
      faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_vmax;
      faux.params = gen; ;
  
      x = _unur_util_find_max(faux, DISTR.BD_LEFT, DISTR.BD_RIGHT, cx);
      if (_unur_FP_is_plusminus_infinity(x)) {
         _unur_error(gen->genid , UNUR_ERR_GENERIC, "Bounding rect (vmax)");  
         return UNUR_ERR_GENERIC;
      }
      GEN.vmax = faux.f(x,faux.params);
    }

    /* additional scaling of boundary rectangle */
    GEN.vmax = GEN.vmax * ( 1+ NROU_RECT_SCALING);
  }

  /* --------------------------------------------------------------------- */
  
  /* calculation of umin and umax */
  if (!(gen->set & NROU_SET_U)) {

    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umin;
    faux.params = gen;

    /* calculating start point for extremum search routine */
    sx = _unur_FP_is_plusminus_infinity(DISTR.BD_LEFT) ? cx-1.: (cx+DISTR.BD_LEFT)/2. ; 
    bx = _unur_FP_is_plusminus_infinity(DISTR.BD_LEFT) ? -BD_MAX: DISTR.BD_LEFT;

    x = (DISTR.BD_LEFT == cx) ? cx: _unur_util_find_max(faux, bx, cx, sx);
          
    while (_unur_FP_is_plusminus_infinity(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       /* _unur_util_find_max() could not yet find a suitable extremum */
       /* trying with a sequence of intervals with decreasing length   */
       bx = bx/10.; sx = bx/2.;  
       x = _unur_util_find_max(faux, bx, cx, sx);
    }
         
    if (_unur_FP_is_plusminus_infinity(x)) {
       /* not able to compute a boundary recangle ...  */ 
       _unur_error(gen->genid , UNUR_ERR_GENERIC, "Bounding rect (umin)");  
       return UNUR_ERR_GENERIC;
    }
    GEN.umin = -faux.f(x,faux.params);

    /* and now, an analogue calculation for umax */

    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umax;
    faux.params = gen;

    /* calculating start point for extremum search routine */
    sx = _unur_FP_is_plusminus_infinity(DISTR.BD_RIGHT) ? cx+1.: (cx+DISTR.BD_RIGHT)/2. ; 
    bx = _unur_FP_is_plusminus_infinity(DISTR.BD_RIGHT) ? BD_MAX: DISTR.BD_RIGHT;

    x = (DISTR.BD_RIGHT == cx) ? cx: _unur_util_find_max(faux, cx, bx, sx);
      
    while (_unur_FP_is_plusminus_infinity(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       /* _unur_util_find_max() could not yet find a suitable extremum */
       /* trying with a sequence of intervals with decreasing length   */
       bx = bx/10.; sx = bx/2.; 
       x = _unur_util_find_max(faux, cx, bx, sx);
    }
           
    if (_unur_FP_is_plusminus_infinity(x)) {
       /* not able to compute a boundary recangle ...  */ 
       _unur_error(gen->genid , UNUR_ERR_GENERIC, "Bounding rect (umax)");  
       return UNUR_ERR_GENERIC;
    }
    GEN.umax = faux.f(x,faux.params);
    
    /* additional scaling of boundary rectangle */
    GEN.umin = GEN.umin - (GEN.umax-GEN.umin)*NROU_RECT_SCALING/2.;
    GEN.umax = GEN.umax + (GEN.umax-GEN.umin)*NROU_RECT_SCALING/2.;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_nrou_rectangle() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_nrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_NROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & NROU_VARFLAG_VERIFY) ? _unur_nrou_sample_check : _unur_nrou_sample;

  gen->destroy = _unur_nrou_free;
  gen->clone = _unur_nrou_clone;

  /* copy some parameters into generator object */
  GEN.umin  = PAR.umin;             /* left u-boundary of bounding rectangle */
  GEN.umax  = PAR.umax;             /* right u-boundary of bounding rectangle */
  GEN.vmax  = PAR.vmax;             /* upper v-boundary of bounding rectangle */
  GEN.center = PAR.center;          /* center of distribution */
  GEN.r = PAR.r;                    /* r-parameter of the generalized rou-method */

  /* initialize parameters */

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_nrou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_nrou_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.nrou

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_nrou_clone() */

/*****************************************************************************/

double
_unur_nrou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( (V = _unur_call_urng(gen->urng)) == 0.);
    V *= GEN.vmax;
    U = GEN.umin + _unur_call_urng(gen->urng) * (GEN.umax - GEN.umin);

    /* compute X */
    X = U/pow(V,GEN.r) + GEN.center;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* accept or reject */
    if (V <= pow(PDF(X), 1./(1.+GEN.r)) )
      return X;
  }

} /* end of _unur_nrou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_nrou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X,fx,sfx,xfx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( (V = _unur_call_urng(gen->urng)) == 0.);
    V *= GEN.vmax;
    U = GEN.umin + _unur_call_urng(gen->urng) * (GEN.umax - GEN.umin);
    
    /* compute x */
    X = U/pow(V,GEN.r) + GEN.center;
    
    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    
    /* evaluate PDF */
    fx = PDF(X);
    
    /* a point on the boundary of the region of acceptance
       has the coordinates ( (X-center) * (fx)^(r/(1+r)), (fx)^(1/(1+r)) ). */
    sfx = pow(fx, 1./(1.+GEN.r));
    xfx = (X-GEN.center) * pow(fx, GEN.r/(1.+GEN.r));
    
    /* check hat */
    if ( ( sfx > (1.+DBL_EPSILON) * GEN.vmax )   /* avoid roundoff error with FP registers */
	 || (xfx < (1.+UNUR_EPSILON) * GEN.umin) 
	 || (xfx > (1.+UNUR_EPSILON) * GEN.umax) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    
    /* accept or reject */
    if (V <= pow(PDF(X), 1./(1.+GEN.r)) )
      return X;
  }

} /* end of _unur_nrou_sample_check() */

/*****************************************************************************/

void
_unur_nrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_nrou_free() */

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
_unur_nrou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = nrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_nrou_sample",gen->genid);
  if (gen->variant & NROU_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(log,"%s: r-parameter = %g",gen->genid, GEN.r);
  _unur_print_if_default(gen,NROU_SET_R);
  fprintf(log,"\n%s:\n",gen->genid);

  /* center */
  fprintf(log,"%s: center = %g\n",gen->genid,GEN.center);
  fprintf(log,"%s:\n",gen->genid);

  /* bounding rectangle */
  fprintf(log,"%s: Rectangle:\n",gen->genid);
  fprintf(log,"%s:    left  upper point = (%g,%g)\n",gen->genid,GEN.umin,GEN.vmax);
  fprintf(log,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN.umax,GEN.vmax);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_nrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
