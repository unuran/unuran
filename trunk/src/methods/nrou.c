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
#include <utils/fminmax_source.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "nrou.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define NROU_VARFLAG_VERIFY   0x002u   /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NROU_SET_U       0x001u     /* set u values of bounding rectangle    */
#define NROU_SET_V       0x002u     /* set v values of bounding rectangle    */
#define NROU_SET_CENTER  0x004u     /* set approximate mode of distribution  */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NROU"         /* type of generator                          */

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

struct unur_gen *_gen; /* generator object for bounding rect calculations */


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.nrou        /* data for parameter object         */
#define GEN       gen->data.nrou        /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

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
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_NROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.umin      = 0.;          /* u-boundary of bounding rectangle (unknown) */
  PAR.umax      = 0.;          /* u-boundary of bounding rectangle (unknown) */
  PAR.vmax      = 0.;          /* v-boundary of bounding rectangle (unknown) */
  PAR.center    = 0.;          /* center of distribution (default: 0)        */

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
unur_nrou_set_rect_u( struct unur_par *par, double umin, double umax )
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
unur_nrou_set_rect_v( struct unur_par *par, double vmax )
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
_unur_aux_bound_vmax(double x, double *p) {
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
   
  return sqrt( _unur_cont_PDF((x),(_gen->distr)) ); 
}

/*---------------------------------------------------------------------------*/

double
_unur_aux_bound_umax(double x, double *p) {
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	

  return (x-p[0]) * sqrt( _unur_cont_PDF((x),(_gen->distr)) );
}

/*---------------------------------------------------------------------------*/

double
_unur_aux_bound_umin(double x, double *p) {
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
  double p[1]; /* parameter for auxiliary functions */
  double x;

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
	  
  /* copy of generator to be used in auxiliary functions */
  _gen = gen;

  /* parameter to be used in auxiliary functions */
  p[0]=GEN.center;
  
  /* calculation of vmax */
  if (!(gen->set & NROU_SET_V)) {
    /* user has not provided any upper bound for v */

    /* we start to check if optional mode is present */
    if ( (gen->distr->set & UNUR_DISTR_SET_MODE) &&
         (DISTR.mode >= DISTR.BD_LEFT) &&
	 (DISTR.mode <= DISTR.BD_RIGHT) ) {
	 
      /* setting vmax to be sqrt(f(mode)) */
      GEN.vmax = sqrt(PDF(DISTR.mode));
    }
    else {
      /* calculating vmax as maximum of sqrt(f(x)) in the domain */
      faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_vmax;
      faux.params = NULL;
  
      x = _unur_util_find_max(faux, DISTR.BD_LEFT, DISTR.BD_RIGHT, p[0]);
      GEN.vmax = faux.f(x,p);
    }
  }

  /* calculation of umin and umax */
  if (!(gen->set & NROU_SET_U)) {
    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umin;
    faux.params = p;

    x = _unur_util_find_max(faux, DISTR.BD_LEFT, p[0], p[0]);
    GEN.umin = -faux.f(x,p);

    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umax;
    faux.params = p;

    x = _unur_util_find_max(faux, p[0], DISTR.BD_RIGHT, p[0]);
    GEN.umax = faux.f(x,p);
  }

  /* TODO : check for (umin,umax) sanity 
            e.g. cauchy-like distributions */


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
    X = U/V + GEN.center;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* accept or reject */
    if (U*U <= PDF(X))
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
    X = U/V + GEN.center;
    
    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    
    /* evaluate PDF */
    fx = PDF(X);
    
    /* the point on the boundary of the region of acceptance
       in direction X = V/U has the coordinates
       ( X * sqrt(fx), sqrt(fx) ). */
    sfx = sqrt(fx);
    xfx = (X-GEN.center) * sfx;
    
    /* check hat */
    if ( ( sfx > (1.+DBL_EPSILON) * GEN.vmax )   /* avoid roundoff error with FP registers */
	 || (xfx < (1.+UNUR_EPSILON) * GEN.umin) 
	 || (xfx > (1.+UNUR_EPSILON) * GEN.umax) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    
    /* accept or reject */
    if (U*U <= PDF(X))
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
  _unur_distr_free(gen->distr);
  _unur_free_genid(gen);

  COOKIE_CLEAR(gen);
  free(gen);

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

  fprintf(log,"%s: center = %g\n",gen->genid,GEN.center);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: Rectangle:\n",gen->genid);
  fprintf(log,"%s:    left  upper point = (%g,%g)\n",gen->genid,GEN.umin,GEN.vmax);
  fprintf(log,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN.umax,GEN.vmax);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_nrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
