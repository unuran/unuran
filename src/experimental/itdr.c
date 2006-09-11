/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      itdr.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    inverse transformed density rejection                        *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and pole.                                                  *
 *      Produce a value x consistent with its density.                       *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      location of pole                                                     *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      splitting point between pole and tail region                         *
 *      c-value for pole and tail region, repectively                        *
 *                                                                           *
 *****************************************************************************
     $Id: nrou.c 2725 2005-12-23 10:52:19Z leydold $
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
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
/* #include <distr/cont.h> */
/* #include <utils/fmax_source.h> */
#include <utils/unur_fp_source.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "itdr.h"
#include "itdr_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants:                                                                */

/* maximum value for parameter c                                             */
#define C_MAX  (-0.1)

/* relative distance for computing derivatives numerically                   */
/* (1.e-6 yieldes better results than 1.e-8 when compared with exact values) */
#define DX (1.e-6)      

/* point near pole for estimating lim_{x->plole} ilc(x)                      */
/* the point x_i * NEAR_POLE is used                                         */
#define NEAR_POLE  (1.e-8)

/* point near pole for testing validitiy of hat function                     */
#define TEST_NEAR_POLE  (1.e-100)

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define ITDR_VARFLAG_VERIFY   0x001u   /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define ITDR_SET_BX      0x001u     /* set splitting point                   */
#define ITDR_SET_CP      0x002u     /* set c-value for pole region           */
#define ITDR_SET_CT      0x004u     /* set c-value for tail region           */

/*---------------------------------------------------------------------------*/

#define GENTYPE "ITDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_sample( struct unur_gen *gen );
static double _unur_itdr_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static void _unur_itdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_itdr_get_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* construct hat function                                                    */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_lc( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* compute local concavity at x                                              */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_ilc( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* compute inverse local concavity at x                                      */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_find_xt( struct unur_gen *gen, double b );
/*---------------------------------------------------------------------------*/
/* solves equation (x-b)*f'(x)+f(x)=0, where f is PDF of distribution        */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_itdr_debug_init( const struct unur_gen *gen, int error );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_itdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_itdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */
#define dPDF(x)   _unur_cont_dPDF((x),(gen->distr))   /* call to derivative of PDF */

/*---------------------------------------------------------------------------*/
/* transformations */

#define T(c,x)   ( -1./pow((x),-(c)) )
#define DT(c,x)  ( -(c)/pow((x),-((c)-1.)) )
#define TI(c,x)  ( 1./pow(-(x),-1./((c))) )
#define FT(c,x)  ( -1./pow(-(x),-((c)+1.)/(c))*((c)/((c)+1.)) )
#define FTI(c,x) ( -1./pow(-(x)*((c)+1.)/(c),-(c)/((c)+1.)) )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_itdr_new( const struct unur_distr *distr )
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

  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode (pole)");
    return NULL; 
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_itdr_par) );
  COOKIE_SET(par,CK_ITDR_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */
  PAR->bx = INFINITY;       /* splitting point betw. pole and tail (unknown) */
  PAR->cp = INFINITY;       /* c-value for pole region (unknown)             */
  PAR->ct = INFINITY;       /* c-value for tail region (unknown)             */
  
  par->method   = UNUR_METH_ITDR;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_itdr_init;

  return par;

} /* end of unur_itdr_new() */

/*****************************************************************************/

int
unur_itdr_set_bx( struct unur_par *par, double bx )
     /*----------------------------------------------------------------------*/
     /* Sets splitting point bx between pole and tail region.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   bx ... splitting point                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  if (! (_unur_FP_greater(bx, par->distr->data.cont.BD_LEFT)) || 
      ! (_unur_FP_less(   bx, par->distr->data.cont.BD_RIGHT)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"bx out of domain");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->bx = bx;

  /* changelog */
  par->set |= ITDR_SET_BX;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_bx() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_cp( struct unur_par *par, double cp )
     /*----------------------------------------------------------------------*/
     /* Sets c-value for transformation T for inverse density in pole region */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cp ... c-value for pole region                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  if ( cp > C_MAX || cp <= -1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->cp = cp;

  /* changelog */
  par->set |= ITDR_SET_CP;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_cp() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_ct( struct unur_par *par, double ct )
     /*----------------------------------------------------------------------*/
     /* Sets c-value for transformation T for density in tail region         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ct ... c-value for tail region                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double range;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  range = ( par->distr->data.cont.BD_RIGHT
	    - par->distr->data.cont.BD_LEFT );
  if ( ct > C_MAX || (ct <= -1. && _unur_isfinite(range)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ct > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->ct = ct;

  /* changelog */
  par->set |= ITDR_SET_CT;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_ct() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, ITDR );

  /* we use a bit in variant */
  par->variant = (verify) 
    ? (par->variant | ITDR_VARFLAG_VERIFY) 
    : (par->variant & (~ITDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_itdr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, ITDR, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= ITDR_VARFLAG_VERIFY;
    SAMPLE = _unur_itdr_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~ITDR_VARFLAG_VERIFY;
    SAMPLE = _unur_itdr_sample;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_itdr_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_itdr_init( struct unur_par *par )
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
  int error;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_ITDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_itdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* create hat function */
  error = _unur_itdr_get_hat(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_itdr_debug_init(gen,error);
#endif
  
  if (error != UNUR_SUCCESS) {
    _unur_itdr_free(gen); return NULL;
  }

  return gen;

} /* end of _unur_itdr_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_itdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_itdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_ITDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & ITDR_VARFLAG_VERIFY) ? _unur_itdr_sample_check : _unur_itdr_sample;

  gen->destroy = _unur_itdr_free;
  gen->clone = _unur_itdr_clone;

  /* copy some parameters into generator object */
  GEN->bx = PAR->bx;       /* splitting point betw. pole and tail */
  GEN->cp = PAR->cp;       /* c-value for pole region             */
  GEN->ct = PAR->ct;       /* c-value for tail region             */

  /* initialize values */
  GEN->xp = INFINITY;      /* design point in pole region         */
  GEN->xt = INFINITY;      /* design point in tail region         */
  GEN->alphap = INFINITY;  /* parameters for hat in pole region   */
  GEN->betap = INFINITY;
  GEN->Tfxt = INFINITY;    /* parameters for hat in tail region   */
  GEN->dTfxt = INFINITY;   /* parameters for hat in tail region   */
  GEN->by = INFINITY;      /* hat of pole region at bx            */
  GEN->Ap = INFINITY;      /* areas in upper pole region          */     
  GEN->Ac = INFINITY;      /* areas in central region             */     
  GEN->At = INFINITY;      /* areas in tail region                */     
  GEN->Atot = INFINITY;    /* total area below hat                */
  
  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_itdr_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_itdr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_itdr_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_itdr_clone() */

/*****************************************************************************/

double
_unur_itdr_sample( struct unur_gen *gen )
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
  double U, X, Y;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on (0,Atot) */
    U = _unur_call_urng(gen->urng) * GEN->Atot;
  
    /* generate pair (X,Y) below hat */
    if (U < GEN->Ap) {
      /* upper pole region */
      Y = ( FTI(GEN->cp, GEN->betap*U + FT(GEN->cp,GEN->alphap+GEN->betap*GEN->by)) 
	    - GEN->alphap ) / GEN->betap;
      X = _unur_call_urng(gen->urng) * TI(GEN->cp, GEN->alphap+GEN->betap*Y);
    }
    else if ((U -= GEN->Ap) < GEN->Ac) {
      /* central region */
      X = U * GEN->bx / GEN->Ac;
      Y = _unur_call_urng(gen->urng) * GEN->by;
    }
    else {
      /* tail region */
      U -= GEN->Ac;
      X = GEN->xt + (FTI(GEN->ct, 
			 GEN->dTfxt*U 
			 + FT(GEN->ct, 
			      GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			 )
		     - GEN->Tfxt) / GEN->dTfxt;
      Y = ( _unur_call_urng(gen->urng) 
	    * TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
    }

    /* accept or reject */
    if (Y <= PDF(X))
      return X;

  }

} /* end of _unur_itdr_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_sample_check( struct unur_gen *gen )
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
/*   double U,V,X,fx,sfx,xfx; */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);

/*   while (1) { */
/*     /\* generate point uniformly on rectangle *\/ */
/*     while ( (V = _unur_call_urng(gen->urng)) == 0.); */
/*     V *= GEN->vmax; */
/*     U = GEN->umin + _unur_call_urng(gen->urng) * (GEN->umax - GEN->umin); */
    
/*     /\* compute x *\/ */
/*     if (GEN->r == 1.) X = U/V + GEN->center; */
/*     else             X = U/pow(V,GEN->r) + GEN->center; */
    
/*     /\* inside domain ? *\/ */
/*     if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) ) */
/*       continue; */
    
/*     /\* evaluate PDF *\/ */
/*     fx = PDF(X); */
    
/*     /\* a point on the boundary of the region of acceptance */
/*        has the coordinates ( (X-center) * (fx)^(r/(1+r)), (fx)^(1/(1+r)) ). *\/ */
/*     if (GEN->r == 1.) { */
/*       /\* normal rou-method with square-root *\/ */
/*       sfx = sqrt(fx); */
/*       xfx = (X-GEN->center) * sfx; */
/*     } */
/*     else { */
/*       /\* generalized rou-method with pow-function *\/ */
/*       sfx = pow(fx, 1./(1.+GEN->r)); */
/*       xfx = (X-GEN->center) * pow(fx, GEN->r/(1.+GEN->r)); */
/*     } */
    
/*     /\* check hat *\/ */
/*     if ( ( sfx > (1.+DBL_EPSILON) * GEN->vmax )   /\* avoid roundoff error with FP registers *\/ */
/* 	 || (xfx < (1.+UNUR_EPSILON) * GEN->umin)  */
/* 	 || (xfx > (1.+UNUR_EPSILON) * GEN->umax) ) */
/*       _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)"); */
    
/*     /\* accept or reject *\/ */
/*     if (GEN->r ==1) { */
/*       /\* normal rou-method with square-root *\/ */
/*       if (V*V <= PDF(X)) */
/*         return X; */
/*     } */
/*     else { */
/*       /\* generalized rou-method with pow-function *\/ */
/*       if (V <= pow(PDF(X), 1./(1.+GEN->r)) ) */
/*         return X; */
/*     } */
/*   } */

  return 0.;

} /* end of _unur_itdr_sample_check() */

/*****************************************************************************/

void
_unur_itdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_ITDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_itdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_itdr_get_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* construct hat function                                               */
     /* compute local concavity at x                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define hp(x)  ( (T(cp,(x)) - GEN->alphap) / GEN->betap )
#define ht(x)  ( TI(ct, GEN->Tfxt + GEN->dTfxt*((x)-xt)) )

  double bx, cp, xp, ct, xt;
  double ilc_near_pole, ilc_bx;
  double lc_bx;
  double near_pole;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,UNUR_ERR_COOKIE);

  /* -------------------- */
  /* Get candidate for bx */

  if ( !(gen->set & ITDR_SET_BX) ) {
    GEN->bx = bx = _unur_itdr_find_xt( gen, 0. );
    if (!_unur_isfinite(bx)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute bx");
      return UNUR_ERR_DISTR_PROP;
    }
  }
  else 
    bx = GEN->bx;

  /* -------------------- */
  /* pole region          */

  /* get cp */
  if ( !(gen->set & ITDR_SET_CP) ) {
    near_pole = bx * NEAR_POLE;
    ilc_near_pole = log(PDF(near_pole)) / log(near_pole);
    ilc_bx = _unur_itdr_ilc(gen, bx);
    cp = min(ilc_near_pole,ilc_bx);
    if (cp > C_MAX) cp = C_MAX;
    
    if (!(gen->set & ITDR_SET_BX) && cp < -0.5) {
      bx *= 2.;
      ilc_bx = _unur_itdr_ilc(gen, bx);
      if (cp > ilc_bx) cp = ilc_bx;
    }

    if (cp <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->cp = cp;
  }
  else
    cp = GEN->cp;

  /* get design point xp */
  while (1) {
    xp = bx * pow(1.+cp, -1./cp);
    if ( !(xp > 0. && xp < bx) ) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: xp");
      return UNUR_ERR_DISTR_PROP;
    }

    GEN->betap = DT(cp,xp) / dPDF(xp);
    GEN->alphap = T(cp,xp) - GEN->betap * PDF(xp);

    if ( hp(TEST_NEAR_POLE) < PDF(TEST_NEAR_POLE) ||
	 hp(bx) < PDF(bx) ) {
      if (gen->set & ITDR_SET_CP) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"inverse pdf not T_cp concave");
	return UNUR_ERR_DISTR_PROP;
      }
      GEN->cp = cp = 0.9 * cp -0.1;
      if (cp < -0.999) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
	return UNUR_ERR_DISTR_PROP;
      }	
    }
    else break;
  }
  GEN->xp = xp;

  /* -------------------- */
  /* tail region          */

  GEN->xt = xt = _unur_itdr_find_xt( gen, bx );

  if ( !(gen->set & ITDR_SET_CT) ) {
    ct = _unur_itdr_lc(gen, 0.5*(bx + xt));
    if (ct > C_MAX) ct = C_MAX;
    if (ct <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: ct");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->ct = ct;
  }
  else
    ct = GEN->ct;

  /* get design point xp */
  lc_bx = _unur_itdr_lc(gen, bx);
  while (1) {
    GEN->Tfxt = T(ct, PDF(xt));
    GEN->dTfxt = DT(ct, PDF(xt)) * dPDF(xt);

    if ( ht(1000.*bx) < PDF(1000.*bx) ||
	 ht(bx) < PDF(bx) ) {
      if (gen->set & ITDR_SET_CT) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"pdf not T_ct concave");
	return UNUR_ERR_DISTR_PROP;
      }
      GEN->ct = ct = 0.5*(ct + lc_bx);
      if (ct < -0.999 || _unur_FP_approx(ct,lc_bx)) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for tail: ct");
	return UNUR_ERR_DISTR_PROP;
      }
    }
    else
      break;
  }

  /* -------------------- */
  /* parameters           */

  GEN->by = hp(bx);
  GEN->Ap = -FT(cp, GEN->alphap + GEN->betap * GEN->by) / GEN->betap;
  GEN->Ac = GEN->by * bx;
  GEN->At = -FT(ct, GEN->Tfxt + GEN->dTfxt * (bx-xt)) / GEN->dTfxt;
  GEN->Atot = GEN->Ap + GEN->Ac + GEN->At;

  /* -------------------- */

  return UNUR_SUCCESS;

#undef hp
#undef ht
} /* end of _unur_itdr_get_hat() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_lc( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* compute local concavity at x                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point x                                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   local concavity                                                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double dx, f, df, ddf;

  dx = x * DX;
  if (x-dx < DISTR.BD_LEFT) dx = x - DISTR.BD_LEFT;
  if (x+dx > DISTR.BD_RIGHT) dx = DISTR.BD_RIGHT - x;

  f = PDF(x);
  df = dPDF(x);
  ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);

  return 1. - ddf*f/(df*df); 
} /* end of _unur_itdr_lc() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_ilc( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* compute inverse local concavity at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point x                                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   inverse local concavity                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double dx, df, ddf;

  dx = x * DX;
  if (x-dx < DISTR.BD_LEFT) dx = x - DISTR.BD_LEFT;
  if (x+dx > DISTR.BD_RIGHT) dx = DISTR.BD_RIGHT - x;

  df = dPDF(x);
  ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);

  return 1.+x*ddf/(df); 
} /* end of _unur_itdr_ilc() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_find_xt( struct unur_gen *gen, double b )
     /*----------------------------------------------------------------------*/
     /* solves equation (x-b)*f'(x)+f(x)=0, where f is PDF of distribution   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   startx ... starting point for finding x                            */
     /*                                                                      */
     /* return:                                                              */
     /*   solution xi                                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
#define FKT(x)  (((x)-b)*dPDF(x) + PDF(x))   /* function for finding root */

  double xl, xu;  /* lower and upper boundary of bracket */
  double xn;      /* new guess for root */

  /* check parameter */
  if (b < 0) return INFINITY;

  /* find bracket for root */
  /* l.h.s. boundary */
  xl = (b>0.) ? b * (1.+DBL_EPSILON) : DBL_MIN;
  while (!_unur_isfinite(FKT(xl))) {
    xl += xl - b;
    if (!_unur_isfinite(xl)) return INFINITY;
  }
  if (FKT(xl)<0.) return INFINITY;

  /* r.h.s. boundary */
  xu = max(2.*b, xl);
  while(FKT(xu) > 0.){
    xl = xu;
    xu += xu - b;
    if (!_unur_isfinite(xu)) return INFINITY;
  }

  /* use bisection to find root */
  while(!_unur_FP_approx(xl,xu)){
    xn = 0.5*(xl+xu);
    if(FKT(xn)>0.) 
      xl = xn;
    else 
      xu = xn;
  }

  /* return point */
  return 0.5*(xl+xu);

#undef FKT
} /* end of _unur_itdr_find_xt() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_itdr_debug_init( const struct unur_gen *gen, int error )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   error ... error code of hat generation                             */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = itdr (inverse transformed density rejection)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_itdr_sample",gen->genid);
  if (gen->variant & ITDR_VARFLAG_VERIFY) fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(log,"%s: bx = %g",gen->genid, GEN->bx);
  fprintf(log,"%s\n", (gen->set & ITDR_SET_BX) ? "" : " [computed]");
  fprintf(log,"%s: cp = %g",gen->genid, GEN->cp);
  fprintf(log,"%s\n", (gen->set & ITDR_SET_CP) ? "" : " [computed]");
  fprintf(log,"%s: ct = %g",gen->genid, GEN->ct);
  fprintf(log,"%s\n", (gen->set & ITDR_SET_CT) ? "" : " [computed]");

  fprintf(log,"%s: xp = %g\n",gen->genid, GEN->xp);
  fprintf(log,"%s: xt = %g\n",gen->genid, GEN->xt);

  fprintf(log,"%s: alphap = %g, betap = %g\n",gen->genid, GEN->alphap, GEN->betap);
  fprintf(log,"%s: Tfxt = %g, dTfxt = %g\n",gen->genid, GEN->Tfxt, GEN->dTfxt);

  fprintf(log,"%s: by = %g\n",gen->genid, GEN->by);

  fprintf(log,"%s: Area = %g + %g + %g = %g\n",gen->genid,
	  GEN->Ap, GEN->Ac, GEN->At, GEN->Atot);

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: **** INIT %s ***\n",gen->genid,
	  (error==UNUR_SUCCESS) ? "successful" : "failed" );   
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_itdr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
