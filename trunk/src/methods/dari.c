/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dari.c                                                       *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    discrete automatic rejection inversion                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PMF of a T(-1/2)-concave distribution;                         *
 *      produce a value x consistent with its PMF.                           *
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
 *   [1] Hoermann W. and G.Derflinger (1997):                                *
 *       An automatic generator for a large class of discrete unimodal       *
 *       distributions, in A.R. Kaylan and A. Lehmann, ESM 97, pp 139-144    *
 *                                                                           *
 *   [2] Hoermann W. and G.Derflinger (1996):                                *
 *       Rejection-inversion to generate variates from monotone discrete     *
 *       distributions, ACM TOMACS 6(3), 169-184                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dari.h"

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DARI_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DARI_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DARI_SET_CFACTOR        0x001u
#define DARI_SET_TABLESIZE      0x002u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DARI"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dari_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute hat.                                                              */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dari_sample( struct unur_gen *generator );
static int _unur_dari_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_dari_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dari_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       par->data.dari        /* data for parameter object         */
#define GEN       gen->data.dari        /* data for generator object         */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */     

#define PMF(x)    _unur_discr_PMF((x),(gen->distr))    /* call to PMF        */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dari_new( const struct unur_distr *distr )
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
     /*                                                                      */
     /* comment:                                                             */
     /*   if the sum over the PMF is not close to 1 it is necessary to       */
     /*   set pmf_sum to an approximate value of its sum (+/- 30 % is ok).  */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;
  
  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );
  
  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  
  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_DARI_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR.c_factor  = 0.664;
  /* optimal value for the normal distribution, which is good for all        */
  /* bell-shaped densities. The minimax approach for that transformation     */
  /* has c_factor = 2.                                                       */

  PAR.squeeze   = 0; /* no squeezes by default as squeezes slow down the     */
  /* sampling for most distributions if PAR.size is big enough. Squeeze is   */
  /* important for the speed only when small samples are required or when    */
  /* the domain of the distribution is very big. (much bigger than PAR.size) */

  PAR.size      = 100; /*size of table that stores the "rejection point" for */
  /* all integers close to the mode when needed the first time while         */
  /* sampling; can speed up the generation considerably.                     */

  par->method   = UNUR_METH_DARI;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dari_init;

  return par;

} /* end of unur_dari_new() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_dari_set_cpfactor( struct unur_par *par, double cpfactor )
     /*----------------------------------------------------------------------*/
     /* set factor for position of left and right construction point         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   cfactor ... factor                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* check new parameter for generator */
  /** TODO: welche werte fuer c sind zulaessig / sinnvoll ? 
  zulaessig ist jedes c>0, man koennte damit falsche Flaechenangaben kompensieren.
  Wenn sum genau bekannt ist, ist ein c > 2 (2 ist der minimax approach) so weit
  ich weiss nie sinnvoll. Ich denke aber, das sollte man besser nicht prinzipiell
  verbieten, hoechstens eine warnung.**/
  if (cpfactor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor <= 0");
    return 0;
  }

  if (cpfactor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp-factor > 2 not recommended. skip");

  /* store date */
  PAR.c_factor = cpfactor;

  /* changelog */
  par->set |= DARI_SET_CFACTOR;

  return 1;

} /* end of unur_dari_set_cpfactor() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_squeeze( struct unur_par *par, int squeeze )
     /*----------------------------------------------------------------------*/
     /* turn on/off using squeezes                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   squeeze ... 0 = no squeeze,  !0 = use squeeze                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );
  
  /* check input */
  _unur_check_par_object( par,DARI );

  /* store data */
  PAR.squeeze = squeeze;

  /* o.k. */
  return 1;

} /* end of unur_dari_set_squeeze() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_tablesize( struct unur_par *par, int size )
     /*----------------------------------------------------------------------*/
     /* set size of table                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   size   ... table size                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* check parameter */
  if (size < 0) {  
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid table size");
    return 0;
  }
  
  /* store data */
  PAR.size = size;

  /* changelog */
  par->set |= DARI_SET_TABLESIZE;

  /* o.k. */
  return 1;
} /* end of unur_dari_set_tablesize() */
  
/*---------------------------------------------------------------------------*/

int
unur_dari_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | DARI_VARFLAG_VERIFY) : (par->variant & (~DARI_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_dari_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= DARI_VARFLAG_VERIFY;
    SAMPLE = _unur_dari_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~DARI_VARFLAG_VERIFY;
    SAMPLE = _unur_dari_sample;
  }

  /* o.k. */
  return 1;

} /* end of unur_dari_chg_verify() */

/*****************************************************************************/

int
unur_dari_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );

  /* set new parameters in distribution object */
  return unur_distr_discr_set_pmfparams(gen->distr,params,n_params);

} /* end of unur_dari_chg_pmfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_dari_chg_domain( struct unur_gen *gen, int left, int right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity.           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_dari_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_mode( struct unur_gen *gen, int mode )
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );
  
  /* copy parameters */
  DISTR.mode = mode;
  
  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return 1;
} /* end of unur_dari_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );

  return unur_distr_discr_upd_mode( gen->distr );
} /* end of unur_dari_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_pmfsum( struct unur_gen *gen, double sum )
     /*----------------------------------------------------------------------*/
     /* change sum over PMF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   sum   ... sum                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );
  
  /* check new parameter for generator */
  if (sum <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"sum <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.sum = sum;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_dari_chg_pmfsum() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_pmfsum( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute sum over PMF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );

  return unur_distr_discr_upd_pmfsum( gen->distr );
} /* end of unur_dari_upd_pmfsum() */

/*****************************************************************************/

struct unur_gen *
_unur_dari_init( struct unur_par *par )
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

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_DARI ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DARI_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dari_create(par);
  if (!gen) return NULL;

  /* create hat and squeeze (setup procedure) */
  if ( !_unur_dari_hat(gen) ) {
    /* error */
    _unur_dari_free(gen);
    return NULL;
  }

  /* hat successfully created */
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dari_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_dari_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_dari_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DARI_PAR,NULL);
  
  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DARI_GEN);

  /* copy distribution object into generator object */
  gen->distr = _unur_distr_clone( par->distr );

  /* check for required data: mode */
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (!unur_distr_discr_upd_mode( gen->distr )) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      _unur_distr_free(gen->distr); free(gen);
      return NULL; 
    }
  }

  /* check for required data: sum over PMF */
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (!unur_distr_discr_upd_pmfsum( gen->distr ))
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF; use default");

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & DARI_VARFLAG_VERIFY) ? _unur_dari_sample_check : _unur_dari_sample;
  gen->destroy = _unur_dari_free;
  gen->clone = _unur_dari_clone;

  /* copy some parameters into generator object */
  GEN.squeeze = PAR.squeeze;        /* squeeze yes/no?                       */
  GEN.c_factor = PAR.c_factor;      /* constant for choice of design point   */

  /* size of auxiliary table; 0 for none
     it cannot be larger than the given domai (avoid overflow) */
  if ((unsigned)DISTR.BD_RIGHT - (unsigned)DISTR.BD_LEFT < INT_MAX)
    GEN.size = min(PAR.size,DISTR.BD_RIGHT-DISTR.BD_LEFT+1);
  else /* length of interval > INT_MAX */
    GEN.size = PAR.size;

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */

  /* allocate */
  GEN.hp = (GEN.size > 0) ? _unur_malloc( GEN.size * sizeof(double) ) : NULL;
  GEN.hb = (GEN.size > 0) ? _unur_malloc( GEN.size * sizeof(char) )   : NULL;

  /* initialize parameters */
  /** TODO: diese initialisierung ist nur zur Sicherheit,
      man sollte kontrollieren, ob das so wirklich noetig ist **/
  GEN.vt=0.;            /* total volume below hat                            */
  GEN.vc=0.;            /* volume below center part                          */
  GEN.vcr=0.;           /* volume center and right together                  */

  GEN.xsq[0]=0.;        /* value necesary for the squeeze computation        */
  GEN.xsq[1]=0.;        /* value necesary for the squeeze computation        */
  GEN.y[0]=0.;          /* value of the transformed density in points of contact */
  GEN.y[1]=0.;          /* value of the transformed density in points of contact */
  GEN.ys[0]=0.;         /* the slope of the transformed hat                  */
  GEN.ys[1]=0.;         /* the slope of the transformed hat                  */
  GEN.ac[0]=0.;         /* left and right starting point of the uniform hat
                           in the center                                     */
  GEN.ac[1]=0.;         /* left and right starting point of the uniform hat
                           in the center                                     */

  GEN.pm=0.;            /* mode probability                                  */
  GEN.Hat[0]=0.;        /* point where the hat starts for the left and
			    the right tail                                   */
  GEN.Hat[1]=0.;        /* point where the hat starts for the left and
			    the right tail                                   */

  GEN.m=0;              /* mode                                              */
  GEN.x[0]=0;           /* points of contact left and right of the mode      */
  GEN.x[1]=0;           /* points of contact left and right of the mode      */
  GEN.s[0]=0;           /* first and last integer of the center part         */
  GEN.s[1]=0;           /* first and last integer of the center part         */
  GEN.n[0]=0;           /* contains the first and the last i 
			   for which values are stored in table              */
  GEN.n[1]=0;           /* contains the first and the last i 
			   for which values are stored in table              */


  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_dari_create() */

/*---------------------------------------------------------------------------*/
#define T(x) (-1./sqrt(x))
#define F(x) (-1./(x))
#define FM(x) (-1./(x))
#define N0 (GEN.n[0])

int
_unur_dari_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute hat                                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  int sign[2] = {-1,1};
  int b[2], d, i, j;
  double v[2], at[2];
  double t0 = 1.;
  char setup = 1;
  char rep = 1;
  
  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_DARI_GEN, 0 );
  
  /* check mode. we assume unimodality 
     since otherwise the PMF would not be T-concave! */
  if (DISTR.BD_LEFT > DISTR.mode)
    DISTR.mode = DISTR.BD_LEFT;
  else if (DISTR.BD_RIGHT < DISTR.mode)
    DISTR.mode = DISTR.BD_RIGHT;

  /* sum must not be zero */
  if (DISTR.sum <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"sum <= 0");
    return 0;
  }
 
  /* Step 0: setup */
  GEN.m = DISTR.mode;
  b[0] = DISTR.BD_LEFT;
  b[1] = DISTR.BD_RIGHT;
  GEN.pm = PMF(GEN.m);
  d = max(2, (int)( GEN.c_factor/(GEN.pm/DISTR.sum)));

  /* step 0.1 */
  do {
    for(i=0; i<=1; i++) {
      GEN.x[i] = GEN.m + sign[i] * d;
      if (sign[i]*GEN.x[i]+1 > sign[i]*b[i]) {
	v[i] = 0; 
	GEN.s[i] = b[i];
      }
      else {
	GEN.y[i] = T( PMF(GEN.x[i]) );
	GEN.ys[i] = sign[i] * (T( PMF(GEN.x[i]+sign[i])) - GEN.y[i]);
	if (GEN.ys[i]*sign[i] > -DBL_EPSILON) {
	  setup = -setup; /* indicate that the hat is not ok */
	  i = 1; 
	}
        else {
	  GEN.s[i] = (int)(0.5+GEN.x[i]+(T(GEN.pm)-GEN.y[i])/GEN.ys[i]);
	  GEN.Hat[i] = ( F(GEN.y[i]+GEN.ys[i]*(GEN.s[i]+sign[i]*1.5-GEN.x[i])) /
			 GEN.ys[i]-sign[i]*PMF(GEN.s[i]+sign[i]) ); 
	  at[i] = GEN.x[i] + (FM(GEN.ys[i]*GEN.Hat[i])-GEN.y[i]) / GEN.ys[i]; 
          if(GEN.squeeze)
	    GEN.xsq[i] = sign[i]*(at[i]-(GEN.s[i]+sign[i]));
	  v[i] = sign[i]*(F(GEN.y[i]+GEN.ys[i]*(b[i]+sign[i]*0.5-GEN.x[i]))/
			  GEN.ys[i]-F(GEN.y[i]+GEN.ys[i]*(at[i]-GEN.x[i]))/GEN.ys[i]);
	}
      }
      if (setup>0)
	GEN.ac[i] = GEN.s[i] + sign[i]*(PMF(GEN.s[i])/GEN.pm-0.5);
    }

    /* step 0.2 */
    if(setup>0) {
      GEN.vc = GEN.pm*(GEN.ac[1]-GEN.ac[0]); 
      GEN.vt = GEN.vc+v[0]+v[1];
      GEN.vcr = GEN.vc+v[1];

      /* step 0.3 */ 
      GEN.n[0] = max(b[0],GEN.m - GEN.size/2);
      GEN.n[1] = GEN.n[0] + GEN.size - 1;
      if (GEN.n[1] > b[1]) {
	GEN.n[1] = b[1];
	GEN.n[0] = GEN.n[1]- GEN.size + 1;
      }
      /* initialize table */
      for (j=0; j<GEN.size; j++)
	GEN.hb[j] = 0;
    }

    /* setup == 1 first try, up to now ok,  ==2 second try, up to now ok */
    /* setup == -1 first try, not ok,  == -2 second try, not ok          */

    if (setup == 1 || setup == -1) {
      t0= 2. * DISTR.sum;
      if (setup==1 && GEN.vt<=t0)
	rep=0;
      else { 
#ifdef UNUR_ENABLE_LOGGING
	/* write info into log file */
	if (gen->debug) _unur_dari_debug_init(gen);
#endif
	setup = 2;
	d = ((int)t0) / GEN.pm;
      }
    }
    else 
      rep=0; 
  } while(rep);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dari_debug_init(gen);
#endif

  if (setup == -2 || GEN.vt > 100.*t0 || GEN.vt <0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"; Area below hat too large or zero!! possible reasons: PDF, mode or area below PMF wrong;  or PMF not T-concave");
    return 0;
  }

  /* o.k. */
  return 1;

} /* end of _unur_dari_hat() */

/*---------------------------------------------------------------------------*/

int
unur_dari_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,gen,0 );
  _unur_check_gen_object( gen,DARI );
  
  /* compute hat  */
  return _unur_dari_hat( gen );
} /* end of unur_dari_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dari_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.dari

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DARI_GEN,NULL);

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  clone->distr = _unur_distr_clone( gen->distr );

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = _unur_gen_clone( gen->gen_aux );

  /* copy additional data */
  if (GEN.size > 0) {
    CLONE.hp = _unur_malloc( GEN.size * sizeof(double) );
    memcpy( CLONE.hp, GEN.hp, GEN.size * sizeof(double) );
    CLONE.hb = _unur_malloc( GEN.size * sizeof(char) );
    memcpy( CLONE.hb, GEN.hb, GEN.size * sizeof(char) );
  }
  
  return clone;

#undef CLONE
} /* end of _unur_dari_clone() */

/*****************************************************************************/

int
_unur_dari_sample( struct unur_gen *gen )
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
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  double U, h;
  double X = 0.;
  int k,i;
  int sign[2] = {-1,1};

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DARI_GEN,0.);

  /* step 1.0 */
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN.vt;

    /* step 1.1 */
    if (U<=GEN.vc) {
      X = U * (GEN.ac[1]-GEN.ac[0]) / GEN.vc + GEN.ac[0]; 
      k = (int)(X+0.5);
      if (k<GEN.m)
	i=0; 
      else
	i=1;
      if (GEN.squeeze && sign[i]*(GEN.ac[i]-GEN.s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN.n[i]) {
	if (!GEN.hb[k-N0]) {
	  GEN.hp[k-N0] = 0.5 - PMF(k)/GEN.pm;
	  GEN.hb[k-N0] = 1;
	}
	h = GEN.hp[k-N0];
      }
      else {
	h = 0.5-PMF(k)/GEN.pm;
      }
      if (h <= sign[i]*(k-X))
	return k;
    }

    /*step 1.2*/ 
    else {
      if (U<= GEN.vcr) {
	i = 1;
	U -= GEN.vc;
      } 
      else {
	i = 0;
	U -= GEN.vcr;
      }

      U = GEN.Hat[i] + sign[i]*U; 
      X = GEN.x[i] + (FM(U*GEN.ys[i])-GEN.y[i]) / GEN.ys[i];
      k = (int)(X+0.5);

      if (GEN.squeeze && (sign[i]*k <= sign[i]*GEN.x[i]+1) && (GEN.xsq[i] <= sign[i]*(X-k))) 
	return k;

      if (sign[i]*k <= sign[i]*GEN.n[i]) {
	if (!GEN.hb[k-N0]) {
	  GEN.hp[k-N0] = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i] - PMF(k);
	  GEN.hb[k-N0] = 1;
	}
	h = GEN.hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i]-PMF(k);
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
} /* end of _unur_dari_sample() */

/*****************************************************************************/

int
_unur_dari_sample_check( struct unur_gen *gen )
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
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  double U, h;
  double X = 0.;
  double hkm05;
  int k,i;
  int sign[2] = {-1,1};
  
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DARI_GEN,0.);
  
  /* step 1.0 */
  while (1) {
    U = _unur_call_urng(gen->urng) * GEN.vt;

    /* step 1.1 */
    if (U <= GEN.vc) {
      X = U * (GEN.ac[1]-GEN.ac[0]) / GEN.vc + GEN.ac[0]; 
      k = (int)(X+0.5);
      if (k<GEN.m)
	i=0; 
      else
	i=1;
      if (GEN.squeeze && sign[i]*(GEN.ac[i]-GEN.s[i]) > sign[i]*(X-k))
	return k;
      if (sign[i]*k <= sign[i]*GEN.n[i]) {
	if (!GEN.hb[k-N0]) {
	  GEN.hp[k-N0] = 0.5 - PMF(k)/GEN.pm;
	  GEN.hb[k-N0] = 1;
	}
	h = GEN.hp[k-N0];
	/* CHECKING HAT */
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_stream_printf(gen->genid,__FILE__,__LINE__,
			      "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN.pm ); 
        }
	/* end CHECKING HAT */
      }
      else {
	h = 0.5 - PMF(k)/GEN.pm;
	/* CHECKING HAT */
	/* here UNUR_EPSILON can be too small for distributions that 
	   have two neighbouring hats. */
	if (h+UNUR_EPSILON*100.<-0.5) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PMF(i) > hat(i) for centerpart");
	  _unur_stream_printf(gen->genid,__FILE__,__LINE__,
			      "i %d PMF(x) %.20e hat(x) %.20e", k,PMF(k),GEN.pm ); 
        }
	/* end CHECKING HAT */
      }
      if (h <= sign[i]*(k-X))
	return k;
    }

    /*step 1.2*/ 
    else {
      if (U<= GEN.vcr) {
	i = 1;
	U -= GEN.vc;
      } 
      else {
	i = 0;
	U -= GEN.vcr;
      }

      U = GEN.Hat[i] + sign[i]*U; 
      X = GEN.x[i] + (FM(U*GEN.ys[i])-GEN.y[i]) / GEN.ys[i];
      k = (int)(X+0.5);
      /* this is for a very rare case that for k of the tail closest to 
	 the mode an x value farer away than 0.5 is geenrated. */
      if(k==GEN.s[i]) 
	k += sign[i];

      if (GEN.squeeze && (sign[i]*k <= sign[i]*GEN.x[i]+1) && (GEN.xsq[i] <= sign[i]*(X-k))) 
	return k;

      if (sign[i]*k <= sign[i]*GEN.n[i]) {
	if(!GEN.hb[k-N0]) {
	  GEN.hp[k-N0] = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i] - PMF(k); 

	  /* CHECKING HAT: (only necessary if(k!=GEN.s+1) as for the border
                            the hat is by construction correct)
	     tests if Hat too low i.e.: (H(k+0.5)-  p_k < H(k-0.5)) */
          if(k != GEN.s[i]+sign[i]) {
            hkm05 = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k-sign[i]*0.5-GEN.x[i])) / GEN.ys[i];
	    if (GEN.hp[k-N0]+UNUR_EPSILON < hkm05) {
	      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			  "for tailpart hat too low, ie hp[k] < H(k-0.5)");
	      _unur_stream_printf(gen->genid,__FILE__,__LINE__,
				  "k %d hp  %.20e H(k-0.5) %.20e ", k,GEN.hp[k-N0],hkm05 ); 
            }
	  }
	  GEN.hb[k-N0] = 1;
	}
	h = GEN.hp[k-N0];
      }
      else {
	h = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i] - PMF(k);
	/* CHECKING HAT:(only necessary if(k!=GEN.s+1) as for the border
                            the hat is by construction correct)
	   tests if Hat too low i.e.: (H(k+0.5)-p_k < H(k-1/2)) */
        hkm05 = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k-sign[i]*0.5-GEN.x[i])) / GEN.ys[i];
        if(k != GEN.s[i]+sign[i]) {
	  if (h+UNUR_EPSILON < hkm05) {
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
			"PMF(i) > hat(i) for tailpart");
	    _unur_stream_printf(gen->genid,__FILE__,__LINE__,
				"k %d h  %.20e H(k-0.5) %.20e ", k,h,hkm05 ); 
	  }
        }
      }
      if (sign[i]*U >= h)
	return k;
    }
  }
  
  
  
#if 0
  if (sign[i]*k <= sign[i]*GEN.n[i]) {
    if(!GEN.hb[k-N0]) {
      GEN.hp[k-N0] = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i] - PMF(k); 
      
      /* CHECKING HAT:
	 tests if Hat too low i.e.: (H(k+0.5)-p_k < H(k-1/2)) */
      if (GEN.hp[k-N0]+UNUR_EPSILON 
	  < sign[i] * F(GEN.y[i]+GEN.ys[i]*(k-sign[i]*0.5-GEN.x[i])) / GEN.ys[i]) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		    "PMF(i) > hat(i) for tailpart");
	_unur_stream_printf(gen->genid,__FILE__,__LINE__,
			    "i %d PMF(x) %e ", k,PMF(k) ); 
      }
      GEN.hb[k-N0] = 1;
    }
    h = GEN.hp[k-N0];
  }
  else {
    h = sign[i] * F(GEN.y[i]+GEN.ys[i]*(k+sign[i]*0.5-GEN.x[i])) / GEN.ys[i] - PMF(k);
    /* CHECKING HAT:
       tests if Hat too low i.e.: (H(k+0.5)-p_k < H(k-1/2)) */
    if (h+UNUR_EPSILON
	< sign[i] * F(GEN.y[i]+GEN.ys[i]*(k-sign[i]*0.5-GEN.x[i]))/GEN.ys[i]) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "PMF(i) > hat(i) for tailpart");
      _unur_stream_printf(gen->genid,__FILE__,__LINE__,
			  "i %d PMF(x) %e ", k,PMF(k) );
    }
  }
  if (sign[i]*U >= h)
    return k;
  
#endif

} /* end of _unur_dari_sample_check() */

#undef N0
#undef T
#undef F
#undef FM

/*****************************************************************************/

void
_unur_dari_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DARI ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);
  
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */
  
  /* free two auxiliary tables */
  if (GEN.hp)   free(GEN.hp);
  if (GEN.hb)   free(GEN.hb);
  
  /* free memory */
  _unur_distr_free(gen->distr);
  _unur_free_genid(gen);

  COOKIE_CLEAR(gen);
  free(gen);

} /* end of _unur_dari_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dari_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DARI_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = dari (discrete automatic rejection inversion)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );

  fprintf(log,"%s: sampling routine = _unur_dari_sample",gen->genid);
  if (gen->variant & DARI_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s:area below hat: total %f center: %f, left tail %f, right tail %f\n", gen->genid,
	  GEN.vt, GEN.vc, GEN.vt-GEN.vcr, GEN.vcr-GEN.vc);
  fprintf(log,"%s: mode %d and mode probability %f\n",gen->genid, GEN.m, GEN.pm); 
  for(i=0;i<=1;i++) {
    fprintf(log,"%s:i=%d: x=%d; Hat=%f; ac=%f; s=%d;\n", gen->genid,
	    i, GEN.x[i], GEN.Hat[i], GEN.ac[i], GEN.s[i]);
    fprintf(log,"%s:i=%d: xsq=%f; y=%f; ys=%f; n:=%d (for aux.table)\n", gen->genid,
	    i, GEN.xsq[i], GEN.y[i], GEN.ys[i], GEN.n[i]);
  }
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_dari_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
