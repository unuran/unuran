/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    numerical inversion of cumulative distribution function      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the c.d.f.                                                *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      c.d.f. at mode                                                       *
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
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define NINV_VARFLAG_NEWTON   0x1u   /* use Newton's method                  */
#define NINV_VARFLAG_REGULA   0x2u   /* use regula falsi (default)           */

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

#define NINV_DEBUG_SAMPLE        0x100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NINV_SET_MAX_ITER     0x001u   /* number of maximal interations      */
#define NINV_SET_X_RESOLUTION 0x002u   /* maximal relative error in x        */
#define NINV_SET_START        0x004u   /* intervals at start (left/right)    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_ninv_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample_regula( struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling (regula falsi).                                            */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.ninv        /* data for parameter object         */
#define GEN       gen->data.ninv        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define CDF(x) ((*(DISTR.cdf))((x),DISTR.params,DISTR.n_params))    /* call to p.d.f. */

/*---------------------------------------------------------------------------*/
/* macros                                                                    */

/* sign of x                                                                 */
#define sgn(x) (((x)<0) ? -1 : 1 )


/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_ninv_new( struct unur_distr *distr )
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
  if (distr==NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return NULL; }

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"c.d.f."); return NULL; }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_NINV_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR.max_iter    = 40;         /* maximal number of iterations              */
  PAR.rel_x_resolution = 1.0e-16;   /* maximal relative error in x           */
  PAR.sl           = -10.;      /* left boundary of interval                 */
  PAR.sr           = 10.;       /* right boundary of interval                */

  par->method      = UNUR_METH_NINV;  /* method and default variant          */
  par->variant     = NINV_VARFLAG_REGULA;  /* default variant                */
  par->set         = 0u;              /* inidicate default parameters        */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  par->genid       = _unur_set_genid(GENTYPE);/* set generator id            */
  par->debug       = UNUR_DEBUGFLAG_DEFAULT;  /* set default debugging flags */

  /* routine for starting generator */
  par->init = unur_ninv_init;

  return par;

} /* end of unur_ninv_new() */

/*****************************************************************************/

int unur_ninv_use_newton( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use Newton's method                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if (!par) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return 0; }

  /* check input */
  _unur_check_par_object( NINV );

  /* check new parameter for generator */
  if (! par->DISTR_IN.pdf) {
    _unur_error(par->genid,UNUR_ERR_DISTR_REQUIRED,"p.d.f."); return 0; }

  /* store date */
  par->variant = NINV_VARFLAG_NEWTON;

  return 1;

} /* end of unur_ninv_use_newton() */

/*---------------------------------------------------------------------------*/

int unur_ninv_use_regula( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use regula falsi                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if (!par) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return 0; }

  /* check input */
  _unur_check_par_object( NINV );

  /* store date */
  par->variant = NINV_VARFLAG_REGULA;

  return 1;

} /* end of unur_ninv_use_regula() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_max_iter( struct unur_par *par, int max_iter )
     /*----------------------------------------------------------------------*/
     /* set number of maximal iterations                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if (!par) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return 0; }

  /* check input */
  _unur_check_par_object( NINV );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(par->genid,UNUR_ERR_PAR_SET,"maximal iterations");
    return 0;
  }

  /* store date */
  PAR.max_iter = max_iter;

  /* changelog */
  par->set |= NINV_SET_MAX_ITER;

  return 1;

} /* end of unur_ninv_set_max_iter() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution)
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   x_resolution ... maximal relative error in x                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if (!par) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return 0; }

  /* check input */
  _unur_check_par_object( NINV );

  /* check new parameter for generator */
  if (x_resolution < DBL_EPSILON) {
    _unur_warning(par->genid,UNUR_ERR_PAR_SET,"x resolution");
    return 0;
  }

  /* store date */
  PAR.rel_x_resolution = x_resolution;

  /* changelog */
  par->set |= NINV_SET_X_RESOLUTION;

  return 1;

} /* end of unur_ninv_set_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_start( struct unur_par *par, double sl, double sr )
     /*----------------------------------------------------------------------*/
     /* set intervals at start (left/right)                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   sl    ... left boundary of interfval                               */
     /*   sr    ... right boundary of interfval                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if (!par) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,""); return 0; }

  /* check input */
  _unur_check_par_object( NINV );

  /* check new parameter for generator */
  if (sl >= sr) {
    _unur_warning(par->genid,UNUR_ERR_PAR_SET,"sl >= sr");
    return 0;
  }

  /* store date */
  PAR.sl = sl;
  PAR.sr = sr;

  /* changelog */
  par->set |= NINV_SET_START;

  return 1;

} /* end of unur_ninv_set_start() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

struct unur_gen *
unur_ninv_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
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
  if (par==NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,"");
    return NULL;
  }

  /* check input */
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(par->genid,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_ninv_create(par);
  if (!gen) { free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_ninv_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_ninv_init() */

/*****************************************************************************/

double
unur_ninv_sample_regula( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use regula falsi)                             */
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
  double x1, x2, a, xtmp;  /* points for RF */
  double f1, f2, ftmp;     /* function values at x1, x2, xtmp */
  double length;           /* (gerichtete) laenge des Intervalls mit ZW */ 
  double lengthrel;        /* relative length of interval */
  int  lengthsgn;    /* "richtung" des Intervalls */
  double step;              /* Vergr"o"sert Startinvervall bis ZW gefunden */
  double dx;               /* RF-Schrittgr"o"se */
  int count = 0;     /* Z"ahler f"ur "keine ZW" */
  int count2 = 0;    /* Z"ahler f"ur "kein Fortschritt" */
  int i;                   /* Schleifenzahler */
    
  double u;     /* uniform random number */

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_NINV_GEN,0.);

  /* initialize starting interval */
  x1 =  GEN.sl;      /* left boudary of interval */
  x2 =  GEN.sr;      /* right boudary of interval */
  if (x1>=x2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); return 0.; }

  /* sample from uniform pseudorandom number generator */
  u = _unur_call_urng(gen);

  /* compute c.d.f. at interval boundaries */
  f1 = CDF(x1) - u;
  f2 = CDF(x2) - u;

  /* search for interval with changing signs */
  step = 1.;  /* Startintervall zu klein -> um 2^n * gap vergr"o"sern */ 
  while ( f1*f2 > 0. ) {
    if ( f1 > 0. ) {/* untere Grenze zu gross */    
      x1 -= step;   
      f1 = CDF(x1) - u;
    }
    else {         /* obere Grenze zu klein */
      x2 += step;
      f2 = CDF(x2) - u;
    }
    /* increase step width */
    step *= 2.; 
  }  
  a = x1;

  /* Sekantensuche, ZW wird beibehalten */
  for (i=0; i < GEN.max_iter; i++) { /* noch zu andernde schleife */
    count++;

    /* f2 immer kleiner (besser), notfalls Tausch */
    if ( fabs(f1) < fabs(f2) ) {   /* Dreieckstausch */
      xtmp = x1; ftmp = f1;
      x1 = x2;   f1 = f2;
      x2 = xtmp; f2 = ftmp;
    }
    
    if ( f2 == 0. ) 
      /* genauer Treffer -- sehr unwahrscheinlich */
      break; /* -> finished */
    
    if ( f1*f2 <= 0) {  /* ZeichenWechsel vorhanden */
      count = 0;   /* zaehler fuer bisektion wird rueckgestellt */
      a = x1;      /* [a, x2] enthaelt ZW */
    }
    
    length = x2 - a;  /* gerichtete laenge */
    //    lengthrel = fabs(length) / ((fabs(x1)<fabs(x2)) ? fabs(x2) : fabs(x1));
    lengthrel = fabs(length);
    lengthsgn = sgn(length);
 
    if ( lengthrel < GEN.rel_x_resolution || count2 > 1 )
      /* x-genauigkeit erreicht -> finished */
      break; /* -> finished */

   /* Sekanten-Schritt */
    dx = ( f1-f2==0. ) ? 0. : f2*(x2-x1)/(f2-f1) ;  

    if (fabs(dx) <= GEN.rel_x_resolution) {
      count2++;       /* 2mal hintereinander kein fortschritt -> abbruch */
      dx = length/2.; /* Bisektionsschritt */
    }    
    else   /* Abbruch-Z"ahler r"uckgesetzt */
      count2 = 0;  
    
    /* kein  ZW  || Schritt fuhrt aus Intervall  */
    if ( count > 1 || (lengthrel-GEN.rel_x_resolution) <= dx*lengthsgn  )  /** TODO **/
      dx = length/2.; /* Bisektionsschritt */
  
    /* Update der Punkte */    
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
 
  }  /* for-schleife ende */

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (gen->debug & NINV_DEBUG_SAMPLE)
      _unur_ninv_debug_sample_regula( gen,u,x2,f2,i );
#endif

  return x2;

} /* end of unur_ninv_sample_regula() */

/*---------------------------------------------------------------------------*/

double
unur_ninv_sample_newton( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use Newton's method)                          */
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
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_NINV_GEN,0.);

  return 0.;

} /* end of unur_ninv_sample_newton() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

void
unur_ninv_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of unur_ninv_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NINV_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* copy generator identifier */
  gen->genid = par->genid;

  /* routines for sampling and destroying generator */
  switch (par->variant) {
  case NINV_VARFLAG_NEWTON:
    SAMPLE = unur_ninv_sample_newton;
    break;
  case NINV_VARFLAG_REGULA: default:
    SAMPLE = unur_ninv_sample_regula;
    break;
  }

  gen->destroy = unur_ninv_free;

  /* copy parameters into generator object */
  GEN.max_iter = PAR.max_iter;  /* maximal number of iterations              */
  GEN.rel_x_resolution = PAR.rel_x_resolution; /* maximal relative error in x*/
  GEN.sl = PAR.sl;         /* left boundary of interval boundaries at start  */
  GEN.sr = PAR.sr;         /* right boundary of interval boundaries at start */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_ninv_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_ninv_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = ninv (numerical inversion of c.d.f.)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = unur_ninv_sample",gen->genid);
  switch (par->variant) {
  case NINV_VARFLAG_NEWTON:
    fprintf(log,"_newton\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    fprintf(log,"_regula\n");
    break;
  }

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_init() */

/*---------------------------------------------------------------------------*/

static void
_unur_ninv_debug_sample_regula( struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (regula falsi)                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  log = unur_get_stream();

  fprintf(log,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d interations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN.max_iter);

} /* end of _unur_ninv_debug_sample_regula() */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
