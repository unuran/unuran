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
/* Starting interval includes this percentage of all univariate rand numbers */
/* must be > 0. and < 1.                                                     */
#define INTERVAL_COVERS  (.9)

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

static double _unur_ninv_regula( struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* algorithm: regula falsi                                                   */
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
#define PDF(x) ((*(DISTR.pdf))((x),DISTR.params,DISTR.n_params))    /* call to p.d.f. */


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
  _unur_check_NULL( GENTYPE,distr,NULL );

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
  par->distr       = distr;        /* pointer to distribution object         */

  /* set default values */
  PAR.max_iter  = 40;              /* maximal number of iterations           */
  PAR.rel_x_resolution = 1.0e-8 ;  /* maximal relative error in x            */
  PAR.s[0]      = 0.;              /* left boundary of interval              */
  PAR.s[1]      = 0.;              /* right boundary of interval             */
  /* if the PAR.s[i] are identical, the PAR.s[i] are set in unur_ninv_init   */
  /* such that INTERVAL_COVERS *100 Percent of the used univariate random    */
  /* numbers will be in that interval                                        */

  par->method   = UNUR_METH_NINV;          /* method and default variant     */
  par->variant  = NINV_VARFLAG_REGULA;     /* default variant                */
  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */

  par->genid    = _unur_set_genid(GENTYPE);/* set generator id               */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

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
  _unur_check_NULL( GENTYPE,par,0 );

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
  _unur_check_NULL( GENTYPE,par,0 );

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
  _unur_check_NULL( GENTYPE,par,0 );

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
  _unur_check_NULL( GENTYPE,par,0 );

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

int
unur_ninv_set_start( struct unur_par *par, double s1, double s2, double s3 )
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /*   Muller/Brent:  s1. s2, s3   starting points                        */
     /* arguments that are used by method are ignored.                       */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*                                                                      */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( NINV );

  /* store date */
  PAR.s[0] = s1;
  PAR.s[1] = s2;
  PAR.s[2] = s3;

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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(par->genid,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_ninv_create(par);
  if (!gen) { free(par); return NULL; }


  /* check arguments */
  switch (par->variant) {

  case NINV_VARFLAG_REGULA:
    if (GEN.s[0] > GEN.s[1]) {
      /* swap interval boundaries */
      double tmp;
      tmp = GEN.s[0]; GEN.s[0] = GEN.s[1]; GEN.s[1] = tmp;
    }

    if (GEN.s[0] == GEN.s[1]) {
      /* length of interval == 0 -> choose bounderies with            */
      /*  INTERVAL_COVERS *100 % chance for sign change in interval   */
      GEN.s[0] = -10.;      /* arbitrary starting value               */
      GEN.s[1] =  10.;      /* arbitrary starting value               */
      GEN.s[0] = _unur_ninv_regula(gen, (1.-INTERVAL_COVERS)/2. );
      GEN.s[1] = GEN.s[0] + 10.;   /* arbitrary interval length       */
      GEN.s[1] = _unur_ninv_regula(gen, (1.+INTERVAL_COVERS)/2. );
    }
    break;

  case NINV_VARFLAG_NEWTON:
    /* nothing to do */
    break;
  }


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
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_regula(gen, _unur_call_urng(gen) ) ;
}


static double 
_unur_ninv_regula( struct unur_gen *gen, double u )
     /*---------------------------------------------------------------------*/
     /*   algorithm: regula falsi                                           */
     /*                                                                     */
     /*   parameters:                                                       */
     /*      gen ... pointer to generator object                            */
     /*      u   ... random number (uniform distribution)                   */
     /*   return:                                                           */
     /*     double (sample from random variate)                             */
     /*                                                                     */
     /*   error:                                                            */
     /*     return 0.                                                       */
     /*---------------------------------------------------------------------*/
{ 
    
  double x1, x2, a, xtmp;/* points for RF                                   */
  double x2abs;          /* absolute value of x2                            */
  double f1, f2, ftmp;   /* function values at x1, x2, xtmp                 */
  double length;         /* oriented length of the interval with sign change*/
  double lengthabs;      /* absolute length of interval                     */
  int  lengthsgn;        /* orientation of the Intervalls                   */
  double step;           /* enlarges interval til sign change found         */
  double dx;             /* RF-stepsize                                     */
  int count = 0;         /* counter for  "no sign change"                   */
  int i;
    

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_NINV_GEN,0.);

  /* initialize starting interval */
   x1 =  GEN.s[0];      /* left boudary of interval */
   x2 =  GEN.s[1];      /* right boudary of interval*/


  if (x1>=x2) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); return 0.; }

 
  /* compute c.d.f. at interval boundaries */
  f1 = CDF(x1) - u;
  f2 = CDF(x2) - u;

  /* search for interval with changing signs */
  step = 1.;     /* interval too small -> make it bigger ( + 2^n * gap ) */ 
  while ( f1*f2 > 0. ) {
    if ( f1 > 0. ) { /* lower boundary too big */    
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1 = CDF(x1) - u;
    }
    else {         /* upper boundary too small */
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2 = CDF(x2) - u;
    }

    /* increase step width */
    step *= 2.;
  }  /* while end -- interval found */ 


  a     = x1;       /* a und x2 soll immer ZW enthalten */

  /* Sekantensuche, ZW wird beibehalten */
  for (i=0; i < GEN.max_iter; i++) { /* noch zu andernde schleife */
    count++;
     
    /* f2 immer kleiner (besser), notfalls Tausch */
    if ( fabs(f1) < fabs(f2) ) {   /* Dreieckstausch */
      xtmp = x1; ftmp = f1;
      x1 = x2;   f1 = f2;
      x2 = xtmp; f2 = ftmp;
    }

    x2abs = fabs(x2);   /* absolute value of x2 */

    
    if ( f2 == 0. )
      /* genauer Treffer  */
      break; /* -> finished */


    if ( f1*f2 <= 0) {  /* ZeichenWechsel vorhanden             */
      count = 0;   /* zaehler fuer bisektion wird rueckgestellt */
      a    = x1;   /* [a, x2] enthaelt ZW                       */
    }
    
    length = x2 - a;  /* gerichtete laenge */
    lengthabs = fabs(length);
    lengthsgn = (length < 0.) ? -1. : 1.;
    
    if ( lengthabs <= GEN.rel_x_resolution * x2abs  )
      /* relative x-genauigkeit erreicht -> finished */
      break; /* -> finished */

  
    /* Sekanten-Schritt  oder Bisektion */
    dx = ( f1-f2==0. ) ? length/2. : f2*(x2-x1)/(f2-f1) ;  
    
    /* minimaler schritt */
    if ( fabs(dx) < GEN.rel_x_resolution * x2abs ){
      dx = lengthsgn * 0.99 * GEN.rel_x_resolution * x2abs;
      while (x2 == x2 - dx){ /* dx zu klein */
	if ( dx != 2.*dx)    /* am genauigkeits-limit des rechners */
	  dx = 2.*dx;
        else
	  dx = length/2.;    /* Bisektion */
      }
    }

       
    /* Bisektionsschritt, wenn:                             */  
    /* kein  ZW  || Schritt fuhrt aus Intervall             */
    if ( count > 1 || 
        (lengthabs-GEN.rel_x_resolution*x2abs)/(dx*lengthsgn) <= 1. )
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
  double x;           /* point for netwon-iteration                   */
  double U;           /* sampled uniform random number                */
  double fx;          /* cdf at x                                     */
  double dfx;         /* pdf at x                                     */
  double fxabs;       /* absolute valuo of fx                         */
  double xtmp, fxtmp; /* temprary variables for x and fx              */
  double xold, fxold; /* remember last values for stopping criterion  */
  double fxtmpabs;    /* fabs of fxtmp                                */
  double damp;        /* damping factor                               */
  double step;        /* helps to escape from flat regions of the cdf */
  int i;              /* counter for for-loop                         */
    
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_NINV_GEN,0.);

  damp = 2.;        /* to be halved at least once */  
  step = 1.;

  /* sample from uniform pseudorandom number generator */
  U = _unur_call_urng(gen);

  /* initialize starting point */
  x     = 1.;    /** TODO: durch s[0] ersetzen ?? **/
  fx    = CDF(x) - U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    /* there is no old value yet */
  fxold = fx;   /* there is no old value yet */ 
  
  /* begin for-loop:  newton-iteration  */
  for (i=0; i < GEN.max_iter; i++) {

    while (dfx == 0.) {   /* function flat at x */
      
      if (fx == 0.)  /* exact hit -> leave while-loopt */
	break; 

      if (fx > 0.)         /* try another x */
        xtmp  = x - step;   
      else 
        xtmp  = x + step;
         
      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);

      if ( fxtmpabs < fxabs ){        /* improvement, update x            */
        step = 1.;     /* set back stepsize */
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmpabs*fxabs < 0. ){/*step was too large, dont update x */
        step /= 2.;                      
      } 
      else{                           /* step was too short, update x     */
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  

      dfx   = PDF(x);
      fxabs = fabs(fx);     
    }   /* end of while-loop, (flat region left) */

   step = 1.;   /* set back stepsize */

   if (fx == 0.)  /* exact hit -> finished */
     break;


    do{    /* newton-step  (damped if nececcary) */
        damp /= 2.;
        xtmp = x - damp * fx/dfx;
        fxtmp = CDF(xtmp) - U;
    } while ( fabs(fxtmp)-fxabs >= fxabs * GEN.rel_x_resolution ); /* no improvement */

    
    /* updation variables according to newton-step      */
    damp  = 2.;       /* set back factor for damping    */
    xold  = x;        /* remember last value of x       */
    fxold = fx;       /* remember last value of fx      */
    x     = xtmp;     /* update x                       */
    fx    = fxtmp;    /* update function value at x     */
    dfx   = PDF(x);   /* update derivative sof fx at x  */
    fxabs = fabs(fx); /* update absolute value of fx    */
 

    /* stopping criterion */
    if ( fabs(x-xold) <= fabs(x) * GEN.rel_x_resolution )        
      break;   /* no improvement with newton-step -> finished */

  }  /* end of for-loop  (MAXITER reached -> finished) */

  return x;

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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

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
  GEN.s[0] = PAR.s[0];     /* staring points                                 */
  GEN.s[1] = PAR.s[1];
  GEN.s[2] = PAR.s[2];

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

  /* check arguments */
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_NINV_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

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

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d interations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN.max_iter);

} /* end of _unur_ninv_debug_sample_regula() */

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/






















