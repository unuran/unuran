/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      stdr.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection with universal bounds          *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and mode of a T_{-1/2}-concave distribution              *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      mode of the density                                                  *
 *      area below p.d.f.                                                    *
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
 *   [1] Leydold J. (2000): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions, preprint.              *
 *                                                                           *
 *   [2] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * This algorithm is based on transformed density rejection (see [2]), were  *
 * universal upper and lower bounds are used. These are derived via the      *
 * ratio-of-uniforms method. See [1] for details and a description of the    *
 * algorithm. It works for any distribution, where -1/sqrt(pdf(x)) is        *
 * concave. This includes all log-concave distributions.                     *
 *                                                                           *
 * (The mirror principle has not been implemented.)                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define STDR_VARFLAG_VERIFY   0x002u   /* run verify mode                    */
#define STDR_VARFLAG_SQUEEZE  0x004u   /* use universal squeeze if possible  */

#define STDR_VARFLAG_MIRROR   0x008u   /* use mirror principle               */
/* not implemented yet! */

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define STDR_SET_FMODE        0x001u   /* cdf at mode is known               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "STDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_stdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_stdr_debug_init( struct unur_par *par, struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.stdr        /* data for parameter object         */
#define GEN       gen->data.stdr        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x) ((*(DISTR.pdf))((x),DISTR.params,DISTR.n_params))    /* call to p.d.f. */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2    1.4142135623731

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_stdr_new( struct unur_distr *distr )
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
  CHECK_NULL(distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.d.f.");
    return NULL;
  }
  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode of p.d.f.");
    return NULL;
  }
  if (!(distr->set & UNUR_DISTR_SET_PDFAREA)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area of p.d.f.");
    return NULL;
  }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_STDR_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR.Fmode        = -1.;     /* c.d.f. at mode (unknown yet)                */

  par->method      = UNUR_METH_STDR;  /* method and default variant          */
  par->variant     = 0u;              /* default variant                     */
  par->set         = 0u;              /* inidicate default parameters        */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  par->debug       = UNUR_DEBUGFLAG_DEFAULT;  /* set default debugging flags */

  /* routine for starting generator */
  par->init = unur_stdr_init;

  return par;

} /* end of unur_stdr_new() */

/*****************************************************************************/

int 
unur_stdr_set_Fmode( struct unur_par *par, double Fmode )
     /*----------------------------------------------------------------------*/
     /* set cdf at mode                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   Fmode ... cdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( STDR );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"cdf(mode) out of range");
    return 0;
  }

  /* store date */
  PAR.Fmode = Fmode;

  /* changelog */
  par->set |= STDR_SET_FMODE;

  return 1;

} /* end of unur_stdr_set_Fmode() */

/*---------------------------------------------------------------------------*/

int
unur_stdr_set_verify( struct unur_par *par, int verify )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( STDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | STDR_VARFLAG_VERIFY) : (par->variant & (~STDR_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_stdr_set_verify() */

/*---------------------------------------------------------------------------*/

int 
unur_stdr_set_usesqueeze( struct unur_par *par, int usesqueeze )
     /*----------------------------------------------------------------------*/
     /* set flag for using universal squeeze (default: off)                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   usesqueeze ... 0 = no squeeze,  !0 = use squeeze                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no squeeze is the default                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( STDR );

  /* we use a bit in variant */
  par->variant = (usesqueeze) ? (par->variant | STDR_VARFLAG_SQUEEZE) : (par->variant & (~STDR_VARFLAG_SQUEEZE));

  /* o.k. */
  return 1;

} /* end of unur_stdr_set_usesqueeze() */

/*****************************************************************************/

struct unur_gen *
unur_stdr_init( struct unur_par *par )
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
  double vm;               /* width of rectangle                             */
  double left,right;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_STDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_STDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_stdr_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute pdf at mode */
  GEN.fm = PDF(DISTR.mode);

  /* fm must be positive */
  if (GEN.fm <= 0.) {
    _unur_error(GENTYPE,UNUR_ERR_INIT,"pdf(mode) <= 0.");
    free(par); unur_stdr_free(gen);
    return NULL;
  }

  /* compute parameters */
  GEN.um = sqrt(GEN.fm);
  vm = DISTR.area / GEN.um;

  if (par->set & STDR_SET_FMODE) {
    /* cdf at mode known */
    GEN.vl = -PAR.Fmode * vm;
    GEN.vr = vm + GEN.vl;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
    GEN.A  = 2 * DISTR.area;
    GEN.al = (DISTR.BD_LEFT  < DISTR.mode) ? (PAR.Fmode * DISTR.area) : 0.;
    GEN.ar = (DISTR.BD_RIGHT > DISTR.mode) ? (GEN.al + DISTR.area) : GEN.A;
    /* Compute areas below hat in left tails and inside domain of pdf */
    if ( (DISTR.BD_LEFT > -INFINITY) &&
	 (DISTR.BD_LEFT < DISTR.mode) )
      GEN.Aleft = GEN.vl * GEN.vl / (DISTR.mode - DISTR.BD_LEFT);
    else
      GEN.Aleft = 0.;
    
    if ( (DISTR.BD_RIGHT < INFINITY) &&
	 (DISTR.BD_RIGHT > DISTR.mode) )
      GEN.Ain = GEN.A - GEN.vr * GEN.vr / (DISTR.BD_RIGHT - DISTR.mode);
    else
      GEN.Ain = GEN.A;
    GEN.Ain -= GEN.Aleft;

    /* it does not make sense to use the mirror principle */
    /* gen->variant &= ~STDR_VARFLAG_MIRROR;  not implemented */
  }

  else {
    /* cdf at mode unknown */
    GEN.vl = -vm;
    GEN.vr = vm;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
    GEN.A  = 4 * DISTR.area;
    GEN.al = DISTR.area;
    GEN.ar = 3 * DISTR.area;
    /* Compute areas below hat in left tails and inside domain of pdf */
    if (DISTR.BD_LEFT > -INFINITY) {
      left = DISTR.BD_LEFT - DISTR.mode;
      GEN.Aleft = (GEN.xl > left) 
	? (GEN.vl * GEN.vl / (-left)) 
	: (GEN.al + GEN.fm * (left - GEN.xl));
    }
    else 
      GEN.Aleft = 0.;
    
    if (DISTR.BD_RIGHT < INFINITY) {
      right = DISTR.BD_RIGHT - DISTR.mode;
      GEN.Ain = (GEN.xr < right) 
	? (GEN.A - GEN.vr * GEN.vr / right)
	: (GEN.ar - GEN.fm * (GEN.xr - right));
    }
    else 
      GEN.Ain = GEN.A;
    GEN.Ain -= GEN.Aleft;
    
    /* we cannot use universal squeeze */
    gen-> variant &= ~STDR_VARFLAG_SQUEEZE;
  }

#if UNUR_DEBUG & UNUR_DB_INFO
    /* write info into log file */
    if (gen->debug) _unur_stdr_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_stdr_init() */

/*****************************************************************************/

double
unur_stdr_sample( struct unur_gen *gen )
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
  double u,v,x,xx,y;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_STDR_GEN,0.);

  while (1) {
    /* uniform ~U(0,1) */
    u = GEN.Aleft + _unur_call_urng(gen) * GEN.Ain;
    if (u==0.) continue;

    if (u < GEN.al) {        /* first part */
      x = - GEN.vl * GEN.vl / u;
      y = (u / GEN.vl);
      y = y*y;
    }
    else if (u <= GEN.ar) {  /* second part */
      x = GEN.xl + (u-GEN.al)/GEN.fm;
      y = GEN.fm;
    }
    else {                   /* third part */
      x = GEN.vr * GEN.vr / (GEN.um * GEN.vr - (u-GEN.ar));
      y = (GEN.A - u) / GEN.vr;
      y = y*y;
    }

    /* accept or reject */
    v = _unur_call_urng(gen);
    y *= v;

    /* evaluate squeeze */
    if (gen->variant & STDR_VARFLAG_SQUEEZE) {
      xx = 2 * x;
      if ( x >= GEN.xl && x <= GEN.xr && y <= GEN.fm/4. )
	return (x + DISTR.mode);
    }

    /* Compute X */
    x += DISTR.mode;

    /* evaluate p.d.f. */
    if (y <= PDF(x))
      return x;
  }
    
} /* end of unur_stdr_sample() */

/*---------------------------------------------------------------------------*/

double
unur_stdr_sample_check( struct unur_gen *gen )
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
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  double u,v,x,xx,fx,y;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_STDR_GEN,0.);

  while (1) {
    /* uniform ~U(0,1) */
    u = GEN.Aleft + _unur_call_urng(gen) * GEN.Ain;
    if (u==0.) continue;

    if (u < GEN.al) {        /* first part */
      x = - GEN.vl * GEN.vl / u;
      y = (u / GEN.vl);
      y = y*y;
    }
    else if (u <= GEN.ar) {  /* second part */
      x = GEN.xl + (u-GEN.al)/GEN.fm;
      y = GEN.fm;
    }
    else {                   /* third part */
      x = GEN.vr * GEN.vr / (GEN.um * GEN.vr - (u-GEN.ar));
      y = (GEN.A - u) / GEN.vr;
      y = y*y;
    }

    /* compute pdf at x */
    fx = PDF(x + DISTR.mode);

    /* verify hat function */
    if (fx > y)
      _unur_error(GENTYPE,UNUR_ERR_SAMPLE,"f(x) > h(x)");

    /* accept or reject */
    v = _unur_call_urng(gen);
    y *= v;

    /* evaluate squeeze */
    if (gen->variant & STDR_VARFLAG_SQUEEZE) {
      xx = 2 * x;
      if ( x >= GEN.xl && x <= GEN.xr && y <= GEN.fm/4. )
	return (x + DISTR.mode);
    }

    /* Compute X */
    x += DISTR.mode;

    /* evaluate p.d.f. */
    if (y <= fx)
      return x;
  }
    
} /* end of unur_stdr_sample_check() */

/*****************************************************************************/

void
unur_stdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_STDR ) {
    _unur_warning(GENTYPE,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_STDR_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  free(gen->genid);
  free(gen);

} /* end of unur_stdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_stdr_create( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_STDR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_STDR_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_make_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & STDR_VARFLAG_VERIFY) ? unur_stdr_sample_check : unur_stdr_sample;
  gen->destroy = unur_stdr_free;

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_INIT,"area and cdf at mode might be wrong");
    DISTR.mode = max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = min(DISTR.mode,DISTR.BD_RIGHT);
  }

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  /* initialize parameters */

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_stdr_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
_unur_stdr_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  fprintf(log,"%s: method  = stdr (simple universal transformed density rection)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = unur_stdr_sample",gen->genid);
  if (par->variant & STDR_VARFLAG_VERIFY)
    fprintf(log,"_check");
  /* else if (par->variant & STDR_VARFLAG_MIRROR)     not implemented */
  /*   fprintf(log,"_mirror"); */
  fprintf(log,"()\n%s:\n",gen->genid);

  if (par->set & STDR_SET_FMODE)
    fprintf(log,"%s: c.d.f. at mode = %g\n",gen->genid,PAR.Fmode);
  else
    fprintf(log,"%s: c.d.f. at mode unknown\n",gen->genid);

  if (gen->variant & STDR_VARFLAG_SQUEEZE)
    fprintf(log,"%s: use universal squeeze\n",gen->genid);
  else
    fprintf(log,"%s: no (universal) squeeze\n",gen->genid);

  if (gen->variant & STDR_VARFLAG_MIRROR)
    fprintf(log,"%s: use mirror principle\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: parts:\n",gen->genid);
  fprintf(log,"%s:\txl = %g\n",gen->genid,GEN.xl);
  fprintf(log,"%s:\txr = %g\n",gen->genid,GEN.xr);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: pdf at mode:\n",gen->genid);
  fprintf(log,"%s:\tfm = %g\n",gen->genid,GEN.fm);
  fprintf(log,"%s:\tum = %g\n",gen->genid,GEN.um);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: areas:\n",gen->genid);
  fprintf(log,"%s:\t    al = %g\n",gen->genid,GEN.al);
  fprintf(log,"%s:\t    ar = %g\n",gen->genid,GEN.ar);
  fprintf(log,"%s:\t Aleft = %g\n",gen->genid,GEN.Aleft);
  fprintf(log,"%s:\t   Ain = %g\n",gen->genid,GEN.Ain);
  fprintf(log,"%s:\tAtotal = %g\n",gen->genid,GEN.A);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_stdr_debug_init() */

#endif

/*****************************************************************************/

