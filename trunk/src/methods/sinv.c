/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      sinv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    Spline approximation for INVerse of CDF                      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the CDF and PDF                                           *
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
 *   [1] Neumaier A. (to be published): Introduction to numerical analysis,  *
 *       Cambridge University Press                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Numerical inversion is a method for generating random variables          *
 *  using the CDF (and in case of newton's method the PDF).                  *
 *                                                                           *
 *  THEOREM:                                                                 *
 *     Let X be a random variable with CDF F(x).                             *
 *     Then the F(X) are  uniformly distributed.                             *
 *                                                                           *
 *  COROLLARY:                                                               *
 *     Starting with uniformly distributed random variables U,               *
 *     the F^(-1)(U) have F(x) as CDF.                                       *
 *                                                                           *
 *  Starting with an U, the task is to find a X fulfilling:                  *
 *    F(X) - U = 0.                                                          *
 *                                                                           *
 *  Numerical algorithms to find zeros that are used in SINV are variants of * 
 *  newton's method (damped newton to guarantee improvement) and             *
 *  the regula falsi ( stabilized regula falsi preserving sign change; at    *
 *  first an interval with sign change is determined).                       *
 *                                                                           *
 *  In both cases it is possible to specify the maximal number of            *
 *  iterations, a desired accuracy in X and starting values for the          *
 *  algorithms.                                                              *
 *  Instead of starting values it is also possible to use a table            *
 *  containing suitable starting values.                                     *
 *  If neither the table nor explicit starting values are used,              *
 *  SINV chooses as starting values:                                         *
 *     newton's method:  x:     CDF(x) = 0.5                                 *
 *     regula falsi:     x1,x2: CDF(x1) = 1 - CDF(x2) = 0.05                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

#define SINV_DEFAULT_COMPUTATION_LIMIT  1.e20

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define SINV_DEBUG_TABLE        0x00000010u   /* print table                 */
#define SINV_DEBUG_CHG          0x00001000u   /* print changed parameters    */
#define SINV_DEBUG_SAMPLE       0x01000000u   /* trace sampling              */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define SINV_SET_BOUNDARY       0x001u  /* boundary for comput. interval     */
#define SINV_SET_U_RESOLUTION   0x002u  /* maximal error in u                */

/*---------------------------------------------------------------------------*/

#define GENTYPE "SINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_sinv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_sinv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_sinv_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create the table with splines                                             */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_splines( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into logfile.               */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_chg_truncated( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* trace changes of the truncated domain.                                    */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.sinv        /* data for parameter object         */
#define GEN       gen->data.sinv        /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */
#define CDF(x)    _unur_cont_CDF((x),(gen->distr))    /* call to CDF         */

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_sinv_new( const struct unur_distr *distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }

  /* if default variant is Newton's method, then we also need the PDF ! */

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_SINV_PAR);

  /* copy input */
  par->distr   = distr;           /* pointer to distribution object          */

  /* set default values */
  PAR.u_resolution = 1.0e-8;      /* maximal error allowed in u              */

  /* default boundary of compution area */
  PAR.bleft     = -SINV_DEFAULT_COMPUTATION_LIMIT;
  PAR.bright    = SINV_DEFAULT_COMPUTATION_LIMIT;

  par->method   = UNUR_METH_SINV; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_sinv_init;

  return par;

} /* end of unur_sinv_new() */

/*****************************************************************************/

int
unur_sinv_set_boundary( struct unur_par *par, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set left and right boundary of computation interval                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   new boundary points must not be +/- INFINITY                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SINV );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return 0;
  }
  if (left <= -INFINITY || right >= INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return 0;
  }

  /* store date */
  PAR.bleft = left;
  PAR.bright = right;

  /* changelog */
  par->set |= SINV_SET_BOUNDARY;

  return 1;

} /* end of unur_sinv_set_boundary() */

/*---------------------------------------------------------------------------*/

int
unur_sinv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... maximal error in u                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SINV );

  /* check new parameter for generator */
  if (u_resolution < DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u resolution");
    return 0;
  }

  /* store date */
  PAR.u_resolution = u_resolution;

  /* changelog */
  par->set |= SINV_SET_U_RESOLUTION;

  return 1;

} /* end of unur_sinv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int 
unur_sinv_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{

  /*** TODO: kann man das so verwenden ??????? ******/


  double Umin, Umax;

  /* check arguments */
  CHECK_NULL(gen, 0);
  _unur_check_gen_object(gen, SINV);

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* set bounds of U -- in respect to given bounds */
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;

  /* check result */
  if (Umin > Umax) {
    /* this is a serios error that should not happen */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

  if (_unur_FP_equal(Umin,Umax)) {
    /* CDF values very close */
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (Umin == 0. || _unur_FP_same(Umax,1.)) {
      /* this is very bad */
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return 0;
    }
  }

  /* copy new boundaries into generator object */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN.Umin = Umin;
  GEN.Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & SINV_DEBUG_CHG) 
    _unur_sinv_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_sinv_chg_truncated() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

struct unur_gen *
_unur_sinv_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_SINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_sinv_create(par);
  if (!gen) { free(par); return NULL; }

  /* free parameters */
  free(par);
  
  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* set bounds of U -- in respect to given bounds                          */
  GEN.CDFmin = GEN.Umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN.CDFmax = GEN.Umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;

  if (_unur_FP_greater(GEN.CDFmin, GEN.CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    _unur_sinv_free(gen); return NULL;
  }

  /* compute splines */
  if (!_unur_sinv_create_table(gen)) {
    _unur_sinv_free(gen);
    return NULL;
  }
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_sinv_debug_init(gen);
#endif

  return gen;

} /* end of _unur_sinv_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_sinv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SINV_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_SINV_GEN);

  /* copy distribution object into generator object */
  gen->distr = _unur_distr_cont_clone( par->distr );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_sinv_sample;

  gen->destroy = _unur_sinv_free;
  gen->clone = _unur_sinv_clone;

  /* copy parameters into generator object */
  GEN.u_resolution = PAR.u_resolution; /* maximal error in u                 */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */

  /* init pointer */
  GEN.splines = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_sinv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_sinv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a table of splines                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* o.k. */
  return 1;
}  /* end of _unur_sinv_create_table() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_sinv_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.sinv

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SINV_GEN,NULL);

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  clone->distr = _unur_distr_cont_clone( gen->distr );

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = unur_gen_clone( gen->gen_aux );

  /* copy additional data for generator object */

  return clone;

#undef CLONE
} /* end of _unur_sinv_clone() */

/*****************************************************************************/

void
_unur_sinv_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_SINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free tables */

  /* free memory */
  _unur_distr_free(gen->distr);
  _unur_free_genid(gen);

  COOKIE_CLEAR(gen);
  free(gen);

} /* end of _unur_sinv_free() */

/*****************************************************************************/

double
_unur_sinv_sample( struct unur_gen *gen )
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
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SINV_GEN,0.);

  return 0.;

} /* end of _unur_sinv_sample() */

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
_unur_sinv_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = SINV (Spline approximation of INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_sinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_sinv_debug_splines(gen);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_sinv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_sinv_debug_splines( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print starting points or table for algorithms into logfile           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_sinv_debug_splines() */

/*---------------------------------------------------------------------------*/

void 
_unur_sinv_debug_chg_truncated( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.Umin,GEN.Umax);

} /* end of _unur_sinv_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
