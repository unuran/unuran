/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cstd.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generators for standard distribution                         *
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
 * CSTD is a wrapper for special generator for Continuous univariate         *
 * STandarD distributions. It only works for distributions in the UNURAN     *
 * library of distributions and will refuse to work otherwise.               *
 * (In detail it rejects a distribution if its id is equal to DISTR_GENERIC, *
 * the id inserted my the unur_distr_cont_new() call.)                       *
 *                                                                           *
 * It calls the initialzation routine provided by the distribution object.   *
 * This routine has to do all setup steps for the special generator.         *
 * if no such routine is given, i.e. distr->init==NULL, then unur_cstd_new() *
 * does not work and the NULL pointer is returned instead of the pointer to  *
 * a parameter object.                                                       *
 *                                                                           *
 * Notice that using a truncated distribution (this can be constructed by    *
 * changing the default domain of a distribution by means of an              *
 * unur_distr_cont_set_domain() call) is only allowed if the inversion       *
 * method is used. Otherwise no parameter object is returned by the          *
 * unur_cstd_new() call.                                                     *
 *                                                                           *
 * Variants (different algorithms for the same distribution) are possible    *
 * and can be selected by unsigned integers using the                        *
 * unur_cstd_set_variant() call.                                             *
 * For possible variants see the generator files for each distribution in    *
 * the distributions directory. However the following are common to all      *
 * distributions:                                                            *
 *                                                                           *
 *    UNUR_STDGEN_DEFAULT   ... the default generator                        *
 *    UNUR_STDGEN_INVERSION ... the inversion method (if available)          *
 *    UNUR_STDGEN_FAST      ... the fasted available special generator       *
 *                                                                           *
 * unur_cstd_set_variant() return 0 if a variant is not implemented, and 1   *
 * otherwise. In the first case the selected variant is not changed.         *
 *                                                                           *
 * It is possible to change the parameters of the chosen distribution        *
 * without building a new generator object by means of the                   *
 * unur_cstd_chg_pdfparams() call.                                           *
 * Notice that it is not possible to change the number of parameters.        *
 * This function only copies the given arguments into the array of           *
 * parameters.                                                               *
 * IMPORTANT: The given parameters are not checked against domain errors     *
 * (as the unur_<distr>_new() calls do).                                     *
 *                                                                           *
 * The domain of a (truncated) distribution can be changed without building  *
 * a new generator object by means of the unur_cstd_chg_truncated() call.    *
 * Notice that this only works when the inversion method is used. Otherwise  *
 * nothing happens to the domain and an error message is produced.           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define CSTD_DEBUG_CHG          0x00001000u   /* print changed parameters    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define CSTD_SET_VARIANT          0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "CSTD"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are sampling routines, since every distribution has its own.        */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_cstd_sample( UNUR_GEN *gen ); does not exist!                */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_chg_pdfparams( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) parameters of distribution                            */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_chg_truncated( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) domain of (truncated) distribution                    */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define CDF(x)    _unur_cont_CDF((x),&(gen->distr))   /* call to c.d.f.            */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_cstd_new( struct unur_distr *distr )
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
  _unur_check_NULL(GENTYPE,distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (distr->id == UNUR_DISTR_GENERIC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"standard distribution");
    return NULL;
  }
  if (DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators");
    return NULL;
  }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_CSTD_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  PAR.sample_routine_name = NULL ;  /* name of sampling routine              */
  PAR.is_inversion = FALSE;         /* method not based on inversion         */

  par->method   = UNUR_METH_CSTD;   /* method                                */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_cstd_init;

  return par;

} /* end of unur_cstd_new() */

/*****************************************************************************/

int 
unur_cstd_set_variant( struct unur_par *par, unsigned variant )
     /*----------------------------------------------------------------------*/
     /* set variant of method                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   variant ... indicator for variant                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  unsigned old_variant;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, 0 );
  _unur_check_NULL( GENTYPE, par->distr, 0 );

  /* check input */
  _unur_check_par_object( par,CSTD );

  /* store date */
  old_variant = par->variant;
  par->variant = variant;

  /* check variant. run special init routine only in test mode */
  if (par->DISTR_IN.init != NULL && par->DISTR_IN.init(par,NULL) ) {
    par->set |= CSTD_SET_VARIANT;    /* changelog */
    return 1;
  }

  /* variant not valid */
  _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
  par->variant = old_variant;
  return 0;

} /* end if unur_cstd_set_variant() */

/*---------------------------------------------------------------------------*/

int 
unur_cstd_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,CSTD );
  if (n_params>0) CHECK_NULL(params,0);
  
  /* check new parameter for generator */
  if (n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* copy parameters */
  DISTR.n_params = n_params;
  memcpy( DISTR.params, params, n_params*sizeof(double) );

  /* changelog */
  /* derived parameters like mode, area, etc. might be wrong now! */
  /* however there is no need for these parameters;               */
  /* thus we do not need:                                         */
  /*    distr->distr.set &= ~UNUR_DISTR_SET_MASK_DERIVED;         */

  /* run special init routine for generator */
  if ( !(DISTR.init(NULL,gen)) ) {
    /* init failed --> could not find a sampling routine */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return 0;
  }

  if ( GEN.is_inversion )
    if( gen->distr.set & UNUR_DISTR_SET_TRUNCATED ) {
      /* truncated domain */
      if (DISTR.cdf == NULL) {
	_unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); return 0; }
      /* compute umin and umax */
      GEN.umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
      GEN.umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug & CSTD_DEBUG_CHG) 
      _unur_cstd_debug_chg_pdfparams( gen );
#endif

  /* o.k. */
  return 1;

} /* end of unur_cstd_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_cstd_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the               */
     /* (truncated) distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,CSTD );

  /* domain can only be changed for inversion method! */
  if ( ! GEN.is_inversion ) { 
    /* this is not the inversion method */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"truncated domain for non inversion method");
    return 0;
  }

  /* c.d.f. required ! */
  if (DISTR.cdf == NULL) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"truncated domain, c.d.f. required");
    return 0;
  }

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    DISTR.trunc[0] = DISTR.domain[0];
  }
  else
    DISTR.trunc[0] = left;

  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    DISTR.trunc[1] = DISTR.domain[1];
  }
  else
    DISTR.trunc[1] = right;

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_TRUNCATED;

  /* indicate that we have a truncated distribution.
     (do not have the standard domain any more) */
  gen->distr.set &= ~UNUR_DISTR_SET_STDDOMAIN;

  /* compute umin and umax */
  GEN.umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN.umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & CSTD_DEBUG_CHG) 
    _unur_cstd_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_cstd_chg_truncated() */

/*****************************************************************************/

struct unur_gen *
_unur_cstd_init( struct unur_par *par )
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
  CHECK_NULL(par,NULL);
  _unur_check_NULL( GENTYPE, par->DISTR_IN.init, NULL );

  /* check input */
  if ( par->method != UNUR_METH_CSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);
  
  /* create a new empty generator object */
  gen = _unur_cstd_create(par);
  if (!gen) { free(par); return NULL; }
  
  /* check for initializing routine for special generator */
  _unur_check_NULL( gen->genid, DISTR.init, (free(par),NULL) );
  
  /* reset flag for inversion method */
  PAR.is_inversion = FALSE;

  /* run special init routine for generator */
  if ( !DISTR.init(par,gen) ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    free(par); _unur_cstd_free(gen); return NULL; 
  }
  
  /* copy information about type of special generator */
  GEN.is_inversion = PAR.is_inversion;
  
  /* domain valid for special generator ?? */
  if (!(gen->distr.set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* domain has been modified */
    gen->distr.set &= UNUR_DISTR_SET_TRUNCATED;
    DISTR.trunc[0] = DISTR.domain[0];
    DISTR.trunc[1] = DISTR.domain[1];

    if ( ! GEN.is_inversion ) { 
      /* this is not the inversion method */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed for non inversion method");
      free(par); _unur_cstd_free(gen); return NULL; 
    }

    if (DISTR.cdf == NULL) {
      /* using a truncated distribution requires a CDF */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed, c.d.f. required");
      free(par); _unur_cstd_free(gen); return NULL; 
    }

    /* compute umin and umax */
    GEN.umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
    GEN.umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_cstd_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_cstd_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_cstd_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_CSTD_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;      /* will be set in _unur_cstd_init() */
  gen->destroy = _unur_cstd_free;

  /* defaults */
  GEN.gen_param = NULL;  /* parameters for the generator                     */
  GEN.n_gen_param = 0;   /* (computed in special GEN.init()                  */

  /* copy some parameters into generator object */
  GEN.umin        = 0;    /* cdf at left boundary of domain                  */
  GEN.umax        = 1;    /* cdf at right boundary of domain                 */

  /* GEN.is_inversion is set in _unur_cstd_init() */

  gen->method = par->method;        /* indicates used method                 */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_cstd_create() */

/*****************************************************************************/

/** 
    double _unur_cstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

/*****************************************************************************/

void
_unur_cstd_free( struct unur_gen *gen )
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

  /* magic cookies */
  COOKIE_CHECK(gen,CK_CSTD_GEN,/*void*/);

  /* check input */
  if ( gen->method != UNUR_METH_CSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.gen_param);
  if (gen->gen_aux)   _unur_free(gen->gen_aux);
  if (gen->gen_aux_2) _unur_free(gen->gen_aux_2);
  free(gen);

} /* end of _unur_cstd_free() */

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
_unur_cstd_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_CSTD_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_CSTD_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  /* sampling routine */
  fprintf(log,"%s: sampling routine = ",gen->genid);
  if (PAR.sample_routine_name)
    fprintf(log,"%s()",PAR.sample_routine_name);
  else
    fprintf(log,"(Unknown)");
  if (PAR.is_inversion)
    fprintf(log,"   (Inversion)");
  fprintf(log,"\n%s:\n",gen->genid);

  if (!(par->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    fprintf(log,"%s: domain has been changed. U in (%g,%g)\n",gen->genid,GEN.umin,GEN.umax);
    fprintf(log,"%s:\n",gen->genid);
  }

} /* end of _unur_cstd_info_init() */

/*---------------------------------------------------------------------------*/

static void 
_unur_cstd_debug_chg_pdfparams( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) parameters of distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_CSTD_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: parameters of distribution changed:\n",gen->genid);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\tparam[%d] = %g\n",gen->genid,i,DISTR.params[i]);
  if (gen->distr.set & UNUR_DISTR_SET_TRUNCATED)
    fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.umin,GEN.umax);

} /* end of _unur_cstd_debug_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

static void 
_unur_cstd_debug_chg_truncated( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_CSTD_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: domain of truncated distribution changed:\n",gen->genid);
  fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.umin,GEN.umax);

} /* end of _unur_cstd_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

