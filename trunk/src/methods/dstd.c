/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dstd.c                                                       *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
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
 * DSTD is a wrapper for special generator for Discrete univariate           *
 * STandarD distributions. It only works for distributions in the UNURAN     *
 * library of distributions and will refuse to work otherwise.               *
 * (In detail it rejects a distribution if its id is equal to DISTR_GENERIC, *
 * the id inserted my the unur_distr_discr_new() call.)                      *
 *                                                                           *
 * It calls the initialzation routine provided by the distribution object.   *
 * This routine has to do all setup steps for the special generator.         *
 * if no such routine is given, i.e. distr->init==NULL, then unur_dstd_new() *
 * does not work and the NULL pointer is returned instead of the pointer to  *
 * a parameter object.                                                       *
 *                                                                           *
 * Notice that using a truncated distribution (this can be constructed by    *
 * changing the default domain of a distribution by means of an              *
 * unur_distr_cont_set_domain() call) is not allowed.                        *
 *                                                                           *
 * Variants (different algorithms for the same distribution) are possible    *
 * and can be selected by unsigned integers using the                        *
 * unur_dstd_set_variant() call.                                             *
 * For possible variants see the generator files for each distribution in    *
 * the distributions directory. However the following are common to all      *
 * distributions:                                                            *
 *    0                     ... the default generator                        *
 *    UNUR_STDGEN_INVERSION ... the inversion method (if available)          *
 * unur_dstd_set_variant() return 0 if a variant is not implemented, and 1   *
 * otherwise. In the first case the selected variant is not changed.         *
 *                                                                           *
 * It is possible to change the parameters of the chosen distribution        *
 * without building a new generator object by means of the                   *
 * unur_dstd_chg_pmfparams() call.                                           *
 * Notice that it is not possible to change the number of parameters.        *
 * This function only copies the given arguments into the array of           *
 * parameters.                                                               *
 * IMPORTANT: The given parameters are not checked against domain errors     *
 * (as the unur_<distr>_new() calls do).                                     *
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

#define DSTD_DEBUG_CHG          0x00001000u   /* print changed parameters    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DSTD_SET_VARIANT          0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DSTD"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dstd_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are no sampling routines, since every distribution has its own.     */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_dstd_sample( UNUR_GEN *gen ); does not exist!                */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_debug_chg_pmfparams( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) parameters of distribution                            */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dstd_new( struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

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
  COOKIE_SET(par,CK_DSTD_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  PAR.sample_routine_name = NULL ;  /* name of sampling routine              */
  PAR.is_inversion = FALSE;         /* method not based on inversion         */

  par->method   = UNUR_METH_DSTD;   /* method                                */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_dstd_init;

  return par;

} /* end of unur_dstd_new() */

/*****************************************************************************/

int 
unur_dstd_set_variant( struct unur_par *par, unsigned variant )
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
  _unur_check_NULL( GENTYPE,par,0 );
  _unur_check_NULL( GENTYPE,par->distr,0 );

  /* check input */
  _unur_check_par_object( par,DSTD );

  /* store date */
  old_variant = par->variant;
  par->variant = variant;

  /* check variant. run special init routine only in test mode */
  if (par->DISTR_IN.init != NULL && par->DISTR_IN.init(par,NULL) ) {
    par->set |= DSTD_SET_VARIANT;    /* changelog */
    return 1;
  }

  /* variant not valid */
  _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
  par->variant = old_variant;
  return 0;

} /* end of unur_dstd_set_variant() */

/*---------------------------------------------------------------------------*/

int 
unur_dstd_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
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
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSTD );
  if (n_params>0) CHECK_NULL(params,0);
  
  /* set new parameters in distribution object */
  if (!unur_distr_discr_set_pmfparams(&(gen->distr),params,n_params))
    return 0;

  /* run special init routine for generator */
  if ( !DISTR.init(NULL,gen) ) {
    /* init failed --> could not find a sampling routine */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return 0;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug & DSTD_DEBUG_CHG) 
      _unur_dstd_debug_chg_pmfparams( gen );
#endif

  /* o.k. */
  return 1;

} /* end of unur_dstd_chg_pmfparams() */

/*****************************************************************************/

struct unur_gen *
_unur_dstd_init( struct unur_par *par )
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
  _unur_check_NULL( GENTYPE,par->DISTR_IN.init,NULL );

  /* check input */
  if ( par->method != UNUR_METH_DSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_DSTD_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dstd_create(par);
  if (!gen) { free(par); return NULL; }

  /* check for initializing routine for special generator */
  _unur_check_NULL( gen->genid, DISTR.init, (free(par),NULL) );

  /* run special init routine for generator */
  if ( !DISTR.init(par,gen) ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    free(par); _unur_dstd_free(gen); return NULL; 
  }

  /* domain must not be changed */
  if (!(gen->distr.set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* domain has been modified --> not allowed */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed");
    free(par); _unur_dstd_free(gen); return NULL; 
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dstd_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_dstd_init() */

/*****************************************************************************/

/** 
    double _unur_dstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

/*****************************************************************************/

void
_unur_dstd_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  /* check input */
  if ( gen->method != UNUR_METH_DSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  if (DISTR.pv)        free( DISTR.pv );
  if (GEN.gen_param)   free(GEN.gen_param);
  if (GEN.gen_iparam)  free(GEN.gen_iparam);
  if (gen->gen_aux)   _unur_free(gen->gen_aux);
  if (gen->gen_aux_2) _unur_free(gen->gen_aux_2);
  free(gen);

} /* end of _unur_dstd_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_dstd_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSTD_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DSTD_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;    /* will be set in _unur_dstd_init() */
  gen->destroy = _unur_dstd_free;

  /* defaults */
  GEN.gen_param = NULL;  /* parameters for the generator      */
  GEN.n_gen_param = 0;   /* (computed in special GEN.init()   */
  GEN.gen_iparam = NULL; /* smake for integer parameters      */
  GEN.n_gen_iparam = 0;

  /* copy some parameters into generator object */
  GEN.umin        = 0;    /* cdf at left boundary of domain   */
  GEN.umax        = 1;    /* cdf at right boundary of domain  */

  gen->method = par->method;        /* indicates used method  */
  gen->variant = par->variant;      /* indicates variant      */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags         */
  gen->urng = par->urng;            /* pointer to urng        */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_dstd_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dstd_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_DSTD_PAR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_discr_debug( &(gen->distr), gen->genid, FALSE );

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

} /* end of _unur_dstd_info_init() */

/*---------------------------------------------------------------------------*/

static void 
_unur_dstd_debug_chg_pmfparams( struct unur_gen *gen )
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:parameters of distribution changed:\n",gen->genid);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\tparam[%d] = %g\n",gen->genid,i,DISTR.params[i]);

} /* end of _unur_dstd_debug_chg_pmfparams() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
