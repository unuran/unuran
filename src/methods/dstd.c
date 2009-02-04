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
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dstd.h"
#include "dstd_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DSTD_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */
#define DSTD_DEBUG_CHG       0x00001000u   /* print changed parameters       */

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

static int _unur_dstd_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dstd_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dstd_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are no sampling routines, since every distribution has its own.     */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_dstd_sample( UNUR_GEN *gen ); does not exist!                */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dstd_debug_chg_pmfparams( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) parameters of distribution                            */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dstd_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dstd_new( const struct unur_distr *distr )
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

  if (!(distr->id & UNUR_DISTR_STD) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"standard distribution");
    return NULL;
  }
  if (DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dstd_par) );
  COOKIE_SET(par,CK_DSTD_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  unsigned old_variant;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, par->distr, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSTD );

  /* store date */
  old_variant = par->variant;
  par->variant = variant;

  /* check variant. run special init routine only in test mode */
  if (par->DISTR_IN.init != NULL && par->DISTR_IN.init(par,NULL)==UNUR_SUCCESS ) {
    par->set |= DSTD_SET_VARIANT;    /* changelog */
    return UNUR_SUCCESS;
  }

  /* variant not valid */
  _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
  par->variant = old_variant;
  return UNUR_ERR_PAR_VARIANT;

} /* end of unur_dstd_set_variant() */


/*****************************************************************************/
/**  Private                                                                **/
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
  /* check for required data: initializing routine for special generator */
  if (par->DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,"");
    return NULL;
  }

  /* check input */
  if ( par->method != UNUR_METH_DSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_DSTD_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dstd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_dstd_check_par(gen) != UNUR_SUCCESS) {
    _unur_dstd_free(gen); return NULL;
  }

  /* run special init routine for generator */
  GEN->is_inversion = FALSE;   /* reset flag for inversion method */
  if ( DISTR.init(NULL,gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    _unur_dstd_free(gen); return NULL; 
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_dstd_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_dstd_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dstd_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int rcode;

  /* check parameters */
  if ( (rcode = _unur_dstd_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* run special init routine for generator */
  GEN->is_inversion = FALSE;   /* reset flag for inversion method */
  if ( DISTR.init(NULL,gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return UNUR_ERR_GEN_DATA;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug & DSTD_DEBUG_REINIT)
      _unur_dstd_debug_chg_pmfparams( gen );
#endif

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_dstd_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
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

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dstd_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DSTD_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;    /* will be set in _unur_dstd_init() */
  gen->destroy = _unur_dstd_free;
  gen->clone = _unur_dstd_clone;
  gen->reinit = _unur_dstd_reinit;

  /* defaults */
  GEN->gen_param = NULL;  /* parameters for the generator      */
  GEN->n_gen_param = 0;   /* (computed in special GEN->init()   */
  GEN->gen_iparam = NULL; /* array for integer parameters      */
  GEN->n_gen_iparam = 0;
  GEN->is_inversion = FALSE;    /* method not based on inversion             */
  GEN->sample_routine_name = NULL ;  /* name of sampling routine  */

  /* copy some parameters into generator object */
  GEN->Umin        = 0;    /* cdf at left boundary of domain   */
  GEN->Umax        = 1;    /* cdf at right boundary of domain  */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dstd_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_dstd_create() */

/*---------------------------------------------------------------------------*/

int
_unur_dstd_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* domain must not be changed */
  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* domain has been modified --> not allowed */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed");
    return UNUR_ERR_GEN_DATA;
  }

  return UNUR_SUCCESS;
} /* end of _unur_dstd_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dstd_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dstd_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSTD_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy parameters for special generators */
  if (GEN->gen_param) {
    CLONE->gen_param = _unur_xmalloc( GEN->n_gen_param * sizeof(double) );
    memcpy( CLONE->gen_param, GEN->gen_param, GEN->n_gen_param * sizeof(double) );
  }
  if (GEN->gen_iparam) {
    CLONE->gen_iparam = _unur_xmalloc( GEN->n_gen_iparam * sizeof(int) );
    memcpy( CLONE->gen_iparam, GEN->gen_iparam, GEN->n_gen_iparam * sizeof(int) );
  }

  return clone;

#undef CLONE
} /* end of _unur_dstd_clone() */

/*---------------------------------------------------------------------------*/

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
  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);

  /* check input */
  if ( gen->method != UNUR_METH_DSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->gen_param)   free(GEN->gen_param);
  if (GEN->gen_iparam)  free(GEN->gen_iparam);

  _unur_generic_free(gen);
} /* end of _unur_dstd_free() */

/*****************************************************************************/

/** 
    double _unur_dstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

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
_unur_dstd_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );

  /* sampling routine */
  fprintf(LOG,"%s: sampling routine = ",gen->genid);
  if (GEN->sample_routine_name)
    fprintf(LOG,"%s()",GEN->sample_routine_name);
  else
    fprintf(LOG,"(Unknown)");
  if (GEN->is_inversion)
    fprintf(LOG,"   (Inversion)");
  fprintf(LOG,"\n%s:\n",gen->genid);

  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    fprintf(LOG,"%s: domain has been changed. U in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);
    fprintf(LOG,"%s:\n",gen->genid);
  }

} /* end of _unur_dstd_info_init() */

/*---------------------------------------------------------------------------*/

void 
_unur_dstd_debug_chg_pmfparams( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) parameters of distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: parameters of distribution changed:\n",gen->genid);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\tparam[%d] = %g\n",gen->genid,i,DISTR.params[i]);

} /* end of _unur_dstd_debug_chg_pmfparams() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dstd_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  int samplesize = 10000;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: DSTD (special generator for Discrete STandarD distribution)\n");
  _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
		      (GEN->is_inversion)?"[implements inversion method]" : "");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#urn] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
			(gen->set & DSTD_SET_VARIANT) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_dstd_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
