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
 * ..... beschreibung ....                                                   *
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

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define CSTD_SET_VARIANT          0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "CSTD"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont

#define PAR       par->data.cstd
#define GEN       gen->data.cstd
#define DISTR     gen->distr.data.cont
#define SAMPLE    gen->sample.cont

#define CDF(x) ((*(DISTR.cdf))((x),GEN.pdf_param,GEN.n_pdf_param))

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
  CHECK_NULL(distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (distr->id == UNUR_DISTR_GENERIC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"standard distribution");
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

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  
  /* routine for starting generator */
  par->init = unur_cstd_init;

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
  CHECK_NULL(par,0);
  CHECK_NULL(par->distr,0);

  /* check input */
  _unur_check_par_object( CSTD );

  /* store date */
  old_variant = par->variant;
  par->variant = variant;

  /* check variant. run special init routine only in test mode */
  if (par->DISTR_IN.init != NULL && par->DISTR_IN.init(par,NULL) ) {
    par->set |= CSTD_SET_VARIANT;    /* changelog */
    return 1;
  }

  /* variant not valid */
  _unur_warning(GENTYPE,UNUR_ERR_SET,"unknown variant for special generator");
  par->variant = old_variant;
  return 0;

} /* end if unur_cstd_set_variant() */

/*****************************************************************************/

struct unur_gen *
unur_cstd_init( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);
  CHECK_NULL(par->DISTR_IN.init,NULL);

  /* check input */
  if ( par->method != UNUR_METH_CSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }

  /* create a new empty generator object */
  gen = _unur_cstd_create(par);
  if (!gen) { free(par); return NULL; }

  /* run special init routine for generator */
  if (par->DISTR_IN.init != NULL)
    par->DISTR_IN.init(par,gen);

  /* init successful ?? */
  if (SAMPLE == NULL) {
    /* could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_INIT,"unknown variant for special generator");
    free(par); unur_cstd_free(gen); return NULL; 
  }
  
  /* domain valid for special generator ?? */
  if (!(par->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* domain has been modified */
    if ( ! PAR.is_inversion ) { 
      /* this is not the inversion method */
      _unur_error(GENTYPE,UNUR_ERR_INIT,"domain changed for non inversion method");
      free(par); unur_cstd_free(gen); return NULL; 
    }
    else if (par->DISTR_IN.cdf == NULL) {
      _unur_error(GENTYPE,UNUR_ERR_INIT,"domain changed, c.d.f. required");
      free(par); unur_cstd_free(gen); return NULL; 
    }
    /* compute umin and umax */
    GEN.umin = (par->DISTR_IN.domain[0] > -INFINITY) ? CDF(par->DISTR_IN.domain[0]) : 0.;
    GEN.umax = (par->DISTR_IN.domain[1] < INFINITY)  ? CDF(par->DISTR_IN.domain[1]) : 1.;
  }

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_cstd_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of unur_cstd_init() */

/*****************************************************************************/

/** 
    double unur_cstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

/*****************************************************************************/

void
unur_cstd_free( struct unur_gen *gen )
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
    _unur_warning(GENTYPE,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.gen_param);
  free(gen);

} /* end of unur_cstd_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

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
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_CSTD_GEN);

  /* set generator identifier */
  _unur_set_genid(gen,GENTYPE);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;    /* will be set in unur_cstd_init() */
  gen->destroy = unur_cstd_free;

  /* defaults */
  GEN.gen_param = NULL;  /* parameters for the generator     */
  GEN.n_gen_param = 0;   /* (computed in special GEN.init()  */

  /* copy some parameters into generator object */
  GEN.pdf_param   = DISTR.params;
  GEN.n_pdf_param = DISTR.n_params;
  GEN.umin        = 0;    /* cdf at left boundary of domain  */
  GEN.umax        = 1;    /* cdf at right boundary of domain */

  gen->method = par->method;        /* indicates used method */
  gen->variant = par->variant;      /* indicates variant     */

  _unur_copy_urng_pointer(par,gen);  /* pointer to urng into generator object*/
  _unur_copy_debugflag(par,gen);     /* copy debugging flags into generator object */

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_cstd_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

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
  CHECK_NULL(par,/*void*/);
  COOKIE_CHECK(par,CK_CSTD_PAR,/*void*/);

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

/*****************************************************************************/
#endif
