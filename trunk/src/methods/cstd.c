/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cstd.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generators for standard distribution (from CRAND)            *
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

#include <ctype.h>
#include <string.h>

#include <unur_methods.h>
#include <unur_stdgen.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

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

#define DISTR   distr->data.cont
#define PAR     par->data.cstd
#define GEN     gen->data.cstd
#define SAMPLE  gen->sample.cont

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

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_CSTD_PAR);

  /* check arguments */
  CHECK_NULL(distr,NULL);
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  /* check input */
  if (distr->id == UNUR_DISTR_GENERIC) {
    _unur_error(GENTYPE,UNUR_ERR_GENERIC,"standard distribution required");
    return NULL;
  }

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  PAR.variant   = 0x0u;             /* use default generator                 */
  par->method   = UNUR_METH_CSTD;   /* method and default variant            */
  par->set      = 0UL;              /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default urng               */

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  
  /* routine for starting generator */
  par->init = unur_cstd_init;

  return par;

} /* end of unur_cstd_new() */

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

  /* create a new empty generator object */
  gen = _unur_cstd_create(par);
  if (!gen) { free(par); return NULL; }

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

double
unur_cstd_sample( struct unur_gen *gen )
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
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  switch (gen->distr->id) {
  case UNUR_DISTR_EXPONENTIAL:
    return ( -GEN.pdf_param[0] * log(1. - _unur_call_urng(gen)) );
  case UNUR_DISTR_GAMMA:
    return (gammarand( GEN.pdf_param[0],gen->urng ) * GEN.pdf_param[1] + GEN.pdf_param[2]);
  case UNUR_DISTR_NORMAL:
    return (GEN.pdf_param[0] + GEN.pdf_param[1] * nbm(gen->urng));
  default:
    _unur_error(gen->genid,UNUR_ERR_INIT,"unknown distribution.");
    return 0;
  }

} /* end of unur_cstd_sample() */

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
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen->distr);
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

  /* copy distribution object */
  gen->distr = _unur_malloc( sizeof(struct unur_distr) );
  unur_distr_copy( gen->distr, par->distr );

  /* routines for sampling and destroying generator */
  SAMPLE = unur_cstd_sample;
  gen->destroy = unur_cstd_free;

  /* copy some parameters into generator object */
  GEN.pdf_param   = gen->DISTR.params;
  GEN.n_pdf_param = gen->DISTR.n_params;

  GEN.variant = PAR.variant;         /* generator to be used */
  gen->method = par->method;

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

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_debug_cont( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = unur_cstd_sample()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_cstd_info_init() */

#endif

/*****************************************************************************/

