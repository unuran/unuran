/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      rect.h                                                       *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    uniformly distributed in rectangle                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:  dimension                                                    *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Wed Dec  1 20:59:45 CET 1999                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "RECT"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_rect_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_rect_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.rect
#define GEN     gen->data.rect
#define SAMPLE  gen->sample.vec

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_rect_new( int dim )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... dimension of hypercube                                     */
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
  COOKIE_SET(par,CK_RECT_PAR);

  /* copy input */
  PAR.dim     = dim;

  /* set default values */
  PAR.domain  = NULL;            /* vertices of hypercube                    */

  par->method      = UNUR_METH_RECT;  /* method and default variant          */
  par->set         = 0UL;             /* inidicate default parameters        */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  _unur_set_debugflag_default(par);   /* set default debugging flags         */
  _unur_set_genid(par,GENTYPE);       /* set generator identifier            */

  /* routine for starting generator */
  par->init = unur_rect_init;

  return par;

} /* end of unur_rect_new() */

/*****************************************************************************/

struct unur_gen *
unur_rect_init( struct unur_par *par )
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
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_RECT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_rect_create(par);
  if (!gen) { free(par); return NULL; }

#if UNUR_DEBUG & UNUR_DB_INFO
    /* write info into log file */
    if (gen->debug) _unur_rect_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_rect_init() */

/*****************************************************************************/

void
unur_rect_sample_vec( struct unur_gen *gen, double *vec )
     /*---------------------------------------------------------------------------*/
     /* sample from generator                                                     */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen ... pointer to generator object                                     */
     /*   vec ... array for storing random vector                                 */
     /*---------------------------------------------------------------------------*/
{ 
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_RECT_GEN,/*void*/);

  if (GEN.domain == NULL) {
    for (i=0; i<GEN.dim; i++)
      vec[i] = _unur_call_urng(gen);
  }
  else ;
  /* not implemented yet */

} /* end of unur_rect_sample() */

/*****************************************************************************/

void
unur_rect_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_RECT_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  if (GEN.domain) free(GEN.domain);
  free(gen);

} /* end of unur_rect_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_rect_create( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_RECT_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_RECT_GEN);

  /* routines for sampling and destroying generator */
  SAMPLE = unur_rect_sample_vec;
  gen->destroy = unur_rect_free;

  /* set all pointers to NULL */
  GEN.domain = NULL;

  /* copy some parameters into generator object */
  GEN.dim = PAR.dim;
  _unur_copy_urng_pointer(par,gen);  /* pointer to urng into generator object*/
  _unur_copy_debugflag(par,gen);     /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);         /* copy generator identifier            */

  gen->method = par->method;   /* indicates method and variant */

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_rect_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
_unur_rect_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = uniformly distributed in hypercube (RECT)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_rect_sample_vec()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: dimension = %d\n",gen->genid,GEN.dim);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_rect_debug_init() */

#endif

/*****************************************************************************/

