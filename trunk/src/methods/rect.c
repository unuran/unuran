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
 * Description ...                                                           *
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

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define RECT_SET_DOMAIN           0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "RECT"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_rect_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
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
  PAR.domain    = NULL;            /* vertices of hypercube                  */

  par->method   = UNUR_METH_RECT;  /* method and default variant             */
  par->variant  = 0u;              /* default variant                        */
  par->set      = 0u;              /* inidicate default parameters           */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->genid    = _unur_set_genid(GENTYPE);/* set generator id               */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = unur_rect_init;

  return par;

} /* end of unur_rect_new() */

/*---------------------------------------------------------------------------*/

int
unur_rect_set_domain_vec( struct unur_par *par, double **domain )
     /*----------------------------------------------------------------------*/
     /* set coordinates for domain boundary                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   domain ... coordinates of domain boundary                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,RECT );

  /* check new parameter for generator */
  /** TODO **/

  /* store data */
  par->data.rect.domain = domain;

  /* changelog */
  par->set |= RECT_SET_DOMAIN;

  /* o.k. */
  return 1;

} /* end of unur_rect_set_domain_vec() */

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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_RECT ) {
    _unur_error(par->genid,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_RECT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_rect_create(par);
  if (!gen) { free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_RECT_GEN,/*void*/);

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

  /* check input */
  if ( gen->method != UNUR_METH_RECT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_RECT_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.domain);
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_RECT_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_RECT_GEN);

  /* copy generator identifier */
  gen->genid = par->genid;

  /* routines for sampling and destroying generator */
  SAMPLE = unur_rect_sample_vec;
  gen->destroy = unur_rect_free;

  /* set all pointers to NULL */
  GEN.domain = NULL;

  /* copy some parameters into generator object */
  GEN.dim = PAR.dim;
  gen->method = par->method;        /* indicates used method */
  gen->variant = par->variant;      /* indicates variant     */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_rect_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

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
  
  /* check arguments */
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_RECT_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_RECT_GEN,/*void*/);

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

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
