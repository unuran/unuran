/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      unif.c                                                       *
 *                                                                           *
 *   dummy generator, produces uniform random numbers in UNURAN framework    *
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

/*---------------------------------------------------------------------------*/

#define GENTYPE "UNIF"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_unif_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.unif
#define GEN     gen->data.unif
#define SAMPLE  gen->sample.cont

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_unif_new( int start, int skip )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   start ...  starting point for subsequences (0 = first number)      */
     /*   skip  ...  skip for subsequence (1 = full sequence)                */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   for testing only.                                                  */
     /*   negative number for start and skip are treated as 0 and 1,         */
     /*   respectively.                                                      */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_UNIF_PAR);

  /* copy input */
  PAR.start = start;               /* starting point for subsequence         */
  PAR.skip = skip;                 /* skip for subsequence                   */

  /* set default values */
  par->method   = UNUR_METH_UNIF;  /* method and default variant             */
  par->variant  = 0u;              /* default variant                        */
  par->set      = 0u;              /* inidicate default parameters           */    
  par->urng     = unur_get_default_urng(); /* use default urng               */

  par->genid    = _unur_set_genid(GENTYPE);/* set generator id               */
  par->debug    = UNUR_DEBUGFLAG_DEFAULT;  /* set default debugging flags    */

  /* routine for starting generator */
  par->init = unur_unif_init;

  return par;

} /* end of unur_unif_new() */

/*****************************************************************************/

struct unur_gen *
unur_unif_init( struct unur_par *par )
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
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_UNIF ) {
    _unur_error(par->genid,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_UNIF_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_unif_create(par);
  if (!gen) { free(par); return NULL; }

  /* skip to starting point */
  for (i=0; i<GEN.start; i++)
    _unur_call_urng(gen);

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_unif_init() */

/*****************************************************************************/

double
unur_unif_sample( struct unur_gen *gen )
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
  register int i;
  register double u;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_UNIF_GEN,0.);

  /* sample uniform random number */
  u = _unur_call_urng(gen);

  /* skip to next starting point in sequence */
  for (i=1; i<GEN.skip; i++)
    _unur_call_urng(gen);

  return u;

} /* end of unur_unif_sample() */

/*****************************************************************************/

void
unur_unif_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_UNIF ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_UNIF_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of unur_unif_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_unif_create( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_UNIF_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_UNIF_GEN);

  /* copy generator identifier */
  gen->genid = par->genid;

  /* routines for sampling and destroying generator */
  SAMPLE = unur_unif_sample;
  gen->destroy = unur_unif_free;

  /* copy some parameters into generator object */
  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  GEN.start = PAR.start;    /* starting point for subsequence */
  GEN.skip = PAR.skip;      /* skip for subsequence */

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_unif_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

/** None **/

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
