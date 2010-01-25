/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mixt.c                                                        *
 *                                                                           *
 *   TYPE:      (continuous) univariate random variate                       *
 *   METHOD:    mixture of distributions (meta method)                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given an array of components and array of probabilities              *
 *      produce a value X consistent with the mixture distribution           *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointers to generators of components                                 *
 *      array of probabilities                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2010 Wolfgang Hoermann and Josef Leydold                  *
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
/* #include <distr/distr_source.h> */
#include <distr/cont.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include <methods/dgt.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "mixt.h"
#include "mixt_struct.h"

/* #ifdef UNUR_ENABLE_INFO */
/* #  include <tests/unuran_tests.h> */
/* #endif */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* Variants: none                                                            *\/ */

/* #define MIXT_VARFLAG_VERIFY   0x002u    /\* run verify mode                    *\/ */
/* #define MIXT_VARFLAG_SQUEEZE  0x004u    /\* use universal squeeze if possible  *\/ */

/* #define MIXT_VARFLAG_MIRROR   0x008u    /\* use mirror principle               *\/ */
/* /\* not implemented yet! *\/ */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* Debugging flags                                                           *\/ */
/* /\*    bit  01    ... pameters and structure of generator (do not use here)   *\/ */
/* /\*    bits 02-12 ... setup                                                   *\/ */
/* /\*    bits 13-24 ... adaptive steps                                          *\/ */
/* /\*    bits 25-32 ... trace sampling                                          *\/ */

/* #define MIXT_DEBUG_REINIT    0x00000010u   /\* print parameters after reinit  *\/ */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* Flags for logging set calls                                               *\/ */

/* #define MIXT_SET_CDFMODE      0x001u    /\* CDF at mode is known               *\/ */
/* #define MIXT_SET_PDFMODE      0x002u    /\* PDF at mode is set                 *\/ */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MIXT"          /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

/* static int _unur_mixt_reinit( struct unur_gen *gen ); */
/* /\*---------------------------------------------------------------------------*\/ */
/* /\* Reinitialize generator.                                                   *\/ */
/* /\*---------------------------------------------------------------------------*\/ */

static struct unur_gen *_unur_mixt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_mixt_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_mixt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_mixt_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_indexgen( const double *prob, int n_prob );
/*---------------------------------------------------------------------------*/
/* create generator for index.                                               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mixt_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/* #ifdef UNUR_ENABLE_INFO */
/* static void _unur_mixt_info( struct unur_gen *gen, int help ); */
/* /\*---------------------------------------------------------------------------*\/ */
/* /\* create info string.                                                       *\/ */
/* /\*---------------------------------------------------------------------------*\/ */
/* #endif */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* abbreviations *\/ */

/* #define DISTR_IN  distr->data.cont      /\* data for distribution object      *\/ */

#define PAR       ((struct unur_mixt_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mixt_gen*)gen->datap) /* data for generator object */
/* #define DISTR     gen->distr->data.cont /\* data for distribution in generator object *\/ */

/* #define BD_LEFT   domain[0]             /\* left boundary of domain of distribution *\/ */
/* #define BD_RIGHT  domain[1]             /\* right boundary of domain of distribution *\/ */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

/* #define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /\* call to PDF         *\/ */

#define INDEX     gen_aux

#define PROB      gen_aux->distr->data.discr.pv
#define COMP      gen_aux_list
#define N_COMP    n_gen_aux_list

/* /\*---------------------------------------------------------------------------*\/ */

/* #define _unur_mixt_getSAMPLE(gen) \ */
/*    ( ((gen)->variant & MIXT_VARFLAG_VERIFY) \ */
/*      ? _unur_mixt_sample_check : _unur_mixt_sample ) */

/* /\*---------------------------------------------------------------------------*\/ */
/* /\* constants                                                                 *\/ */

/* #define SQRT2    1.4142135623731 */

/* /\*---------------------------------------------------------------------------*\/ */

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mixt_new( int n, const double *prob, struct unur_gen **comp )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   n    ... number of components                                      */
     /*   prob ... probabilities for components                              */
     /*   comp ... array of pointers to components (generators)              */
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
  _unur_check_NULL( GENTYPE, prob, NULL );
  _unur_check_NULL( GENTYPE, comp, NULL );
  if (n<1) { _unur_error(GENTYPE,UNUR_ERR_DISTR_DOMAIN,"n < 1"); return NULL; }
  /* checking type of generator objects in 'comp' is delayed to init */

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mixt_par) );
  COOKIE_SET(par,CK_MIXT_PAR);

  /* copy input */
  par->distr    = NULL;      /* pointer to distribution object               */

  /* copy data */
  PAR->n_comp   = n;         /* number of components                         */
  PAR->prob     = prob;      /* probabilities for components                 */
  PAR->comp     = comp;      /* array of pointers to components (generators) */

  /* set default values */
  par->method   = UNUR_METH_MIXT;   /* method and default variant            */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* inidicate default parameters          */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_mixt_init;

  return par;

} /* end of unur_mixt_new() */

/*****************************************************************************/

/* int  */
/* unur_mixt_set_usesqueeze( struct unur_par *par, int usesqueeze ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* set flag for using universal squeeze (default: off)                  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   par    ... pointer to parameter for building generator object      *\/ */
/*      /\*   usesqueeze ... 0 = no squeeze,  !0 = use squeeze                   *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*                                                                      *\/ */
/*      /\* comment:                                                             *\/ */
/*      /\*   no squeeze is the default                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   /\* check arguments *\/ */
/*   _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL ); */
/*   _unur_check_par_object( par, MIXT ); */

/*   /\* we use a bit in variant *\/ */
/*   par->variant = (usesqueeze)  */
/*     ? (par->variant | MIXT_VARFLAG_SQUEEZE)  */
/*     : (par->variant & (~MIXT_VARFLAG_SQUEEZE)); */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */

/* } /\* end of unur_mixt_set_usesqueeze() *\/ */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_mixt_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_MIXT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mixt_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* probabilities */
  gen->INDEX = _unur_mixt_indexgen(PAR->prob,PAR->n_comp);

  /* components */
  gen->N_COMP = PAR->n_comp;    /* number of components                         */
  gen->COMP = _unur_xmalloc( gen->N_COMP * sizeof(struct unur_gen *));
  for (i=0; i<gen->N_COMP; i++)
    gen->COMP[i] = unur_gen_clone(PAR->comp[i]);

  /* free parameters */
  _unur_par_free(par);

  /* check parameters */
  if (_unur_mixt_check_par(gen) != UNUR_SUCCESS) {
    _unur_mixt_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_mixt_debug_init(gen);
#endif

  return gen;
} /* end of _unur_mixt_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mixt_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mixt_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MIXT_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* object 'par' does not contain a distribution object */
  /* so we create one.                                   */
  gen->distr = unur_distr_cont_new();

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mixt_sample;
  gen->destroy = _unur_mixt_free;
  gen->clone = _unur_mixt_clone;
  gen->reinit = NULL;    /* reinit not implemented ! */

  /* initialize parameters: none */

/* #ifdef UNUR_ENABLE_INFO */
/*   /\* set function for creating info string *\/ */
/*   gen->info = _unur_mixt_info; */
/* #endif */

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_mixt_create() */

/*---------------------------------------------------------------------------*/

int
_unur_mixt_check_par( struct unur_gen *gen )
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
  int i;
  int type;

  /* check probabilities */
  if (gen->INDEX == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"invalid probabilities");
    return UNUR_ERR_GEN_DATA;
  }

  /* check generator objects */
  for (i=0; i<gen->N_COMP; i++) {
    /* all generators must sample from univariate distributions */
    if (gen->COMP[i] == NULL) {
      _unur_error(gen->genid,UNUR_ERR_NULL,"component is NULL");
      return UNUR_ERR_NULL;
    }
    _unur_check_NULL( gen->genid, gen->COMP[i], UNUR_ERR_NULL);
    type = gen->COMP[i]->method & UNUR_MASK_TYPE;
    if ( type != UNUR_METH_DISCR && 
	 type != UNUR_METH_CONT  &&
	 type != UNUR_METH_CEMP  ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"component not univariate");
      return UNUR_ERR_GEN_INVALID;
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_mixt_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mixt_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_mixt_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MIXT_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_mixt_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_mixt_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_MIXT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_mixt_free() */

/*****************************************************************************/

double
_unur_mixt_sample( struct unur_gen *gen )
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
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
/*   double U,V,X,xx,y; */
  int J;
  double X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_MIXT_GEN,INFINITY);

  /* sample index */
  J = unur_sample_discr(gen->INDEX);

  /* sample from selected component */
  X = unur_sample_cont(gen->COMP[J]);

  /* FIXME: discrete ? */

  return X;
} /* end of _unur_mixt_sample() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

struct unur_gen *
_unur_mixt_indexgen( const double *prob, int n_prob )
/*---------------------------------------------------------------------------*/
/* create generator for index using method DGT.                              */
/*                                                                           */
/* parameters:                                                               */
/*   prob   ... probability vector                                           */
/*   n_prob ... length of probability vector                                 */
/*                                                                           */
/* return:                                                                   */
/*   pointer to generator object                                             */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{
  struct unur_distr *distr;
  struct unur_par *par;
  struct unur_gen *igen;

  /* create generator */
  distr = unur_distr_discr_new();
  unur_distr_discr_set_pv(distr, prob, n_prob);
  par = unur_dgt_new(distr);
  igen = unur_init(par);

  /* clear working space */
  unur_distr_free(distr);
  
  return igen;

} /* end of _unur_mixt_indexgen() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_mixt_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  struct unur_gen *comp;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = MIXT (MIXTure of distributions -- meta method)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: sampling routine = _unur_mixt_sample\n",gen->genid);
/*   if (gen->variant & MIXT_VARFLAG_VERIFY) */
/*     fprintf(LOG,"_check"); */
/*   /\* else if (gen->variant & MIXT_VARFLAG_MIRROR)     not implemented *\/ */
/*   /\*   fprintf(LOG,"_mirror"); *\/ */
/*   fprintf(LOG,"()\n%s:\n",gen->genid); */
  fprintf(LOG,"%s:\n",gen->genid);

  /* probabilities */
  fprintf(LOG,"%s: probabilities (%d) = \n",gen->genid, gen->N_COMP);
  fprintf(LOG,"%s:   %g",gen->genid, (gen->PROB)[0]);
  for (i=1; i<gen->N_COMP; i++)
    fprintf(LOG,", %g", (gen->PROB)[i]);
  fprintf(LOG,"\n%s:\n",gen->genid);

  /* components */
  fprintf(LOG,"%s: components (%d):\n",gen->genid, gen->N_COMP);
  for (i=0; i<gen->N_COMP; i++) {
    comp = gen->COMP[i];
    fprintf(LOG,"%s:   [%d]: %s\n",gen->genid, i, comp->genid);
    fprintf(LOG,"%s:\t type = ",gen->genid); 
    switch (comp->distr->type) {
    case UNUR_DISTR_CONT:
    case UNUR_DISTR_CEMP:
      fprintf(LOG,"continuous\n");
      break;
    case UNUR_DISTR_DISCR:
      fprintf(LOG,"discrete\n");
      break;
    default:
      fprintf(LOG,"[unknown]\n");
    }
    fprintf(LOG,"%s:\t name = %s\n",gen->genid, comp->distr->name);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  

} /* end of _unur_mixt_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/* /\*---------------------------------------------------------------------------*\/ */
/* #ifdef UNUR_ENABLE_INFO */
/* /\*---------------------------------------------------------------------------*\/ */

/* void */
/* _unur_mixt_info( struct unur_gen *gen, int help ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* create character string that contains information about the          *\/ */
/*      /\* given generator object.                                              *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen  ... pointer to generator object                               *\/ */
/*      /\*   help ... whether to print additional comments                      *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   struct unur_string *info = gen->infostr; */
/*   struct unur_distr *distr = gen->distr; */
/*   int samplesize = 10000; */
/*   double rc, rc_approx; */

/*   /\* generator ID *\/ */
/*   _unur_string_append(info,"generator ID: %s\n\n", gen->genid); */
  
/*   /\* distribution *\/ */
/*   _unur_string_append(info,"distribution:\n"); */
/*   _unur_distr_info_typename(gen); */
/*   _unur_string_append(info,"   functions = PDF\n"); */
/*   _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]); */
/*   _unur_string_append(info,"   mode      = %g   %s\n", DISTR.mode, */
/* 		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : ""); */
/*   _unur_string_append(info,"   area(PDF) = %g\n", DISTR.area); */
/*   if (gen->set & MIXT_SET_CDFMODE) */
/*     _unur_string_append(info,"   F(mode)   = %g\n", GEN->Fmode);  */
/*   else */
/*     _unur_string_append(info,"   F(mode)   = [unknown]\n");  */

/*   if (help) { */
/*     if ( distr->set & UNUR_DISTR_SET_MODE_APPROX )  */
/*       _unur_string_append(info,"\n[ Hint: %s ]\n", */
/* 			  "You may provide the \"mode\""); */
/*   } */
/*   _unur_string_append(info,"\n"); */

/*   /\* method *\/ */
/*   _unur_string_append(info,"method: MIXT (Simple Ratio-Of-Uniforms)\n"); */
/*   if (gen->set & MIXT_SET_CDFMODE) */
/*     _unur_string_append(info,"   use CDF at mode\n"); */
/*   if (gen->variant & MIXT_VARFLAG_SQUEEZE) */
/*     _unur_string_append(info,"   use squeeze\n"); */
/*   _unur_string_append(info,"\n"); */

/*   /\* performance *\/ */
/*   _unur_string_append(info,"performance characteristics:\n"); */
/*   rc = (gen->set & MIXT_SET_CDFMODE) ? 2. : 4.; */
/*   if (_unur_isfinite(DISTR.BD_RIGHT) || _unur_isfinite(DISTR.BD_LEFT)) { */
/*     rc_approx = unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize); */
/*     _unur_string_append(info,"   rejection constant <= %g  [approx. = %.2f]\n", rc,rc_approx); */
/*   } */
/*   else { */
/*     _unur_string_append(info,"   rejection constant = %g\n", rc); */
/*   } */
/*   _unur_string_append(info,"\n"); */

/*   /\* parameters *\/ */
/*   if (help) { */
/*     _unur_string_append(info,"parameters:\n"); */
/*     if (gen->set & MIXT_SET_CDFMODE) */
/*       _unur_string_append(info,"   cdfatmode = %g\n", GEN->Fmode);  */
/*     else */
/*       _unur_string_append(info,"   cdfatmode = [not set]\n");  */

/*     if (gen->variant & MIXT_VARFLAG_SQUEEZE) */
/*       _unur_string_append(info,"   usesqueeze\n"); */

/*     if (gen->variant & MIXT_VARFLAG_VERIFY) */
/*       _unur_string_append(info,"   verify = on\n"); */

/*     _unur_string_append(info,"\n"); */

/*     /\* Not displayed: */
/*        int unur_mixt_set_pdfatmode( UNUR_PAR *parameters, double fmode ); */
/*     *\/ */
/*   } */

/*   /\* Hints *\/ */
/*   if (help) { */
/*     if ( !(gen->set & MIXT_SET_CDFMODE)) */
/*       _unur_string_append(info,"[ Hint: %s ]\n", */
/* 			  "You can set \"cdfatmode\" to reduce the rejection constant."); */
/*     _unur_string_append(info,"\n"); */
/*   } */

/* } /\* end of _unur_mixt_info() *\/ */

/* /\*---------------------------------------------------------------------------*\/ */
/* #endif   /\* end UNUR_ENABLE_INFO *\/ */
/* /\*---------------------------------------------------------------------------*\/ */
