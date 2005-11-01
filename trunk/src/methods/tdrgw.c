/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdrgw.c                                                      *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection - Gilks & Wild variant         *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given logPDF of a log-concave distribution                           *
 *      produce a value x consistent with its density                        *
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
 *                                                                           *
 *   REFERENCES and DESCRIPTION:                                             *
 *                                                                           *
 *   See tdr.c                                                               *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "tdrgw.h"
#include "tdrgw_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TDRGW_VARFLAG_VERIFY     0x0100u   /* flag for verifying mode        */
#define TDRGW_VARFLAG_PEDANTIC   0x0800u   /* whether pedantic checking is used */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TDRGW_DEBUG_IV           0x00000010u
#define TDRGW_DEBUG_SPLIT        0x00010000u
/* #define TDRGW_DEBUG_SAMPLE       0x01000000u */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TDRGW_SET_CPOINTS        0x004u
#define TDRGW_SET_N_CPOINTS      0x008u
#define TDRGW_SET_GUIDEFACTOR    0x010u
#define TDRGW_SET_MAX_SQHRATIO   0x040u
#define TDRGW_SET_MAX_IVS        0x080u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TDRGW"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdrgw_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdrgw_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_tdrgw_sample( struct unur_gen *generator );
static double _unur_tdrgw_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdrgw_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_starting_cpoints( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute intervals from given starting construction points.                */
/*---------------------------------------------------------------------------*/

inline static int _unur_tdrgw_interval_parameter( struct unur_gen *gen, struct unur_tdrgw_interval *iv );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for interval.                                  */
/*---------------------------------------------------------------------------*/

static struct unur_tdrgw_interval *_unur_tdrgw_interval_new( struct unur_gen *gen,
							     double x, double logfx );
/*---------------------------------------------------------------------------*/
/* make a new interval with left construction point x.                       */
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_tangent_intersection_point( struct unur_gen *gen,
						   struct unur_tdrgw_interval *iv, double *ipt );
/*---------------------------------------------------------------------------*/
/* compute cutting point of interval into left and right part.               */
/*---------------------------------------------------------------------------*/

static double _unur_tdrgw_interval_area( struct unur_gen *gen, struct unur_tdrgw_interval *iv,
					 double slope, double x );
static double _unur_tdrgw_interval_logarea( struct unur_gen *gen, struct unur_tdrgw_interval *iv,
					    double slope, double x );
/*---------------------------------------------------------------------------*/
/* compute area below piece of hat or squeeze in interval.                   */
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_interval_split( struct unur_gen *gen,
				       struct unur_tdrgw_interval *iv_old, double x, double logfx );
/*---------------------------------------------------------------------------*/
/* split am interval point x. return 0 if not successful.                    */                                           
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_improve_hat( struct unur_gen *gen, struct unur_tdrgw_interval *iv,
				    double x, double logfx);
/*---------------------------------------------------------------------------*/
/* improve hat function by splitting interval                                */
/*---------------------------------------------------------------------------*/

static int _unur_tdrgw_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_debug_init_start( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after (almost empty generator) object has been created.             */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_tdrgw_debug_split_start( const struct unur_gen *gen,
					    const struct unur_tdrgw_interval *iv,
					    double x, double logfx );
static void _unur_tdrgw_debug_split_stop( const struct unur_gen *gen,
					   const struct unur_tdrgw_interval *iv_left,
					   const struct unur_tdrgw_interval *iv_right );
/*---------------------------------------------------------------------------*/
/* print before and after an interval has been split (not / successfully).   */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_tdrgw_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_tdrgw_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define logPDF(x)  _unur_cont_logPDF((x),(gen->distr))   /* call to logPDF   */
#define dlogPDF(x) _unur_cont_dlogPDF((x),(gen->distr))  /* call to derivative of log PDF */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_tdrgw_new( const struct unur_distr* distr )
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

  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF"); return NULL; }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of logPDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_tdrgw_par) );
  COOKIE_SET(par,CK_TDRGW_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR->guide_factor        = 2.;     /* size of guide table / number of intervals */

  PAR->starting_cpoints    = NULL;   /* pointer to array of starting points  */
  PAR->n_starting_cpoints  = 2;      /* number of starting points            */
  PAR->max_ivs             = 100;    /* maximum number of intervals          */
  PAR->max_ratio           = 0.99;   /* bound for ratio  Atotal / Asqueeze   */
  PAR->bound_for_adding    = 0.5;    /* do not add a new construction point in an interval,
				       where ambigous region is too small, i.e. if
				       area / ((A_hat - A_squeeze)/number of segments) < bound_for_adding */
 
  par->method   = UNUR_METH_TDRGW;   /* method                               */
  par->variant  = 0u;                /* default variant                      */

  par->set      = 0u;               /* inidicate default parameters          */
  par->urng     = unur_get_default_urng(); /* use default URNG               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_tdrgw_init;

  return par;

} /* end of unur_tdrgw_new() */

/*****************************************************************************/

int
unur_tdrgw_set_max_sqhratio( struct unur_par *par, double max_ratio )
     /*----------------------------------------------------------------------*/
     /* set bound for ratio A(squeeze) / A(hat)                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ratio ... upper bound for ratio to add a new construction point*/
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1.+DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ratio = max_ratio;

  /* changelog */
  par->set |= TDRGW_SET_MAX_SQHRATIO;

  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tdrgw_get_sqhratio( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio    ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDRGW, INFINITY );

  return (GEN->Asqueeze / GEN->Atotal);

} /* end of unur_tdrgw_get_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tdrgw_get_hatarea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below hat                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDRGW, INFINITY );

  return GEN->Atotal;

} /* end of unur_tdrgw_get_hatarea() */

/*---------------------------------------------------------------------------*/

double
unur_tdrgw_get_squeezearea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below squeeze                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TDRGW, INFINITY );

  return GEN->Asqueeze;

} /* end of unur_tdrgw_get_squeezearea() */

/*---------------------------------------------------------------------------*/

int
unur_tdrgw_set_max_intervals( struct unur_par *par, int max_ivs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ivs   ... maximum number of intervals                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= TDRGW_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int 
unur_tdrgw_set_cpoints( struct unur_par *par, int n_cpoints, const double *cpoints )
     /*----------------------------------------------------------------------*/
     /* set construction points for hat function                             */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   n_cpoints ... number of starting construction points               */
     /*   cpoints   ... pointer to array of starting construction points     */
     /*                 (NULL for changing only the number of default points)*/
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_cpoints < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }

  if (cpoints)
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_cpoints; i++ )
      if (cpoints[i] <= cpoints[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }

  /* store date */
  PAR->starting_cpoints = cpoints;
  PAR->n_starting_cpoints = n_cpoints;

  /* changelog */
  par->set |= TDRGW_SET_N_CPOINTS | ((cpoints) ? TDRGW_SET_CPOINTS : 0);

  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tdrgw_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= TDRGW_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_tdrgw_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TDRGW_VARFLAG_VERIFY) : (par->variant & (~TDRGW_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tdrgw_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TDRGW, UNUR_ERR_GEN_INVALID );

  /* we use a bit in variant */
  gen->variant = (verify) ? (gen->variant | TDRGW_VARFLAG_VERIFY) : (gen->variant & (~TDRGW_VARFLAG_VERIFY));

  /* sampling routines */
  SAMPLE = (gen->variant & TDRGW_VARFLAG_VERIFY) ? _unur_tdrgw_sample_check : _unur_tdrgw_sample;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdrgw_chg_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tdrgw_set_pedantic( struct unur_par *par, int pedantic )
     /*----------------------------------------------------------------------*/
     /* turn pedantic mode on/off                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   pedantic ... 0 = no pedantic mode, !0 = use pedantic mode          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   pedantic is the default                                            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TDRGW );

  /* we use a bit in variant */
  par->variant = (pedantic) ? (par->variant | TDRGW_VARFLAG_PEDANTIC) : (par->variant & (~TDRGW_VARFLAG_PEDANTIC));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tdrgw_set_pedantic() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_tdrgw_init( struct unur_par *par )
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
/*   int i,k; */

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_TDRGW ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TDRGW_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tdrgw_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdrgw_debug_init_start(par,gen);
#endif

  /* get starting points */
  if (_unur_tdrgw_starting_cpoints(par,gen)!=UNUR_SUCCESS) {
    _unur_par_free(par); _unur_tdrgw_free(gen);
    return NULL;
  }

  /* compute intervals for given starting points */
  if (_unur_tdrgw_starting_intervals(par,gen)!=UNUR_SUCCESS) {
    _unur_par_free(par); _unur_tdrgw_free(gen);
    return NULL;
  }

  /* update maximal number of intervals */
  if (GEN->n_ivs > GEN->max_ivs) {
    GEN->max_ivs = GEN->n_ivs;
  }

  /* make initial guide table */
  _unur_tdrgw_make_guide_table(gen);

  /* free parameters */
  _unur_par_free(par);

  /* is there any hat at all ? */
  if (GEN->Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdrgw_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdrgw_debug_init_finished(gen);
#endif
  
  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

  /* o.k. */
  return gen;

} /* end of _unur_tdrgw_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdrgw_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TDRGW_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_tdrgw_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDRGW_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  gen->destroy = _unur_tdrgw_free;
  gen->clone = _unur_tdrgw_clone;

  /* sampling routines */
  SAMPLE = (par->variant & TDRGW_VARFLAG_VERIFY) ? _unur_tdrgw_sample_check : _unur_tdrgw_sample;
  
  /* set all pointers to NULL */
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN->max_ivs = max(2*PAR->n_starting_cpoints,PAR->max_ivs);  /* maximum number of intervals */
  GEN->max_ratio = PAR->max_ratio;    /* bound for ratio  Atotal / Asqueeze    */

  /* copy variant */
  gen->variant = par->variant;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_tdrgw_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdrgw_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_tdrgw_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_tdrgw_interval *iv, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = iv->next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_tdrgw_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tdrgw_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
    }
    /* next step */
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;


  /* make new guide table */
  CLONE->guide = NULL;
  _unur_tdrgw_make_guide_table(clone);

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tdrgw_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TDRGW ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdrgw_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_tdrgw_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free table */
  if (GEN->guide)  free(GEN->guide);

  /* free other memory not stored in list */
  _unur_generic_free(gen);

} /* end of _unur_tdrgw_free() */

/*---------------------------------------------------------------------------*/

double
_unur_tdrgw_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (original variant by Gilks & Wild)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*   squeeze(x) = f(x0) * exp(sq * (x-x0))                              */
     /*                                                                      */
     /*   left hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                   */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = f(x1) * exp( (Tf)'(x1) *  (x-x1) )                  */
     /*   generation:                                                        */
     /*      X = x1 + 1/(Tf)'(x1) * \log( (Tf)'(x1)/f(x1) * U + 1 )          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_tdrgw_interval *iv, *pt;
  double U, logV;                   /* (log of) uniform random number        */
  double X;                         /* generated point                       */
  double logfx, logsqx, loghx;      /* log of density, squeeze, and hat at X */
  double x0, logfx0, dlogfx0, fx0;  /* construction point and logPDF at x0   */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDRGW_GEN,INFINITY);

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(gen->urng);

    /* look up in guide table and search for segment */
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* l.h.s. or r.h.s. of hat */
    if (-U < iv->Ahatr) { /* right */
      pt = iv->next;
      /* U unchanged */
    }
    else {                /* left */
      pt = iv;
      U += iv->Ahat;
    }

    /* PDF at x0 */
    x0 = pt->x;
    logfx0 = pt->logfx;
    dlogfx0 = pt->dlogfx;
    fx0 = exp(logfx0);

    /* random variate */
    if (dlogfx0 == 0.)
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      /* x = x0 + log(t + 1.) / dlogfx0; is cheaper but numerical unstable */
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    
    /* log of hat at x */
    loghx = logfx0 + dlogfx0*(X - x0);

    /* log of a random point between 0 and hat at x */
    logV = log(_unur_call_urng(gen->urng)) + loghx;
    
    /* log of spueeze at x */
    logsqx = (iv->Asqueeze > 0.) ? iv->logfx + iv->sq*(X - iv->x) : -INFINITY;
 
    /* below squeeze ? */
    if (logV <= logsqx)
      return X;
    
    /* log of PDF at x */
    logfx = logPDF(X);

    /* below PDF ? */
    if (logV <= logfx)
      return X;

    /* being above PDF is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdrgw_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & TDRGW_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* else reject and try again */

  }

} /* end of _unur_tdrgw_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_tdrgw_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results                             */
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
  struct unur_tdrgw_interval *iv, *pt;
  double U, logV;                   /* (log of) uniform random number        */
  double X;                         /* generated point                       */
  double logfx, logsqx, loghx;      /* log of density, squeeze, and hat at X */
  double x0, logfx0, dlogfx0, fx0;  /* construction point and logPDF at x0   */

  int error = 0;                    /* indicates error                       */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDRGW_GEN,INFINITY);

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(gen->urng);

    /* evaluate inverse of hat CDF */

    /* look up in guide table and search for segment */
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* l.h.s. or r.h.s. of hat */
    if (-U < iv->Ahatr) { /* right */
      pt = iv->next;
      /* U unchanged */
    }
    else {                /* left */
      pt = iv;
      U += iv->Ahat;
    }

    /* PDF at x0 */
    x0 = pt->x;
    logfx0 = pt->logfx;
    dlogfx0 = pt->dlogfx;
    fx0 = exp(logfx0);

    /* random variate */
    if (dlogfx0 == 0.)
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      /* x = x0 + log(t + 1.) / dlogfx0; is cheaper but numerical unstable */
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    
    /* log of hat at x */
    loghx = logfx0 + dlogfx0*(X - x0);

    /* log of spueeze at x */
    logsqx = (iv->Asqueeze > 0.) ? iv->logfx + iv->sq*(X - iv->x) : -INFINITY;
 
    /* log of PDF at x */
    logfx = logPDF(X);

    /* check result */
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(logfx,loghx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not log-concave!");
      error = 1;
    }
    if (_unur_FP_less(logfx,logsqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not log-concave!");
      error = 1;
    }

    /* log of a random point between 0 and hat at x */
    logV = log(_unur_call_urng(gen->urng)) + loghx;

    /* below squeeze ? */
    if (logV <= logsqx)
      return X;

    /* below PDF ? */
    if (logV <= logfx)
      return X;

    /* being above PDF is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdrgw_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & TDRGW_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* reject and try again */

  }

} /* end of _unur_tdrgw_sample_check() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_improve_hat( struct unur_gen *gen, struct unur_tdrgw_interval *iv,
			  double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* improve hat function by splitting interval                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   iv         ... pointer to interval that has to be split            */
     /*   x          ... splitting point                                     */
     /*   logfx      ... value of logPDF at splitting point                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... improving hat successful                       */
     /*   others          ... error: PDF not monotone in interval            */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* is there any reason to improve hat ? */
  if (! (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) ) {
    /* no more construction points (avoid calling this function any more) */
    GEN->max_ivs = GEN->n_ivs;
    return UNUR_SUCCESS;
  }

  /* add construction point */
  result = _unur_tdrgw_interval_split(gen, iv, x, logfx);
  if (result!=UNUR_SUCCESS && result!=UNUR_ERR_SILENT) {
    /* condition for PDF is violated! */
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    if (gen->variant & TDRGW_VARFLAG_PEDANTIC) {
      /* replace sampling routine by dummy routine that just returns INFINITY */
      SAMPLE = _unur_sample_cont_error;
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /* splitting successful --> update guide table */
  /** TODO: it is not necessary to update the guide table every time. 
      But then (1) some additional bookkeeping is required and
      (2) the guide table method requires a acc./rej. step. **/
  _unur_tdrgw_make_guide_table(gen);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_tdrgw_improve_hat() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting intervals.                  */
     /* if not provided as arguments compute these                           */
     /* by means of the "equiangular rule" from AROU.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdrgw_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, logfx, logfx_last;
  int is_increasing;
  int i;
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDRGW_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);
  
  /* reset counter of intervals */
  GEN->n_ivs = 0;

  /* prepare for computing construction points */
  if (!PAR->starting_cpoints) {
    /* angles of boundary of domain */
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = DISTR.BD_LEFT;
  is_increasing = TRUE;
    
  logfx = logfx_last = _unur_isfinite(x) ? logPDF(x) : -INFINITY;
  iv = GEN->iv = _unur_tdrgw_interval_new( gen, x, logfx );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  /* logPDF(x) overflow */

  /* now all the other points */
  for( i=0; i<=PAR->n_starting_cpoints; i++ ) {

    /* construction point */
    if (i < PAR->n_starting_cpoints) {
      if (PAR->starting_cpoints) {
	/* construction points provided by user */
	x = PAR->starting_cpoints[i];
	/* check starting point */
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle );
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of PDF at starting point */
    logfx = _unur_isfinite(x) ? logPDF(x) : -INFINITY;

    /* check value of PDF at starting point */
    if (!is_increasing && logfx > logfx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* check whether we are outside of support of PDF */
    if (!_unur_isfinite(logfx) && !_unur_isfinite(logfx_last) ) {
      /* we do not need two such point */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<PAR->n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the PDF is constant 0 on all construction points.
	     then we need both boundary points. */
	  iv->x = x;  /* we only have to change x, everything else remains unchanged */
	  x_last = x;
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with logPDF(x) > -INFINITY */
	break;
    }
    
    /* need a new interval */
    iv->next = _unur_tdrgw_interval_new( gen, x, logfx );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;    /* logPDF(x) overflow */
    
    /* skip pointer to current interval */
    iv = iv->next;

    /* PDF still increasing ? */
    if (is_increasing && logfx < logfx_last)
      is_increasing = FALSE;

    /* store last computed values */
    x_last = x;
    logfx_last = logfx;

  }

  /* we have left the loop with the right boundary of the support of PDF
     make shure that we will never use iv for sampling. */
  iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
  iv->Acum = INFINITY;
  iv->ip = iv->x;
  iv->next = NULL;         /* terminate list */
  --(GEN->n_ivs);          /* we do not count the terminating interval */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdrgw_starting_cpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdrgw_interval *iv, *iv_new, *iv_tmp;
  double x, logfx;              /* construction point, value of logPDF at x  */
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDRGW_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_TDRGW_IV,UNUR_ERR_COOKIE);
  
  /* compute paramters for all intervals */
  for( iv=GEN->iv; iv->next != NULL; ) {

    /* compute parameters for interval */
    switch (_unur_tdrgw_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:   /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN->n_ivs);
      
      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make sure that we will never use this segment */
	iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      continue;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area below hat infinite.
       insert new construction point. */
    x = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */

    /* value of logPDF at x */
    logfx = logPDF(x);

    /* add a new interval, but check if we had to used too many intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_tdrgw_interval_new( gen, x, logfx );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* logPDF(x) overflow */


    /* if fx is 0, then we can cut off the tail of the distribution
       (since it must be T-concave)  */
    if (!_unur_isfinite(logfx) ) {

      if (!_unur_isfinite(iv->logfx) ) {
	/* cut off left tail */
	iv_new->next = iv->next;
	free(iv);
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	iv = iv_new;
      }

      else if (!_unur_isfinite(iv->next->logfx) ) {
	/* cut off right tail */
	free(iv->next);
	--(GEN->n_ivs);
	iv->next = iv_new;
      }

      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    else {
      /* insert new interval into linked list */
      iv_new->next = iv->next;
      iv->next = iv_new;
    }

  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdrgw_starting_intervals() */

/*---------------------------------------------------------------------------*/

struct unur_tdrgw_interval *
_unur_tdrgw_interval_new( struct unur_gen *gen, double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x       ... left point of new interval                             */
     /*   logfx   ... value of logPDF at x                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdrgw_interval *iv;
/*   double dfx; */

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,NULL);

  /* first check logfx */
  if (!(logfx < INFINITY)) {
    /* overflow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"logPDF(x) overflow");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_xmalloc( sizeof(struct unur_tdrgw_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN->n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_TDRGW_IV);

  /* avoid uninitialized variables */
  iv->Acum = iv->Ahat = iv->Ahatr = iv->Asqueeze = 0.;
  iv->ip = iv->sq = 0.;

  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->logfx = logfx;      /* value of logPDF at x */

  /* derivative of transformed density */
  iv->dlogfx = _unur_isfinite(logfx) ? dlogPDF(x) : INFINITY;
  
  /* the program requires dlogPDF > -INFINITY */
  if ( !(iv->dlogfx > -INFINITY))
    iv->dlogfx = INFINITY;

  return iv;

} /* end of _unur_tdrgw_interval_new() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_interval_parameter( struct unur_gen *gen, struct unur_tdrgw_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (Gilks & Wild variant)                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... do not add this construction point             */
     /*   UNUR_ERR_INF    ... area = INFINITY                                */
     /*   others          ... error (PDF not T-concave)                      */
     /*----------------------------------------------------------------------*/
{
  double Ahatl;    /* area below hat at left side of intersection point */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDRGW_IV,UNUR_ERR_COOKIE); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_TDRGW_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents.
     used to partition interval into left hand part (construction point of tangent
     on the left hand boundary) and right hand part (construction point of tangent
     on the left hand boundary). */
  if ( _unur_tdrgw_tangent_intersection_point(gen,iv,&(iv->ip))!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;

  /* squeeze and area below squeeze */
  if (_unur_isfinite(iv->logfx) && _unur_isfinite(iv->next->dlogfx) ) {

    /* we do not compute the slope when the construction points
       are too close. at least 8 significant digits should remain. */
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   /* construction points too close */

    /* slope of transformed squeeze */
    iv->sq = (iv->next->logfx - iv->logfx) / (iv->next->x - iv->x);

    /* check squeeze */
    /* we have to take care about round off error.
       the following accepts PDFs with might be a little bit not T_concave */
    if ( ( (iv->sq > iv->dlogfx       && !_unur_FP_approx(iv->sq,iv->dlogfx)) ||
	   (iv->sq < iv->next->dlogfx && !_unur_FP_approx(iv->sq,iv->next->dlogfx)) )
	 && iv->next->dlogfx < INFINITY ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* volume below squeeze */
    /* always integrate from point with greater value of transformed density
       to the other point */
    iv->Asqueeze = (iv->logfx > iv->next->logfx)
      ? _unur_tdrgw_interval_area( gen, iv, iv->sq, iv->next->x)
      : _unur_tdrgw_interval_area( gen, iv->next, iv->sq, iv->x);
    
    /* check for fatal numerical errors */
    if (!_unur_isfinite(iv->Asqueeze))
      iv->Asqueeze = 0.;
  }

  else {  /* no squeeze */
    iv->sq = 0.;
    iv->Asqueeze = 0.;
  }

  /* volume below hat */
  Ahatl = _unur_tdrgw_interval_area( gen, iv, iv->dlogfx, iv->ip);
  iv->Ahatr = _unur_tdrgw_interval_area( gen, iv->next, iv->next->dlogfx, iv->ip);

  /* areas below head unbounded ? */
  if (! (_unur_isfinite(Ahatl) && _unur_isfinite(iv->Ahatr)) )
    return UNUR_ERR_INF;

  /* total area */
  iv->Ahat = iv->Ahatr + Ahatl;

  /* check area */
  /* we cannot be more accurate than in the `check squeeze' section */
  if ( iv->Asqueeze > iv->Ahat && !_unur_FP_approx(iv->Asqueeze, iv->Ahat) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"A(squeeze) > A(hat). PDF not T-concave!");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdrgw_interval_parameter() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_interval_split( struct unur_gen *gen, struct unur_tdrgw_interval *iv_oldl, double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* split interval iv_oldl into two intervals at point x                 */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv_oldl ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   logfx   ... value of logPDF at x                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdrgw_interval *iv_newr;  /* pointer to new interval */
  struct unur_tdrgw_interval iv_bak;    /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_TDRGW_IV,UNUR_ERR_COOKIE);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDRGW_DEBUG_SPLIT)
    _unur_tdrgw_debug_split_start( gen,iv_oldl,x,logfx );
#endif

  /* the splitting point must be inside the interval */
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* back up data */
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_tdrgw_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (!_unur_isfinite(logfx)) {
    
    /* one of the two boundary points must be 0, too! */
    if (!_unur_isfinite(iv_oldl->logfx)) {
      /* chop off left part (it's out of support) */
      iv_oldl->x = x;
    }
    else if (!_unur_isfinite(iv_oldl->next->logfx)) {
      /* chop off right part (it's out of support) */
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* compute parameters for chopped interval */
    success = _unur_tdrgw_interval_parameter(gen, iv_oldl);
    
    /* we did not add a new interval */
    iv_newr = NULL;
  }

  else {
    
    /* we need a new interval */
    iv_newr = _unur_tdrgw_interval_new( gen, x, logfx );
    if (iv_newr == NULL) {
      /* logPDF(x) overflow */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* insert into linked list */
    iv_newr->next = iv_oldl->next;
    iv_oldl->next = iv_newr;
    
    /* compute parameters for interval */
    success   = _unur_tdrgw_interval_parameter(gen, iv_oldl);
    success_r = _unur_tdrgw_interval_parameter(gen, iv_newr);
    
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }

  /* successfull ? */
  if (success!=UNUR_SUCCESS) {
    /* cannot split interval at given point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");

    /* the case of unbounded hat is treated as round-off error for
       very steep tangents. so we simply do not add this construction point. */

    /* restore old interval */
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_tdrgw_interval));

    /* decrement counter for intervals and free unused interval */
    if (iv_newr) {
      --(GEN->n_ivs);
      free( iv_newr );
    }

  return ( (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
	   ? UNUR_ERR_GEN_CONDITION : UNUR_SUCCESS );
  }

  /* successful */

  /* update total area below hat and squeeze */
  GEN->Atotal   = ( GEN->Atotal - iv_bak.Ahat
		   + iv_oldl->Ahat + ((iv_newr) ? iv_newr->Ahat : 0.) );
  GEN->Asqueeze = ( GEN->Asqueeze - iv_bak.Asqueeze
		   + iv_oldl->Asqueeze + ((iv_newr) ? iv_newr->Asqueeze : 0. ) );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDRGW_DEBUG_SPLIT)
    _unur_tdrgw_debug_split_stop( gen,iv_oldl,iv_newr );
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdrgw_interval_split() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_tangent_intersection_point( struct unur_gen *gen, struct unur_tdrgw_interval *iv, double *ipt )
     /*----------------------------------------------------------------------*/
     /* compute cutting point of interval into left and right part.          */
     /* (1) use intersection point of tangents of transformed hat.           */
     /* (2) use mean point if (1) is unstable due to roundoff errors.        */
     /* (3) use boundary point which is closer to the mode. this is          */
     /*     important when the transformed tagents are extremely steep.      */
     /*     (This might cause a serious roundoff error while sampling.)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   ipt ... pointer to intersection point                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDRGW_IV,UNUR_ERR_COOKIE); 

  /*
     case: there is no tangent at one of the boundary points of the interval
           (then the slope is INFINITY)
     or
     case: the tangents are too steep  (--> case (3))
  */
  if ( iv->dlogfx > 1.e+140 ) {
    *ipt = iv->x;        /* intersection point = left boundary of interval */
    return UNUR_SUCCESS;
  }
  if ( iv->next->dlogfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dlogfx)) {
    *ipt = iv->next->x;   /* intersection point = right boundary of interval */
    return UNUR_SUCCESS;
  }
  /** TODO: 1.e+140 (= sqrt(DBL_MAX) / 1.e15) is arbitrary  **/

  /* test for T-concavity */
  if ( _unur_FP_less( iv->dlogfx, iv->next->dlogfx ) ) {

    /* it might happen because of round-off errors
       that iv->next->dTfx is almost zero although it should be large.
       thus we ignore this case. */
    if ( fabs(iv->dlogfx) < DBL_EPSILON * fabs(iv->next->dlogfx) ) {
      *ipt = iv->x;        /* intersection point = left boundary of interval */
      iv->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else if ( fabs(iv->next->dlogfx) < DBL_EPSILON * fabs(iv->dlogfx) ) {
      *ipt = iv->next->x;   /* intersection point = right boundary of interval */
      iv->next->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not log-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /** TODO: the following test is too sensitve to roundoff errors **/
  /*    if (iv->next->logfx > iv->x + iv->dlogfx*(iv->next->x - iv->x)) { */
  /*      _unur_warning(gen->genid,UNUR_ERR_INIT,"tangent below PDF, not log-concave!"); */
  /*      return UNUR_ERR_INIT; */
  /*    } */
  
  /* case (2): computing intersection of tangents is unstable */
  if (_unur_FP_approx(iv->dlogfx, iv->next->dlogfx)) {
    /* use mean point */
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }

  /* case (1): compute intersection point of tangents (regular case) */
  *ipt = ( (iv->next->logfx - iv->logfx - iv->next->dlogfx * iv->next->x + iv->dlogfx * iv->x) /
	   (iv->dlogfx - iv->next->dlogfx) );

  /* check position of intersection point */
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    /* intersection point of tangents not in interval.
       This is mostly the case for numerical reasons.
       Thus we is the center of the interval instead.
       if the PDF not T-concave, it will catched at a later
       point when we compare slope of tangents and squeeze. */
    *ipt = 0.5 * (iv->x + iv->next->x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdrgw_tangent_intersection_point() */

/*---------------------------------------------------------------------------*/

double
_unur_tdrgw_interval_area( struct unur_gen *gen, struct unur_tdrgw_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute area below piece of hat or squeeze in                             */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   area                                                                    */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(t-x0)) dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDRGW_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDRGW_IV,INFINITY); 

  return exp(_unur_tdrgw_interval_logarea(gen, iv, slope, x ));

} /* end of _unur_tdrgw_interval_area() */

/*---------------------------------------------------------------------------*/

double
_unur_tdrgw_interval_logarea( struct unur_gen *gen, struct unur_tdrgw_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute log of area below piece of hat or squeeze in                      */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   log of area                                                             */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(t-x0)) dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double x0, logfx0;
  double logxdiff;
  double t, logt;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDRGW_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDRGW_IV,INFINITY); 

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return -INFINITY;

  /* if the construction point is at infinity, we cannot compute an area.
     (in this case we should have x == iv->x == INFINITY). */
  if (!_unur_isfinite(iv->x)) 
    return INFINITY;

  /* unbounded? */
  if ( !_unur_isfinite(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  /* construction point x0 of tangent and log of PDF at x0 */
  x0 = iv->x;
  logfx0 = iv->logfx;

  /* log of |x - x=| */
  logxdiff = log(fabs(x - x0));

  /* case: hat/squeeze constant --> area = f(x0) * |x - x0| */
  if (slope == 0.)
    return (_unur_isfinite(x) ? logfx0 + logxdiff : INFINITY);

  /* case: domain unbounded --> area = f(x0) / |slope| */
  if (!_unur_isfinite(x))
    return (logfx0 - log(fabs(slope)));

  /* case bounded domain --> area = | f(x0)/slope * (\exp(slope*(x-x0))-1) | */
  /* have to deal with numerical problems when x \approx x0                  */
  t = slope * (x - x0);
  logt = log(fabs(slope)) + logxdiff;

  if (fabs(t) > 1.e-6) {
    if (t > MAXLOG / 10.)
      return ( logfx0 + logxdiff + t - logt );
    else
      return ( logfx0 + logxdiff + log( fabs(exp(t) - 1.) ) - log(fabs(t)) );
  }

  else  /* use Taylor series */
    return (logfx0 + logxdiff + log1p(t/2. + t*t/6.));

} /* end of _unur_tdrgw_interval_logarea() */

/*---------------------------------------------------------------------------*/

int
_unur_tdrgw_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdrgw_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDRGW_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    int max_guide_size = (GEN->guide_factor > 0.) ? (GEN->max_ivs * GEN->guide_factor) : 1;
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tdrgw_interval*) );
  }

  /* first we need cumulated areas in intervals */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDRGW_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN->guide_size = (int)(GEN->n_ivs * GEN->guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDRGW_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   /* this is the last virtual intervall --> do not use */
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;

  return UNUR_SUCCESS;
} /* end of _unur_tdrgw_make_guide_table() */

/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_debug_init_start( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print after (almost empty generator) object has been created.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TDRGW_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = Gilks & Wild variant of transformed density rejection\n",gen->genid);
  fprintf(log,"%s: transformation T_c(x) = log(x)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tdrgw_sample",gen->genid);
  if (par->variant & TDRGW_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(par,TDRGW_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,PAR->max_ratio*100.);
  _unur_print_if_default(par,TDRGW_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s: Derandomized ARS disabled ",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR->guide_factor);
  _unur_print_if_default(par,TDRGW_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(par,TDRGW_SET_N_CPOINTS);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & TDRGW_SET_CPOINTS)
    for (i=0; i<PAR->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",PAR->starting_cpoints[i]);
    }
  else
    fprintf(log," use \"equdistribution\" rule [default]");
  fprintf(log,"\n%s:\n",gen->genid);
  
  fflush(log);

} /* end of _unur_tdrgw_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_tdrgw_debug_intervals(gen,"INIT completed",TRUE);

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdrgw_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void 
_unur_tdrgw_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile (orig. variant by Gilks & Wild) */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen         ... pointer to generator object                        */
     /*   header      ... header for table                                   */
     /*   print_areas ... whether table of areas should be printed           */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tdrgw_interval *iv;
  double sAsqueeze, sAhatl, sAhatr, Atotal;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (header) fprintf(log,"%s:%s\n",gen->genid,header);

  fprintf(log,"%s:Intervals: %d\n",gen->genid,GEN->n_ivs);
  if (GEN->iv) {
    if (gen->debug & TDRGW_DEBUG_IV) {
      fprintf(log,"%s: Nr.            tp            ip       logf(tp)     dlogf(tp)       squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDRGW_IV,RETURN_VOID); 
	fprintf(log,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->ip, iv->logfx, iv->dlogfx, iv->sq);
      }
      COOKIE_CHECK(iv,CK_TDRGW_IV,RETURN_VOID); 
      fprintf(log,"%s:[...]: %#12.6g                %#12.6g  %#12.6g\n", gen->genid,
	      iv->x, iv->logfx, iv->dlogfx);
    }
    fprintf(log,"%s:\n",gen->genid);
  }
  else
    fprintf(log,"%s: No intervals !\n",gen->genid);

  if (!print_areas || GEN->Atotal <= 0.) return;

  /* print and sum areas below squeeze and hat */
  Atotal = GEN->Atotal;
  if (gen->debug & TDRGW_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat (left and right)\t\t   cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
    if (GEN->iv) {
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDRGW_IV,RETURN_VOID); 
	sAsqueeze += iv->Asqueeze;
	sAhatl += iv->Ahat - iv->Ahatr;
	sAhatr += iv->Ahatr;
	fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
		iv->Ahat-iv->Ahatr, iv->Ahatr, iv->Ahat * 100. / Atotal, 
		iv->Acum, iv->Acum * 100. / Atotal);
      }
      fprintf(log,"%s:       ----------  ---------  |  ------------------------  ---------  +\n",gen->genid);
      fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)            %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAsqueeze, sAsqueeze * 100. / Atotal,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / Atotal);
      fprintf(log,"%s:\n",gen->genid);
    }
  }

  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_tdrgw_debug_intervals() */

/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    _unur_tdrgw_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(log,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
    _unur_tdrgw_debug_intervals(gen,"Intervals after failure:",FALSE);
  }
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdrgw_debug_free() */

/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_debug_split_start( const struct unur_gen *gen, 
				const struct unur_tdrgw_interval *iv,
				double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting interval                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   iv    ... pointer to interval                                      */
     /*   x     ... split at this point                                      */
     /*   logfx ... value of logPDF at x                                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_TDRGW_IV,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: split interval at x = %g \t\tlogf(x) = %g\n",gen->genid,x,logfx);
  fprintf(log,"%s: old interval:\n",gen->genid);
  fprintf(log,"%s:   left  construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->x,iv->logfx);
  fprintf(log,"%s:   right construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->logfx);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv->Asqueeze,iv->Asqueeze*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv->Ahat - iv->Asqueeze),(iv->Ahat - iv->Asqueeze)*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv->Ahat - iv->Ahatr, iv->Ahatr, iv->Ahat*100./GEN->Atotal);

  fflush(log);

} /* end of _unur_tdrgw_debug_split_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdrgw_debug_split_stop( const struct unur_gen *gen, 
			       const struct unur_tdrgw_interval *iv_left, 
			       const struct unur_tdrgw_interval *iv_right )
     /*----------------------------------------------------------------------*/
     /* write info about new splitted intervals                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv_left  ... pointer to new left hand interval                     */
     /*   iv_right ... pointer to new right hand interval                    */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);       COOKIE_CHECK(gen,CK_TDRGW_GEN,RETURN_VOID);
  CHECK_NULL(iv_left,RETURN_VOID);   COOKIE_CHECK(iv_left,CK_TDRGW_IV,RETURN_VOID);

  if (iv_right == NULL) iv_right = iv_left;

  log = unur_get_stream();

  fprintf(log,"%s: inserted point:\n",gen->genid);
  fprintf(log,"%s: x = %g, logf(x) = %g, dlogf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->logfx, iv_right->dlogfx, iv_right->sq);
  fprintf(log,"%s: new intervals:\n",gen->genid);
  fprintf(log,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(log,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(log,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);

  fprintf(log,"%s: left interval:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv_left->Asqueeze,
	  iv_left->Asqueeze*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv_left->Ahat - iv_left->Asqueeze),
	  (iv_left->Ahat - iv_left->Asqueeze) * 100./GEN->Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv_left->Ahat - iv_left->Ahatr,
	  iv_left->Ahatr,
	  iv_left->Ahat * 100./GEN->Atotal);

  if (iv_left == iv_right)
    fprintf(log,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(log,"%s: right interval:\n",gen->genid);
    fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,
	    iv_right->Asqueeze*100./GEN->Atotal);
    fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    (iv_right->Ahat - iv_right->Asqueeze),
	    (iv_right->Ahat - iv_right->Asqueeze) * 100./GEN->Atotal);
    fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	    iv_right->Ahat - iv_right->Ahatr,
	    iv_right->Ahatr,
	    iv_right->Ahat * 100./GEN->Atotal);
  }

  fprintf(log,"%s: total areas:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (GEN->Atotal - GEN->Asqueeze) * 100./GEN->Atotal);
  fprintf(log,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);

  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdrgw_debug_split_stop() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
