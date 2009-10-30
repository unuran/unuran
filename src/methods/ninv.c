/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    numerical inversion of cumulative distribution function      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the CDF                                                   *
 *      newton's method: additional pointer to the PDF                       *
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
 *   REFERENCES:                                                             *
 *   [1] Neumaier A. (to be published): Introduction to numerical analysis,  *
 *       Cambridge University Press                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Numerical inversion is a method for generating random variables          *
 *  using the CDF (and in case of newton's method the PDF).                  *
 *                                                                           *
 *  THEOREM:                                                                 *
 *     Let X be a random variable with CDF F(x).                             *
 *     Then the F(X) are  uniformly distributed.                             *
 *                                                                           *
 *  COROLLARY:                                                               *
 *     Starting with uniformly distributed random variables U,               *
 *     the F^(-1)(U) have F(x) as CDF.                                       *
 *                                                                           *
 *  Starting with an U, the task is to find a X fulfilling:                  *
 *    F(X) - U = 0.                                                          *
 *                                                                           *
 *  Numerical algorithms to find zeros that are used in NINV are variants of * 
 *  newton's method (damped newton to guarantee improvement) and             *
 *  the regula falsi ( stabilized regula falsi preserving sign change; at    *
 *  first an interval with sign change is determined).                       *
 *                                                                           *
 *  In both cases it is possible to specify the maximal number of            *
 *  iterations, a desired accuracy in X and starting values for the          *
 *  algorithms.                                                              *
 *  Instead of starting values it is also possible to use a table            *
 *  containing suitable starting values.                                     *
 *  If neither the table nor explicit starting values are used,              *
 *  NINV chooses as starting values:                                         *
 *     newton's method:  x:     CDF(x) = 0.5                                 *
 *     regula falsi:     x1,x2: CDF(x1) = 1 - CDF(x2) = 0.05                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "ninv.h"
#include "ninv_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* Starting interval for Regula Falsi includes this percentage of all        */
/* univariate random numbers (must be > 0. and < 1.)                         */
#define INTERVAL_COVERS  (0.5)

/* maximum number of steps to find sign change in Regula Falsi               */
#define MAX_STEPS (100)

/* STEPFAC* (s[1]-s[0]) is used as first step length for finding sign change */
#define STEPFAC  (0.4)

/* for i > I_CHANGE_TO_BISEC Regula Falsi is always replaced by bisection    */
#define I_CHANGE_TO_BISEC (50)

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define NINV_VARFLAG_NEWTON   0x1u   /* use Newton's method                  */
#define NINV_VARFLAG_REGULA   0x2u   /* use regula falsi (default)           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NINV_DEBUG_REINIT    0x00000002u   /* print parameters after reinit  */
#define NINV_DEBUG_TABLE     0x00000010u   /* print table                    */
#define NINV_DEBUG_CHG       0x00001000u   /* print changed parameters       */
#define NINV_DEBUG_SAMPLE    0x01000000u   /* trace sampling                 */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NINV_SET_MAX_ITER     0x001u   /* number of maximal iterations       */
#define NINV_SET_X_RESOLUTION 0x002u   /* maximal tolerated relative x-error */
#define NINV_SET_U_RESOLUTION 0x004u   /* maximal tolerated (abs.) u-error   */
#define NINV_SET_START        0x008u   /* intervals at start (left/right)    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_sample_regula( struct unur_gen *gen );
static double _unur_ninv_sample_newton( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_regula( const struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* algorithm: regula falsi                                                   */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_newton( const struct unur_gen *gen, double u);
/*---------------------------------------------------------------------------*/
/* algorithm: newton method                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_compute_start( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* get starting points for numerical inversion                               */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create the table with starting points                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into LOG file               */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample_regula( const struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling (regula falsi).                                            */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample_newton( const struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling (newton's method).                                         */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_chg_truncated( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* trace changes of the truncated domain.                                    */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_ninv_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_ninv_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_ninv_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */
#define CDF(x)    _unur_cont_CDF((x),(gen->distr))    /* call to CDF         */

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CONT *
_unur_ninv_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    return _unur_ninv_sample_newton;
  case NINV_VARFLAG_REGULA:
  default:
    return _unur_ninv_sample_regula;
  }
} /* end of _unur_ninv_getSAMPLE() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_ninv_new( const struct unur_distr *distr )
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

  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }

  /* if default variant is Newton's method, then we also need the PDF ! */

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_ninv_par) );
  COOKIE_SET(par,CK_NINV_PAR);

  /* copy input */
  par->distr   = distr;            /* pointer to distribution object         */

  /* set default values */
  PAR->max_iter  = 100;            /* maximal number of iterations           */
  PAR->x_resolution = 1.0e-8;      /* maximal tolerated relative x-error     */
  PAR->u_resolution = -1.;         /* maximal tolerated u-error -- DISABLED  */

  /* starting points for numerical inversion */
  PAR->s[0]      = 0.0;     /* regula falsi: left boundary of starting interval
			      newton: starting point                         */
  PAR->s[1]      = 0.0;     /* regula falsi: right boundary of starting interval
			      newton: not used                               */
  /* If s1 and s2 are equal a defaults are used, see below */

  PAR->table_on  = FALSE;   /* Do not use a table for starting points
			      by default.                                    */
 
  par->method   = UNUR_METH_NINV;          /* method and default variant     */
  par->variant  = NINV_VARFLAG_REGULA;     /* Use regula falsi as default 
					      method                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_ninv_init;

  return par;

} /* end of unur_ninv_new() */

/*****************************************************************************/

int
unur_ninv_set_usenewton( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use Newton's method                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   /* use regula falsi instead  */
    return UNUR_ERR_DISTR_REQUIRED;
 }

  /* store date */
  par->variant = NINV_VARFLAG_NEWTON;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_usenewton() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_useregula( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use regula falsi                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* store date */
  par->variant = NINV_VARFLAG_REGULA;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_useregula() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_max_iter( struct unur_par *par, int max_iter )
     /*----------------------------------------------------------------------*/
     /* set number of maximal iterations                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximal iterations");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_iter = max_iter;

  /* changelog */
  par->set |= NINV_SET_MAX_ITER;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_max_iter() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_max_iter( struct unur_gen *gen, int max_iter )
     /*----------------------------------------------------------------------*/
     /* change number of maximal iterations                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(gen->genid, UNUR_ERR_PAR_SET, "maximal iterations");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  GEN->max_iter = max_iter;

  /* changelog */
  gen->set |= NINV_SET_MAX_ITER;

  /* o.k.  */
  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_max_iter() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_x_resolution( struct unur_par *par, double x_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated relative x-error                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   x_resolution ... x-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (x_resolution > 0. && x_resolution < 2.*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }

  /* store date */
  PAR->x_resolution = x_resolution;

  /* changelog */
  par->set |= NINV_SET_X_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_x_resolution( struct unur_gen *gen, double x_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated relative x-error                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   x_resolution ... x-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (x_resolution > 0. && x_resolution < DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"x-resolution too small");
    x_resolution = 2.*DBL_EPSILON;
  }

  /* store date */
  GEN->x_resolution = x_resolution;

  /* changelog */
  gen->set |= NINV_SET_X_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated u-error                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... u-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* check new parameter for generator */
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }

  /* store date */
  PAR->u_resolution = u_resolution;

  /* changelog */
  par->set |= NINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_u_resolution( struct unur_gen *gen, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal tolerated u-error                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   u_resolution ... u-error                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (u_resolution > 0. && u_resolution < 5*DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"u-resolution too small");
    u_resolution = 1.e-15;
  }

  /* store date */
  GEN->u_resolution = u_resolution;

  /* changelog */
  gen->set |= NINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_start( struct unur_par *par, double s1, double s2 )
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /* arguments that are not used by method are ignored.                   */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   s1    ... left starting point                                      */
     /*   s2    ... right starting point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  /* store date */
  if ( s1 <= s2 ) {
     PAR->s[0] = s1;
     PAR->s[1] = s2;
  }
  else {
     PAR->s[0] = s2;
     PAR->s[1] = s1;
  }

  /* changelog */
  par->set |= NINV_SET_START;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_start() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_start( struct unur_gen *gen, double s1, double s2 )
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /* arguments that are used by method are ignored.                       */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   s1    ... left starting point                                      */
     /*   s2    ... right starting point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* store date */
  if ( s1 <= s2 ) {
     GEN->s[0] = s1;
     GEN->s[1] = s2;
  }
  else {
     GEN->s[0] = s2;
     GEN->s[1] = s1;
  }

  /* disable table (we now want to use only two starting points) */
  GEN->table_on = FALSE;

  /* compute these points */
  _unur_ninv_compute_start(gen);
  /* this call should not fail */

  /* changelog */
  gen->set |= NINV_SET_START;

  return UNUR_SUCCESS;

} /* end of unur_ninv_chg_start() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_table( struct unur_par *par, int tbl_pnts )
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* the function unur_ninv_set_start() is overruled                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   tbl_pnts ... number of table points                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NINV );

  PAR->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  PAR->table_on = TRUE;

  return UNUR_SUCCESS;

} /* end of unur_ninv_set_table() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_table( struct unur_gen *gen, int tbl_pnts )
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* set somewhere else will be ignored                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   tbl_pnts ... number of table points                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /*  free(GEN->table);   not freed, because realloc() is used */ 
  /*  free(GEN->f_table); not freed, because realloc() is used */
  GEN->table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;

  result = _unur_ninv_create_table(gen); 

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_CHG) 
    if (result==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
#endif
  
  return result;

} /* end of unur_ninv_chg_table() */

/*---------------------------------------------------------------------------*/

int 
unur_ninv_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  double Umin, Umax;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object( gen, NINV, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* set bounds of U -- in respect to given bounds */
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;

  /* check result */
  if (Umin > Umax) {
    /* this is a serios error that should not happen */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  if (_unur_FP_equal(Umin,Umax)) {
    /* CDF values very close */
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      /* this is very bad */
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }

  /* copy new boundaries into generator object */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_CHG) 
    _unur_ninv_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_ninv_chg_truncated() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_ninv_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* check variant */
  if (par->variant == NINV_VARFLAG_NEWTON && ! par->DISTR_IN.pdf) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    par->variant = NINV_VARFLAG_REGULA;   /* use regula falsi instead  */
  }

  /* create a new empty generator object */    
  gen = _unur_ninv_create(par);
  _unur_par_free(par);
  if (!gen) { return NULL; }
  
  /* check parameters */
  if (_unur_ninv_check_par(gen) != UNUR_SUCCESS) {
    _unur_ninv_free(gen); return NULL;
  }

  /* compute starting points for numerical inversion */
  if (GEN->table_on) {
    /* use a table */
    if (_unur_ninv_create_table(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }
  else {
    /* use given points or percentiles */
    if (_unur_ninv_compute_start(gen)!=UNUR_SUCCESS) {
      _unur_ninv_free(gen); return NULL;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ninv_debug_init(gen);
#endif

  return gen;

} /* end of _unur_ninv_init() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_reinit( struct unur_gen *gen )
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
  if ( (rcode = _unur_ninv_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* compute normalization constant for standard distribution */
  if (DISTR.upd_area != NULL)
    if ((DISTR.upd_area)(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cannot compute normalization constant");
      return UNUR_ERR_GEN_DATA;
    }

  /* regenerate table */
  if (GEN->table != NULL)
    rcode = _unur_ninv_create_table(gen);

  else /* or compute starting points */
    rcode = unur_ninv_chg_start( gen, 0., 0. );

  /* (re)set sampling routine */
  SAMPLE = _unur_ninv_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NINV_DEBUG_REINIT) {
    _unur_distr_cont_debug( gen->distr, gen->genid );
    if (rcode==UNUR_SUCCESS) _unur_ninv_debug_start( gen );
  }
#endif

  return UNUR_SUCCESS;
} /* end of _unur_ninv_reinit() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_ninv_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NINV_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_ninv_getSAMPLE(gen);
  gen->destroy = _unur_ninv_free;
  gen->clone = _unur_ninv_clone;
  gen->reinit = _unur_ninv_reinit;

  /* copy parameters into generator object */
  GEN->max_iter = PAR->max_iter;      /* maximal number of iterations          */
  GEN->x_resolution = PAR->x_resolution; /* maximal tolerated relative x-error */
  GEN->u_resolution = PAR->u_resolution; /* maximal tolerated u-error          */
  GEN->table_on = PAR->table_on;      /* useage of table for starting points   */
  GEN->table_size = PAR->table_size;  /* number of points for table            */
  GEN->s[0] = PAR->s[0];              /* starting points                       */
  GEN->s[1] = PAR->s[1];

  /* init pointer */
  GEN->table = NULL;
  GEN->f_table = NULL;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_ninv_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_ninv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_check_par( struct unur_gen *gen )
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
  /* checking x-error or checking u-error must be enabled */
  if ( GEN->x_resolution < 0. && GEN->u_resolution < 0. ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"both x-resolution and u-resolution negativ. using defaults.");
    GEN->x_resolution = 1.e-8;
  }

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* set bounds of U -- in respect to given bounds                          */
  GEN->CDFmin = GEN->Umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN->CDFmax = GEN->Umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;

  if (_unur_FP_greater(GEN->CDFmin, GEN->CDFmax)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    return UNUR_ERR_GEN_DATA;
  }

  return UNUR_SUCCESS;
} /* end of _unur_ninv_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ninv_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_ninv_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NINV_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy additional data for generator object */
  if (GEN->table) {
    CLONE->table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->table, GEN->table, GEN->table_size * sizeof(double) );
    CLONE->f_table = _unur_xmalloc( GEN->table_size * sizeof(double) );
    memcpy( CLONE->f_table, GEN->f_table, GEN->table_size * sizeof(double) );
  }

  return clone;

#undef CLONE
} /* end of _unur_ninv_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free tables */
  if (GEN->table)   free(GEN->table);
  if (GEN->f_table) free(GEN->f_table);

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_ninv_free() */

/*****************************************************************************/

double
_unur_ninv_sample_regula( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use regula falsi)                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_regula( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
} /* end of _unur_ninv_sample_regula() */

/*---------------------------------------------------------------------------*/

double 
_unur_ninv_sample_newton( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use newtons method)                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_newton( gen, 
         GEN->Umin + (_unur_call_urng(gen->urng)) * (GEN->Umax - GEN->Umin) );
}

/*---------------------------------------------------------------------------*/

double
unur_ninv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* get approximate value of inverse CDF at u approximately              */
     /* (user call)                                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF (0<=u<=1)                         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double x;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);

  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  /* = NaN */
  }
  
  /* compute inverse CDF */
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    x = _unur_ninv_newton(gen,u);
    break;
  case NINV_VARFLAG_REGULA:
  default:
    x = _unur_ninv_regula(gen,u);
    break;
  }

  /* validate range */
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];

  return x;

} /* end of unur_hinv_eval_approxinvcdf() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

double 
_unur_ninv_regula( const struct unur_gen *gen, double u )
     /*---------------------------------------------------------------------*/
     /*   algorithm: regula falsi with bisection steps                      */
     /*                                                                     */
     /*   parameters:                                                       */
     /*      gen ... pointer to generator object                            */
     /*      u   ... random number (uniform distribution)                   */
     /*   return:                                                           */
     /*     double (sample from random variate)                             */
     /*                                                                     */
     /*   error:                                                            */
     /*     return INFINITY                                                 */
     /*                                                                     */
     /*   Remark:                                                           */
     /*     The routine computes the root of CDF(x)-u                       */
     /*---------------------------------------------------------------------*/
{ 
  double x1, x2, a, xtmp;/* points for regular falsi                        */
  double f1, f2,fa, ftmp;/* function values at x1, x2, xtmp                 */
  double length;         /* oriented length of the interval with sign change*/
  double lengthabs;      /* absolute length of interval                     */
  double lengthsgn;      /* orientation of the Intervalls                   */
  double step;           /* enlarges interval til sign change found         */
  double dx;             /* RF-stepsize                                     */
  int count_nosc = 0;    /* counter for  "no sign change occured"           */
  int i;                 /* loop variable, index                            */
  int step_count;        /* counts number of steps finding sign change      */
  double min_step_size;  /* minimal step size for regula falsi              */
  double rel_u_resolution; /* relative u resolution                         */
  double abs_x_resolution; /* absolute x resolution                         */
  int x_goal, u_goal;    /* whether precision goal is reached               */

  /* check arguments */
  CHECK_NULL(gen, INFINITY);  COOKIE_CHECK(gen, CK_NINV_GEN, INFINITY);

  /* compute relative u resolution */
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
		       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
		       INFINITY );

  /* compute absolute x resolution eps^2 used close to 0 */
  abs_x_resolution = GEN->x_resolution * GEN->x_resolution;
  
  /* -- 1. initialize starting interval -- */

  if (GEN->table_on) {
    /* -- 1a. use table -- */

    /* 0 <= i <= table_size-2  */
    if ( _unur_FP_same(GEN->CDFmin, GEN->CDFmax) ) {
      /* CDF values in table too close, so we use median point since */
      /* there is no difference between CDF values.  */
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (u - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }

    /* set starting point for regular falsi */
    if ( ! _unur_FP_is_minus_infinity(GEN->table[i]) ){
      x1 = GEN->table[i];
      f1 = GEN->f_table[i]; 
    }
    else{
      x1 = GEN->table[i+1] + (GEN->table[i+1] - GEN->table[i+2]);
      f1 = CDF(x1);
    }

    if( ! _unur_FP_is_infinity(GEN->table[i+1]) ){
      x2 = GEN->table[i+1];
      f2 = GEN->f_table[i+1];
    }
    else{
      x2 = GEN->table[i] + (GEN->table[i] - GEN->table[i-1]);
      f2 = CDF(x2);
    }
  }

  else { 
    /* 1b. -- no table available -- */
    x1 =  GEN->s[0];      /* left boudary of interval */
    f1 =  GEN->CDFs[0];
    x2 =  GEN->s[1];      /* right boudary of interval*/   
    f2 =  GEN->CDFs[1];
  }

  /*  -- 1c. check for ordering of starting points -- */

  if ( x1 >= x2 ) { 
    xtmp = x1; ftmp = f1;
    x1   = x2; f1   = f2;
    x2 = xtmp + fabs(xtmp)*DBL_EPSILON;
    f2 = CDF(x2); 
  }

  /* -- 1d. check for boundary of truncated domain -- */
 
  /* in case of truncated domain there might be better starting points */
  /* ! no problems with INFINITY !  */
  if ( x1 < DISTR.trunc[0] || x1 >= DISTR.trunc[1] ){
    x1 = DISTR.trunc[0];
    f1 = GEN->Umin;    /* = CDF(x1) */
  }
  if ( x2 > DISTR.trunc[1] || x2 <= DISTR.trunc[0] ){
    x2 = DISTR.trunc[1];
    f2 = GEN->Umax;    /* = CDF(x2) */
  }

  /* -- 1z. compute function values at interval boundaries -- */
  f1 -= u;  f2 -= u;


  /* -- 2. search for enclosing bracket (interval where f changes signs) -- */
 
  step = (GEN->s[1]-GEN->s[0]) * STEPFAC;
  step_count = 0;
  while ( f1*f2 > 0. ) {
    /* interval too small -> make it bigger ( + 2^n * gap ) -- */
    if ( f1 > 0. ) { /* lower boundary too big */    
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1  = CDF(x1) - u;
    }
    else {         /* upper boundary too small */
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2  = CDF(x2) - u;
    }

    /* increase step width */
    if (step_count < MAX_STEPS) {
      ++step_count;
      step *= 2.;
      /* safe guard for the case where (GEN->s[1]-GEN->s[0]) is very small*/
      if( step_count > 20 && step < 1.) step = 1.; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
                  "Regula Falsi cannot find interval with sign change");
      return ( (f1>0.) ? DISTR.trunc[0] : DISTR.trunc[1] );
    }
  }

  /* -- 2a. sign changes always within [a, x2] -- */
  a = x1; fa = f1; 


  /* -- 3. secant step, preserve sign change -- */

  /* secant step, preserve sign change */
  for (i=0; TRUE; i++) {

    /* -- 3a. check sign of f1 and f2 -- */

    if ( f1*f2 < 0.) { 
      /* sign change found     */
      count_nosc = 0;   /* reset counter */
      /* f2 always less (better), otherwise exchange */
      if ( fabs(f1) < fabs(f2) ) {
	xtmp = x1; ftmp = f1;
	x1 = x2;   f1 = f2;
	x2 = xtmp; f2 = ftmp;
      }
      /* sign changes within [a, x2] */
      a = x1; fa = f1;
    }
    else {
      /* increment counter for "no sign change occured" */
      count_nosc++;
      /* must not update 'a' */
    }

    /* length of enclosing bracket */
    length = x2 - a;                       /* oriented length */
    lengthabs = fabs(length);              /* absolute length */
    lengthsgn = (length < 0.) ? -1. : 1.;

    /* -- 3b. check stopping criterions -- */

    if (i >= GEN->max_iter)
      /* maximum number of iterations reached --> abort */
      break;

    if ( GEN->x_resolution > 0. ) {
      /* check x-error */
      if ( _unur_iszero(f2) ||                          /* exact hit */ 
	   lengthabs <= GEN->x_resolution * fabs(x2) || /* relative x resolution reached */
	   lengthabs <= abs_x_resolution ) {            /* absolute x resolution reached close to 0 */
	x_goal = TRUE;
      }
      else if ( _unur_FP_same(fa, f2) ) {
	/* flat region  */
	_unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		      "flat region: accuracy goal in x cannot be reached");
	x_goal = TRUE;
      }
      else
	x_goal = FALSE;
    }
    else {
      /* no check */
      x_goal = TRUE;
    }

    if ( GEN->u_resolution > 0. ) {
      /* check u-error */
      u_goal = (fabs(f2) <= rel_u_resolution);    /* relative u resolution */
      if ( _unur_FP_same(a, x2) ) {
	/* sharp peak or pole */
	_unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		      "sharp peak or pole: accuracy goal in u cannot be reached");
	u_goal = TRUE;
      }
    }
    else {
      u_goal = TRUE;
    }

    /* goal reached ? */
    if (x_goal && u_goal)
      /*finished*/
      break;

    /* -- 3c. try secant step -- */

    /* step size: secant or bisection step */
    dx = (_unur_FP_same(f1,f2)) ? length/2. : f2*(x2-x1)/(f2-f1) ;  

    /* minimal step size */
    if (GEN->u_resolution < 0.) 
      /* we only look at the x-error, so we do not need shorter steps */
      min_step_size = fabs(x2) * GEN->x_resolution;
    else
      min_step_size = lengthabs * DBL_EPSILON; 

    if ( fabs(dx) < min_step_size ) {
      dx = lengthsgn * 0.99 * min_step_size;

      while ( x2 == x2 - dx ){ /* dx too small  */
        if ( dx != 2.*dx)    /* near limit of calculations */
          dx = 2.*dx;
        else
          dx = length/2.;    /* bisection step   */
      }
    }

    /* bisection step if:                             */  
    /* no sign change   || step leads out of interval */
    if ( count_nosc > 1 || i > I_CHANGE_TO_BISEC ||
	 (lengthabs-GEN->x_resolution*fabs(x2))/(dx*lengthsgn) <= 1. )
      dx = length/2.; /* bisection step        */
  
    /* -- 3c. update point -- */    
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
    
  }  /* end of for-loop */

  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");
 
  /* ensure location within given (truncated) domain */
  x2 = _unur_max( x2, DISTR.trunc[0]);
  x2 = _unur_min( x2, DISTR.trunc[1]);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file (in case error) */
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample_regula( gen,u,x2,f2,i );
#endif

  /* finished */
  return x2;

} /* end of _unur_ninv_sample_regula()  */


/*---------------------------------------------------------------------------*/

double
_unur_ninv_newton( const struct unur_gen *gen, double U )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use Newton's method)                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*     U   ... random number (uniform distribution)                     */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double x;           /* point for netwon-iteration                   */
  double fx;          /* cdf at x                                     */
  double dfx;         /* pdf at x                                     */
  double fxabs;       /* absolute valuo of fx                         */
  double xtmp, fxtmp; /* temprary variables for x and fx              */
  double xold;        /* remember last values for stopping criterion  */
  double fxtmpabs;    /* fabs of fxtmp                                */
  double damp;        /* damping factor                               */
  double step;        /* helps to escape from flat regions of the cdf */
  int i;              /* counter for for-loop, index                  */
  int flat_count;     /* counter of steps in flat region              */
  double rel_u_resolution; /* relative u resolution                   */
  double abs_x_resolution; /* absolute x resolution                   */
  int x_goal, u_goal; /* whether precision goal is reached            */

  /* maximal number of steps to leave flat region */
  const int MAX_FLAT_COUNT = 40;
    
  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NINV_GEN,INFINITY);

  /* compute relative u resolution */
  rel_u_resolution = ( (GEN->u_resolution > 0.) ? 
                       (GEN->Umax - GEN->Umin) * GEN->u_resolution :
                       INFINITY );

  /* compute absolute x resolution eps^2 used close to 0 */
  abs_x_resolution = GEN->x_resolution * GEN->x_resolution;
  
  /* -- 1. initialize starting interval -- */

  if (GEN->table_on) {
    /* -- 1a. use table -- */

    /* 0 <= i <= table_size-2  */
    if ( _unur_FP_same(GEN->CDFmin,GEN->CDFmax) ) {
      /* CDF values in table too close, so we use median point since */
      /* there is no difference between CDF values.                  */
      i = GEN->table_size/2;
    }
    else {
      i = (int) ( GEN->table_size * (U - GEN->CDFmin) / (GEN->CDFmax - GEN->CDFmin) );
      if (i<0) i = 0;
      else if (i > GEN->table_size - 2) i = GEN->table_size - 2;
    }

    if (_unur_FP_is_infinity(GEN->table[i+1])) {
      x  = GEN->table[i];
      fx = GEN->f_table[i];
    }
    else {
      x  = GEN->table[i+1];
      fx = GEN->f_table[i+1];
    }

  }

  else { 
    /* 1b. -- no table available -- */
    x  = GEN->s[0];
    fx = GEN->CDFs[0];
  }

  /* -- 1c. check for boundary of truncated domain -- */

  if ( x < DISTR.trunc[0] ){
    x  = DISTR.trunc[0];
    fx = GEN->Umin;    /* = CDF(x) */
  }
  else if ( x > DISTR.trunc[1] ){
    x  = DISTR.trunc[1];
    fx = GEN->Umax;    /* = CDF(x) */
  }

  /* -- 1z. compute values for starting point -- */

  fx   -= U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    /* there is no old value yet */

  damp = 2.;        /* to be halved at least once */  
  step = 1.;

  /* -- 2. Newton iteration -- */

  for (i=0; i < GEN->max_iter; i++) {

    flat_count = 0;
    while (_unur_iszero(dfx)) {   /* function flat at x */
      /* printf("step: %g, x: %g, fx: %g, dfx: %g\n",step, x, fx, dfx); */

      if (_unur_iszero(fx))  /* exact hit -> leave while-loop */
	break; 

      if (fx > 0.) {         /* try another x */
        xtmp = x - step; 
	xtmp = _unur_max( xtmp, DISTR.domain[0] );
      }
      else {
        xtmp  = x + step;
	xtmp = _unur_min( xtmp, DISTR.domain[1] );
      }

      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);

      if ( fxtmpabs < fxabs ) {       /* improvement, update x               */
	/* printf("fxabs: %g tmpabs: %g\n", fxabs, fxtmpabs); */
        step = 1.;     /* set back stepsize */
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmp*fx < 0. ) {     /* step was too large, don't update x  */
        step /= 2.;                      
      } 
      else {                          /* step was too short, update x        */
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  

      dfx   = PDF(x);
      fxabs = fabs(fx);

      if (flat_count < MAX_FLAT_COUNT)
	flat_count++;
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_SAMPLING,
		    "Newton's method cannot leave flat region");
	x = _unur_max( x, DISTR.trunc[0]);
	x = _unur_min( x, DISTR.trunc[1]);
	return x;
      }
    }   /* end of while-loop, (leaving flat region) */
    
    step = 1.;   /* set back stepsize */

    if (_unur_iszero(fx))  /* exact hit -> finished */
      break;


    if (_unur_isfinite(dfx)) {
      do {    /* newton-step  (damped if nececcary) */
	damp /= 2.;
	xtmp = x - damp * fx/dfx;
	/* make sure that new point is inside (truncated) domain */
	xtmp = _unur_min( xtmp, DISTR.trunc[1] );
	xtmp = _unur_max( xtmp, DISTR.trunc[0] );
	fxtmp = CDF(xtmp) - U;
      } while (fabs(fxtmp) > fxabs * (1.+UNUR_SQRT_DBL_EPSILON));   /* no improvement */
    }
    else {
      /* we cannot use Newton's rule if the derivative is not finite. */
      /* this happens when we hit a pole of the PDF.                  */
      /* use a bisection step instead.                                */
      xtmp = 0.5*(x + xold);
      fxtmp = CDF(xtmp) - U;
    }
    
    /* updation variables according to newton-step      */
    damp  = 2.;       /* set back factor for damping    */
    xold  = x;        /* remember last value of x       */
    x     = xtmp;     /* update x                       */
    fx    = fxtmp;    /* update function value at x     */
    dfx   = PDF(x);   /* update derivative sof fx at x  */
    fxabs = fabs(fx); /* update absolute value of fx    */
 
    /* -- 2z. check stopping criterions -- */

    if ( GEN->x_resolution > 0. ) {
      /* check x-error */
      if ( _unur_iszero(fx) ||                            /* exact hit */ 
	   fabs(x-xold) <= GEN->x_resolution * fabs(x) || /* relative x resolution reached */
           fabs(x-xold) <= abs_x_resolution ) {           /* absolute x resolution reached close to 0 */
 	x_goal = TRUE;
      }
      else
        x_goal = FALSE;
    }
    else {
      /* no check */
      x_goal = TRUE;
    }

    if ( GEN->u_resolution > 0. ) {
      /* check u-error */
      if ( fabs(fx) <= rel_u_resolution ) {    /* relative u resolution */
      	u_goal = TRUE;
      }
      else if ( _unur_FP_same(xold, x) ) {
        /* sharp peak or pole */
        _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
                      "sharp peak or pole: accuracy goal in u cannot be reached");
        u_goal = TRUE;
      }
      else
        u_goal = FALSE;

    }
    else {
      u_goal = TRUE;
    }

    /* goal reached ? */
    if (x_goal && u_goal)
      /*finished*/
      break;
  }  /* end of for-loop  (MAXITER reached -> finished) */


  if (i >= GEN->max_iter)
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,
		  "max number of iterations exceeded: accuracy goal might not be reached");

  /* make sure that result is within boundaries of (truncated) domain */
  x = _unur_max( x, DISTR.trunc[0]);
  x = _unur_min( x, DISTR.trunc[1]);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file (in case error) */
  if (gen->debug & NINV_DEBUG_SAMPLE)
    _unur_ninv_debug_sample_newton(gen, U, x, fx, i);
#endif
  
  return x;

} /* end of _unur_ninv_sample_newton() */

/*---------------------------------------------------------------------------*/
int
_unur_ninv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a table for starting points and a table with                  */
     /* the corresponding function values                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;
  double x;
  int table_size = GEN->table_size;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);

  GEN->table    = _unur_xrealloc( GEN->table,   table_size * sizeof(double));
  GEN->f_table  = _unur_xrealloc( GEN->f_table, table_size * sizeof(double));

  /* get arbitrary points */
  GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
  GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
  GEN->CDFs[0]  = CDF(GEN->s[0]);
  GEN->CDFs[1]  = CDF(GEN->s[1]);

  /* table can't be used to calculate itself */
  GEN->table_on = FALSE;
  
  /* calculation of the tables   */

  /* left and right boundary */
  GEN->table[0]              = DISTR.domain[0];
  GEN->f_table[0]            = GEN->CDFmin;    /* CDF(DISTR.domain[0]) */
  GEN->table[table_size-1]   = DISTR.domain[1];
  GEN->f_table[table_size-1] = GEN->CDFmax;    /* CDF(DISTR.domain[1]) */
  
  /* all the other points. we compute these points from boundary to center */
  for (i=1; i<table_size/2; i++){

    /* compute table point and CDF at table point */
    x = GEN->CDFmin + i * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[i]   = _unur_ninv_regula(gen,x);
    GEN->f_table[i] = CDF(GEN->table[i]);

    x = GEN->CDFmin + (table_size-i-1) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size-1-i] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size-1-i] = CDF(GEN->table[table_size-1-i]);

    /* set new starting points for computing table points in next step */
    if (GEN->table[i] > -INFINITY) {
      GEN->s[0] = GEN->table[i];
      GEN->CDFs[0] = GEN->f_table[i];
    }
    if (GEN->table[table_size-1-i] < INFINITY) {
      GEN->s[1] = GEN->table[table_size-1-i];
      GEN->CDFs[1] = GEN->f_table[table_size-1-i];
    }

  }  /* end of for()                                                     */

  /* the median point (only if table_size is odd) */
  if (table_size & 1) { 
    x = GEN->CDFmin + (table_size/2) * (GEN->CDFmax - GEN->CDFmin) / (table_size-1.);  
    GEN->table[table_size/2] = _unur_ninv_regula(gen,x);
    GEN->f_table[table_size/2] = CDF(GEN->table[table_size/2]);
  }  

  /* calculation of tables finished  */

  GEN->table_on = TRUE;

  /* o.k. */
  return UNUR_SUCCESS;

}  /* end of _unur_ninv_create_table() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_compute_start( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting points for numerical inversion.                     */
     /* use 5% percentiles.                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double u;

  /* check arguments */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  _unur_check_gen_object(gen, NINV, UNUR_ERR_GEN_INVALID);

  if( GEN->table_on )
    /* we have a table --> nothing to do */
    return UNUR_SUCCESS;

  if( !_unur_FP_same(GEN->s[0],GEN->s[1]) ) {
    /* use given starting points (indicated by s[0] != s[1]) --> nothing to do */
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    return UNUR_SUCCESS;
  }

  switch (gen->variant) {

  case NINV_VARFLAG_REGULA:

    /* get arbitrary points */
    GEN->s[0] = _unur_max( DISTR.domain[0], -10.);
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]);
    GEN->CDFs[1] = CDF(GEN->s[1]);    

    /* left percentile */
    u = GEN->CDFmin + 0.5*(1.-INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[0] = _unur_ninv_regula(gen,u);
    GEN->CDFs[0] = CDF(GEN->s[0]);

    /* right percentile */
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    u = GEN->CDFmin + 0.5*(1.+INTERVAL_COVERS)*(GEN->CDFmax-GEN->CDFmin);
    GEN->s[1] = _unur_ninv_regula(gen,u);
    GEN->CDFs[1] = CDF(GEN->s[1]);
    
    break;    /* case REGULA end */

  case NINV_VARFLAG_NEWTON:

    /* get arbitrary points */
    GEN->s[0] = _unur_max( DISTR.domain[0], -9.987655 );
    GEN->s[1] = _unur_min( DISTR.domain[1], GEN->s[0]+20. );
    GEN->CDFs[0] = CDF(GEN->s[0]); 
    GEN->CDFs[1] = CDF(GEN->s[1]);

    /* median */
    u = 0.5 * (GEN->CDFmin + GEN->CDFmax);
    GEN->s[0] = _unur_ninv_regula(gen, u); /* reglua more stable */
    GEN->CDFs[0] = CDF(GEN->s[0]);

    break;    /* case NEWTON end */

  default:

    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;

  }  /* end of switch  */

  /* o.k. */
  return UNUR_SUCCESS;

}  /* end of _unur_ninv_compute_start() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ninv (numerical inversion of CDF)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_ninv_sample",gen->genid);
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    fprintf(LOG,"_newton\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    fprintf(LOG,"_regula\n");
    break;
  }
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_ninv_debug_start(gen);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print starting points or table for algorithms into LOG file          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (GEN->table_on) {
    fprintf(LOG,"%s: use table (size = %d)\n",gen->genid,GEN->table_size);
    if (gen->debug & NINV_DEBUG_TABLE)
      for (i=0; i<GEN->table_size; i++)
	fprintf(LOG,"%s:\tx = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->table[i],GEN->f_table[i]);
  }
  else { /* no table */
    fprintf(LOG,"%s: starting points:\n",gen->genid);
    fprintf(LOG,"%s:\ts[0] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[0],GEN->CDFs[0]);
    if (gen->variant & NINV_VARFLAG_REGULA)
      fprintf(LOG,"%s:\ts[1] = %12.6g, F(x) = %10.8f\n",gen->genid,GEN->s[1],GEN->CDFs[1]);
  }

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_start() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_sample_regula( const struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (regula falsi)                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d iterations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN->max_iter);

} /* end of _unur_ninv_debug_sample_regula() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_sample_newton( const struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (newton's method)                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d iterations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN->max_iter);

} /* end of _unur_ninv_debug_sample_newton() */

/*---------------------------------------------------------------------------*/

void 
_unur_ninv_debug_chg_truncated( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NINV_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(LOG,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);

} /* end of _unur_ninv_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_ninv_info( struct unur_gen *gen, int help )
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
  double n_iter;
  int samplesize = 10000;

  int use_newton = (gen->variant==NINV_VARFLAG_NEWTON) ? TRUE : FALSE;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = CDF");
  if (use_newton) 
    _unur_string_append(info," PDF");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n\n");
      
  /* method */
  _unur_string_append(info,"method: NINV (Numerical INVersion)\n");
  if (use_newton) 
    _unur_string_append(info,"   Newton method\n");
  else
    _unur_string_append(info,"   Regula falsi\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  n_iter = unur_test_count_pdf(gen,samplesize,FALSE,NULL)/(2.*samplesize);
  if (!use_newton) n_iter *= 2.;
  _unur_string_append(info,"   average number of iterations = %.2f  [approx.]\n", n_iter);

  if (GEN->table_on) {
    _unur_string_append(info,"   starting points = table of size %d\n", GEN->table_size);
  }
  else {
    _unur_string_append(info,"   starting points = ");
    if (use_newton) {
      _unur_string_append(info,"%g (CDF = %g)  %s\n", GEN->s[0], GEN->CDFs[0],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
    else {
      _unur_string_append(info,"%g, %g  (CDF = %g, %g)   %s\n",
			  GEN->s[0],GEN->s[1], GEN->CDFs[0],GEN->CDFs[1],
			  (gen->set & NINV_SET_START) ? "" : "[default]");
    }
  }
  _unur_string_append(info,"\n");


  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    if (use_newton) 
      _unur_string_append(info,"   usenewton\n");
    else
      _unur_string_append(info,"   useregula  [default]\n");

    _unur_string_append(info,"   x_resolution = %g  %s\n", GEN->x_resolution,
			(gen->set & NINV_SET_X_RESOLUTION) ? "" : "[default]");

    _unur_string_append(info,"   max_iter = %d  %s\n", GEN->max_iter,
			(gen->set & NINV_SET_MAX_ITER) ? "" : "[default]");

    /* Not displayed:
       int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
       int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
    */
    _unur_string_append(info,"\n");
  }


  /* Hints */
  if (help) {
    if (! (gen->set & NINV_SET_X_RESOLUTION) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase accuracy by decreasing \"x_resolution\".");
    if (! (gen->set & NINV_SET_MAX_ITER) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase \"max_iter\" if you encounter problems with accuracy.");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
