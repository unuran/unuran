/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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
 *   REFERENCES:                                                             *
 *   [1] Ahrens J. H. (1993): Sampling from general distributions by         *
 *       suboptimal division of domains,                                     *
 *       Grazer Math. Berichte 319, 30pp.                                    *
 *                                                                           *
 *   [2] Ahrens J. H. (1995): An one-table method for sampling from          *
 *       continuous and discrete distributions,                              *
 *       Computing 54(2), pp. 127-146                                        *
 *                                                                           *
 *   SEE ALSO:                                                               *
 *   [3] Gilks, W. R. and Wild,  P. (1992):                                  *
 *       Adaptive rejection sampling for Gibbs sampling,                     *
 *       Applied Statistics 41, pp. 337-348                                  *
 *                                                                           *
 *   [4] Zaman, A. (1996), Generation of Random Numbers from an Arbitrary    *
 *       Unimodal Density by Cutting Corners, unpublished manuskript         *
 *       available at http://chenab.lums.edu.pk/~arifz/                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <uniform/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "tabl.h"
#include "tabl_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

#define TABL_DEFAULT_COMPUTATION_LIMIT  1.e20

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TABL_VARFLAG_VERIFY       0x001u   /* flag for verifying mode           */

/* indicate how to split interval */
#define TABL_VARMASK_SPLIT        0x0f0u  /* split at        computation     convergence of hat */
#define TABL_VARFLAG_SPLIT_POINT  0x010u  /* sampled point    none            slowest          */
#define TABL_VARFLAG_SPLIT_MEAN   0x020u  /* mean point       slower          better           */
#define TABL_VARFLAG_SPLIT_ARC    0x040u  /* "arcmean"        very slow       very good for almost unbounded domain */

/* indicate if starting intervals have to be split */
#define TABL_VARFLAG_STP_A        0x100u  /* use equal area rule (SPLIT A in [1])   */

#define TABL_VARFLAG_USEDARS      0x200u  /* use main subdivisions (SPLIT B in [1]) 
                                             (= derandomized ARS)                   */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TABL_DEBUG_IV    0x00000100u /* show intervals                       */
#define TABL_DEBUG_A_IV  0x00000200u /* show intervals after split A, before split B */
#define TABL_DEBUG_DARS  0x00020000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TABL_SET_GUIDEFACTOR      0x001u
#define TABL_SET_SLOPES           0x004u
#define TABL_SET_AREAFRACTION     0x008u
#define TABL_SET_MAX_IVS          0x010u
#define TABL_SET_MAX_SQHRATIO     0x020u
#define TABL_SET_N_STP            0x040u
#define TABL_SET_BOUNDARY         0x080u
#define TABL_SET_USE_DARS         0x100u
#define TABL_SET_DARS_FACTOR      0x200u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TABL"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_tabl_sample( struct unur_gen *gen );
static double _unur_tabl_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals.                                               */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_starting_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals, use given slopes                              */
/*---------------------------------------------------------------------------*/
static int _unur_tabl_get_starting_intervals_from_mode( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals, use mode and domain                           */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_split_a_starting_intervals( struct unur_par *par, struct unur_gen *gen, struct unur_tabl_interval *iv_slope );
/*---------------------------------------------------------------------------*/
/* split starting intervals according to [1]                                 */
/* SPLIT A (equal areas rule)                                                */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_run_dars( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* run derandomized adaptive rejection sampling.                             */
/* (split starting intervals according to [1] SPLIT B,                       */
/* but instead of the iteration in [1] use "arcmean".                        */
/*---------------------------------------------------------------------------*/

static int
_unur_tabl_split_interval( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			   double x, double fx, unsigned split_mode );
/*---------------------------------------------------------------------------*/
/* split interval (replace old one by two new ones in same place)            */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed I.               */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init_finished( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed II.              */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print header before runniung derandomized adaptive rejection sampling.    */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_dars( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has run derandomized adaptive rejection sampling.   */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_intervals( const struct unur_gen *gen, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals.                                                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_tabl_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_tabl_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_tabl_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_tabl_par) );
  COOKIE_SET(par,CK_TABL_PAR);

  /* copy input */
  par->distr        = distr;     /* pointer to distribution object           */

  /* set default values */
  PAR->slopes        = NULL;      /* pointer to slopes of PDF                 */
  PAR->n_slopes      = 0;         /* number of slopes                         */

  PAR->n_starting_cpoints = 30;   /* number of starting points                */
  PAR->area_fract    = 0.1;       /* parameter for equal area rule (default from [1] ) */

  PAR->max_ivs       = 1000;      /* maximum number of intervals              */
  PAR->max_ratio     = 0.90;      /* bound for ratio  Atotal / Asqueeze       */

  PAR->guide_factor  = 1.; /* guide table has same size as array of intervals */

  PAR->darsfactor    = 0.99;   /* factor for (derandomized) ARS.
				 do not add a new construction point in a interval
				 where abiguous region is too small          */

  /* default boundary of compution area */
  PAR->bleft     = -TABL_DEFAULT_COMPUTATION_LIMIT;
  PAR->bright    = TABL_DEFAULT_COMPUTATION_LIMIT;

  par->method   = UNUR_METH_TABL;              /* indicate method            */
  par->variant  = (TABL_VARFLAG_SPLIT_MEAN |   /* variant: split at arc_mean */
		   TABL_VARFLAG_STP_A      |   /* run SPLIT A on slopes      */
		   TABL_VARFLAG_USEDARS    );  /* run DARS (SPLIT B) on slopes */


  par->set      = 0u;                      /* inidicate default parameters   */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_tabl_init;

  return par;

} /* end of unur_tabl_new() */

/*****************************************************************************/

int
unur_tabl_set_usedars( struct unur_par *par, int usedars )
     /*----------------------------------------------------------------------*/
     /* set flag for using DARS (derandomized adaptive rejection sampling).  */
     /* additionally the rule for splitting intervals can be set.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usedars   ... 0 = do not use,  1 = use DARS                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using using DARS is the default                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* we use a bit in variant */
  if (usedars)
    par->variant |= TABL_VARFLAG_USEDARS;
  else
    par->variant &= ~TABL_VARFLAG_USEDARS;

  /* changelog */
  par->set |= TABL_SET_USE_DARS;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_set_usedars() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_darsfactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for derandomized adaptive rejection sampling              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... parameter for DARS                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
    
  /* store date */
  PAR->darsfactor = factor;

  /* changelog */
  par->set |= TABL_SET_DARS_FACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_darsfactor() */

/*---------------------------------------------------------------------------*/

int 
unur_tabl_set_variant_splitmode( struct unur_par *par, unsigned splitmode )
     /*----------------------------------------------------------------------*/
     /* set setup variant for adaptive rejection sampling                    */
     /*                                                                      */
     /* There are three variants for adaptive rejection sampling. These      */
     /* differ in the way how an interval is split:                          */
     /*    splitmode 1: use the generated point to split the interval.       */
     /*    splitmode 2: use the mean point of the interval.                  */
     /*    splitmode 3: use the arcmean point.                               */
     /* Default is splitmode 3.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   splitmode ... indicator for variant                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* store date */
  par->variant &= ~TABL_VARMASK_SPLIT;
  switch (splitmode) {
  case 1:
    par->variant |= TABL_VARFLAG_SPLIT_POINT;
    return UNUR_SUCCESS;
  case 2:
    par->variant |= TABL_VARFLAG_SPLIT_MEAN;
    return UNUR_SUCCESS;
  case 3:
    par->variant |= TABL_VARFLAG_SPLIT_ARC;
    return UNUR_SUCCESS;
  default:
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid variant");
    return UNUR_ERR_PAR_SET;
  }
} /* end if unur_tabl_set_variant_splitmode() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_max_sqhratio( struct unur_par *par, double max_ratio )
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
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ratio = max_ratio;

  /* changelog */
  par->set |= TABL_SET_MAX_SQHRATIO;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_sqhratio( const struct unur_gen *gen )
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
  _unur_check_gen_object( gen, TABL, INFINITY );

  return (GEN->Asqueeze / GEN->Atotal);

} /* end of unur_tabl_get_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_hatarea( const struct unur_gen *gen )
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
  _unur_check_gen_object( gen, TABL, INFINITY );

  return GEN->Atotal;

} /* end of unur_tabl_get_hatarea() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_squeezearea( const struct unur_gen *gen )
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
  _unur_check_gen_object( gen, TABL, INFINITY );

  return GEN->Asqueeze;

} /* end of unur_tabl_get_squeezearea() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_max_intervals( struct unur_par *par, int max_ivs )
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
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= TABL_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_get_n_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get current number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   number of intervals     ... on success                             */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, TABL, 0 );

  return GEN->n_ivs;

} /* end of unur_tabl_get_n_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_areafraction( struct unur_par *par, double fraction )
     /*----------------------------------------------------------------------*/
     /* set parameter for equal area rule                                    */
     /* (each bar has size fraction * area below PDF)                        */           
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   fraction  ... fraction of area for bar                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (fraction < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"area factor < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->area_fract = fraction;

  /* changelog */
  par->set |= TABL_SET_AREAFRACTION;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_areafraction() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_nstp( struct unur_par *par, int n_stp )
     /*----------------------------------------------------------------------*/
     /* set number of construction points for hat at initialization          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TABL_SET_N_STP;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_nstp() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_slopes( struct unur_par *par, const double *slopes, int n_slopes )
     /*----------------------------------------------------------------------*/
     /* set slopes of PDF                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   slopes   ... pointer to list of slopes                             */
     /*   n_slopes ... number of slopes                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a]                        */
     /*   such that PDF(a) >= PDF(b).                                        */
     /*   slopes must be decreasing, non-overlapping and sorted              */
     /*----------------------------------------------------------------------*/
{
  int i;
  double al, bl;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if( n_slopes <= 0 ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"number of slopes <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* check slopes */
  al = slopes[0];
  bl = slopes[1];
  for( i=1; i<n_slopes; i++ ) {
    /* we do not check here if f(a) >= f(b), since we make no calculations heres */
    if( al > slopes[2*i] || bl > slopes[2*i+1] ) {
      _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes (overlapping or not in ascending order)");
      return UNUR_ERR_PAR_SET;
    }
    al = slopes[2*i];
    bl = slopes[2*i+1];
  }

  /* INFINITY is not allowed */
  if (_unur_FP_is_minus_infinity(slopes[0]) || _unur_FP_is_infinity(slopes[2*n_slopes-1])) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes must be bounded");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->slopes = slopes;
  PAR->n_slopes = n_slopes;

  /* changelog */
  par->set |= TABL_SET_SLOPES;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_slopes() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_guidefactor( struct unur_par *par, double factor )
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
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= TABL_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_boundary( struct unur_par *par, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set left and right boundary of computation interval                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   new boundary points must not be +/- INFINITY                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (left <= -INFINITY || right >= INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->bleft = left;
  PAR->bright = right;

  /* changelog */
  par->set |= TABL_SET_BOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_boundary() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, TABL );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TABL_VARFLAG_VERIFY) : (par->variant & (~TABL_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, TABL, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= TABL_VARFLAG_VERIFY;
    SAMPLE = _unur_tabl_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~TABL_VARFLAG_VERIFY;
    SAMPLE = _unur_tabl_sample;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_tabl_init( struct unur_par *par )
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
  int i,k;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_TABL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tabl_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* get slopes for starting generator */
  if (_unur_tabl_get_starting_intervals(par,gen)!=UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make hat function.");
    _unur_par_free(par); _unur_tabl_free(gen);
    return NULL;
  }

  if (par->variant & TABL_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */
 
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TABL_DEBUG_DARS) {
      /* make initial guide table (only necessary for writing debug info) */
      _unur_tabl_make_guide_table(gen);
      /* write info into log file */
      _unur_tabl_debug_init(par,gen);
      _unur_tabl_debug_dars_start(par,gen);
    }
#endif

    for (i=0; i<3; i++) {
      /* we make several tries */
      
      /* run DARS */
      if (_unur_tabl_run_dars(par,gen)!=UNUR_SUCCESS) {
	_unur_par_free(par); _unur_tabl_free(gen);
	return NULL;
      }

      /* make initial guide table */
      _unur_tabl_make_guide_table(gen);
  
      /* check if DARS was completed */
      if (GEN->n_ivs < GEN->max_ivs) {
	/* ran ARS instead */
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
      if (gen->debug) {
        if (gen->debug & TABL_DEBUG_DARS)
	  _unur_tabl_debug_dars(par,gen);
	else
	  _unur_tabl_debug_init(par,gen);
	_unur_tabl_debug_init_finished(par,gen);
      }
#endif
  }

  else { /* do not run DARS */
    /* make initial guide table */
    _unur_tabl_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) {
      _unur_tabl_debug_init(par,gen);
      _unur_tabl_debug_init_finished(par,gen);
    }
#endif
  }

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_tabl_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_tabl_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_tabl_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TABL_GEN);

  /* check for required data: area */
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA))
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS)
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF, use default instead");

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & TABL_VARFLAG_VERIFY) ? _unur_tabl_sample_check : _unur_tabl_sample;
  gen->destroy = _unur_tabl_free;
  gen->clone = _unur_tabl_clone;

  /* set all pointers to NULL */
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;

  /* the boundaries for our computation limits are intersection of the       */
  /* domain of the distribution and the given computation boundaries.        */
  if (par->distr->set & UNUR_DISTR_SET_DOMAIN) {
    PAR->bleft  = max(PAR->bleft, DISTR.BD_LEFT);
    PAR->bright = min(PAR->bright,DISTR.BD_RIGHT);
  }
  GEN->bleft       = PAR->bleft;         /* left boundary of domain            */
  GEN->bright      = PAR->bright;        /* right boundary of domain           */

  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN->max_ivs   = PAR->max_ivs;         /* maximum number of intervals        */
  GEN->max_ratio = PAR->max_ratio;       /* bound for ratio  Atotal / Asqueeze */
  GEN->darsfactor = PAR->darsfactor;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_tabl_create() */

/*****************************************************************************/

struct unur_gen *
_unur_tabl_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_tabl_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_tabl_interval *iv,*next, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_tabl_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tabl_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
    }
    /* next step */
    next = iv->next;
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;

  /* make new guide table */
  CLONE->guide = NULL;
  _unur_tabl_make_guide_table(clone);

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tabl_clone() */

/*****************************************************************************/

double
_unur_tabl_sample( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,INFINITY);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen->urng);

    /* look up in guide table and search for interval */
    iv =  GEN->guide[(int) (u * GEN->guide_size)];
    u *= GEN->Atotal;

    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,INFINITY);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u < iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
    }

    else {
      /* between spueeze and hat --> have to valuate PDF */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(x);

      /* being above squeeze is bad. split interval. */
      if (GEN->n_ivs < GEN->max_ivs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  switch (_unur_tabl_split_interval( gen, iv, x, fx,(gen->variant & TABL_VARMASK_SPLIT)) ) {
	  case UNUR_SUCCESS:
	  case UNUR_ERR_SILENT:
	    _unur_tabl_make_guide_table(gen);
	    /** TODO: it is not necessary to update the guide table every time. 
		But then (1) some additional bookkeeping is required and
		(2) the guide table method requires a acc./rej. step. **/
	    break;
	  default:
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  }
	}
	else {
	  /* no more construction points (avoid to many second if statement above */
	  GEN->max_ivs = GEN->n_ivs;
	}
      }

      /* now accept or reject */
      u = _unur_call_urng(gen->urng);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	return x;
    }
  }

} /* end of _unur_tabl_sample() */

/*****************************************************************************/

double
_unur_tabl_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TABL_GEN,INFINITY);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen->urng);

    /* look up in guide table and search for interval */
    iv =  GEN->guide[(int) (u * GEN->guide_size)];
    u *= GEN->Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,INFINITY);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      x = iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze;
      /* test whether PDF is monotone */
      fx = PDF(x);
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");
      /* at last return number */
      return x;
    }

    else {
      /* between spueeze and hat --> have to valuate PDF */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      fx = PDF(x);

      /* test whether PDF is monotone */
      if (_unur_FP_greater(fx,iv->fmax))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. PDF not monotone in interval");
      if (_unur_FP_less(fx,iv->fmin))
	_unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. PDF not monotone in interval");

      /* being above squeeze is bad. split interval. */
      if (GEN->n_ivs < GEN->max_ivs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  switch (_unur_tabl_split_interval( gen, iv, x, fx,(gen->variant & TABL_VARMASK_SPLIT)) ) {
	  case UNUR_SUCCESS:
	  case UNUR_ERR_SILENT:
	    _unur_tabl_make_guide_table(gen);
	    /** TODO: it is not necessary to update the guide table every time. 
		But then (1) some additional bookkeeping is required and
		(2) the guide table method requires a acc./rej. step. **/
	    break;
	  default:
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  }
	}
	else {
	  /* no more construction points (avoid to many second if statement above */
	  GEN->max_ivs = GEN->n_ivs;
	}
      }

      /* now accept or reject */
      u = _unur_call_urng(gen->urng);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	return x;
    }
  }

} /* end of _unur_tabl_sample_check() */

/*****************************************************************************/

void
_unur_tabl_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TABL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tabl_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_tabl_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free table */
  if (GEN->guide)  free(GEN->guide);

  /* free other memory */
  _unur_generic_free(gen);

} /* end of _unur_tabl_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static int
_unur_tabl_get_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   PDF(a) >= PDF(b)                                                   */
     /*----------------------------------------------------------------------*/
{

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* we have two cases: 
     (1) we are given slopes --> check these, compute domain if necessary
     (2) we are given domain and mode --> compute slopes */

  if (PAR->n_slopes > 0 )
    /* slopes are given */
    return _unur_tabl_get_starting_intervals_from_slopes(par,gen);

  if (par->distr->set & UNUR_DISTR_SET_MODE)
    /* no slopes given. need domain and mode */
    /* compute slopes */
    return _unur_tabl_get_starting_intervals_from_mode(par,gen);

  /* else */
  _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"number of slopes <= 0, domain or mode not given.");
  return UNUR_ERR_GEN_DATA;

} /* end of _unur_tabl_get_starting_intervals() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_get_starting_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals, slopes are given by user.                */
     /* estimate domain when not given.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   PDF(a) >= PDF(b)                                                   */
     /*----------------------------------------------------------------------*/
{
  /** TODO: check for slopes out of support !! **/

  struct unur_tabl_interval *iv;
  int i;

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* init counter of intervals */
  GEN->n_ivs = 0;
  iv = NULL;

  /* boundary of computational interval are reset by boundaries of slopes */
  GEN->bleft = INFINITY;
  GEN->bright = -INFINITY;

  /* compute initial intervals */
  for ( i=0; i < 2*PAR->n_slopes; i+=2 ) {
    /* get a new interval and link into list */
    if (i==0)
      iv = GEN->iv = _unur_xmalloc(sizeof(struct unur_tabl_interval));    /* the first interval */
    else
      iv = iv->next = _unur_xmalloc(sizeof(struct unur_tabl_interval));  /* all the other intervals */
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);

    /* max and min of PDF in interval */
    iv->xmax = PAR->slopes[i];      
    iv->fmax = PDF(iv->xmax);
    iv->xmin = PAR->slopes[i+1];    
    iv->fmin = PDF(iv->xmin);

    /* check for overflow */
    if (_unur_FP_is_infinity(iv->fmax) || _unur_FP_is_infinity(iv->fmin)) {
      /* overflow */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      iv->next = NULL;  /* terminate list (freeing list) */
      return UNUR_ERR_GEN_DATA;
    }

    /* check slopes */
    if (iv->fmax < iv->fmin) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"slopes non-decreasing");
      iv->next = NULL;  /* terminate list (freeing list) */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area of slope and sign of slope (increasing/decreasing) */
    iv->slope = (iv->xmax > iv->xmin) ? 1 : -1;
    iv->Ahat = iv->slope * (iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = iv->slope * (iv->xmax - iv->xmin) * iv->fmin;
    /** TODO: possible overflow/underflow ?? **/
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* estimate domain */
    if (iv->slope > 0) {
      GEN->bleft = min(GEN->bleft,iv->xmin);
      GEN->bright = max(GEN->bright,iv->xmax);
    }
    else {
      GEN->bleft = min(GEN->bleft,iv->xmax);
      GEN->bright = max(GEN->bright,iv->xmin);
    }

    /* split interval following [1], split A */
    if (par->variant & TABL_VARFLAG_STP_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return UNUR_ERR_GEN_DATA;
    }
  }

  /* terminate list */
  iv->next = NULL;

  /* reset domain of distribution */
  DISTR.trunc[0] = DISTR.BD_LEFT = GEN->bleft;
  DISTR.trunc[1] = DISTR.BD_RIGHT = GEN->bright;

  /* reset area below distribution */
  gen->distr->set &= ~UNUR_DISTR_SET_PDFAREA;
  unur_distr_cont_upd_pdfarea( gen->distr );

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_get_starting_intervals_from_slopes() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_get_starting_intervals_from_mode( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals                                           */
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
  /** TODO: check for slopes out of support !! **/

  struct unur_tabl_interval *iv;

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* init linked list of intervals */
  GEN->n_ivs = 0;

  /* compute initial intervals */
  while (1) {
    /* the first interval */
    iv = GEN->iv = _unur_xmalloc(sizeof(struct unur_tabl_interval));
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);

    if (DISTR.mode <= GEN->bleft) {
      /* only one ascending interval <a,b> = [a,b] */
      iv->xmax = GEN->bleft;
      iv->xmin = GEN->bright;
      break;
    }

    if (DISTR.mode >= GEN->bright) {
      /* only one descending interval <a,b> = [b,a] */
      iv->xmax = GEN->bright;
      iv->xmin = GEN->bleft;
      break;
    }

    /* one descending and one ascending interval */
    iv->xmax = DISTR.mode;
    iv->xmin = GEN->bleft;

    /* the second interval */
    iv = iv->next = _unur_xmalloc(sizeof(struct unur_tabl_interval));  /* all the other intervals */
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);

    iv->xmax = DISTR.mode;
    iv->xmin = GEN->bright;
    break;
  }

  /* terminate list */
  iv->next = NULL;

  /* compute parameters */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);

    /* max and min of PDF in interval */
    iv->fmax = PDF(iv->xmax);
    iv->fmin = PDF(iv->xmin);

    /* check for overflow */
    if (_unur_FP_is_infinity(iv->fmax) || _unur_FP_is_infinity(iv->fmin)) {
      /* overflow */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      return UNUR_ERR_GEN_DATA;
    }

    /* area of slope and sign of slope (increasing/decreasing) */
    iv->slope = (iv->xmax > iv->xmin) ? 1 : -1;
    iv->Ahat = iv->slope * (iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = iv->slope * (iv->xmax - iv->xmin) * iv->fmin;
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* split interval following [1], split A */
    if (par->variant & TABL_VARFLAG_STP_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return UNUR_ERR_GEN_DATA;
    }

  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_get_starting_intervals_from_mode() */

/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_split_a_starting_intervals( struct unur_par *par, 
				       struct unur_gen *gen, 
				       struct unur_tabl_interval *iv_slope )
     /*----------------------------------------------------------------------*/
     /* split starting intervals according to [1]                            */
     /* SPLIT A (equal areas rule)                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter list                            */
     /*   gen       ... pointer to generator object                          */
     /*   iv_slope  ... pointer to interval of slope                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to last interval in list of splitted slope                 */
     /*   NULL on error                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv, *iv_last;
  double bar_area, x, fx;
  
  /* check arguments */
  CHECK_NULL(par,NULL);       COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  CHECK_NULL(gen,NULL);       COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  CHECK_NULL(iv_slope,NULL);  COOKIE_CHECK(iv_slope,CK_TABL_IV,NULL);

  if (GEN->n_ivs >= GEN->max_ivs) 
    /* maximal number of intervals reached */
    return NULL;

  if (iv_slope->slope != 1 && iv_slope->slope != -1 ) {
    /* this should not happen:
       invalid slope.          */
    _unur_warning( gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return iv_slope;
  }

  iv = iv_slope;        /* pointer to actual interval */
  iv_last = iv_slope;   /* pointer to last interval in list */
  /* (maximal) area of bar (= hat in one interval) */
  bar_area = DISTR.area * PAR->area_fract;

  while (iv->Ahat > bar_area) {
    /* compute splitting point:
       slope == +1 --> move from right to left
       slope == -1 --> move from left to right */
    x = iv->xmax - iv->slope * bar_area / iv->fmax;

    /* compute PDF */
    fx = PDF(x);

    /* check for overflow */
    if (_unur_FP_is_infinity(fx)) {
      /* overflow */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      return NULL;
    }

    /* now split interval at x */
    switch (_unur_tabl_split_interval( gen, iv, x, fx, TABL_VARFLAG_SPLIT_POINT )) {
    case UNUR_SUCCESS:  /* splitting succesful */
      if (iv->slope > 0) {
	if (iv_last == iv_slope)
	  iv_last = iv->next;
      }
      else { /* iv->slope < 0 */
	iv = iv->next; break;
      }
      break;
    case UNUR_ERR_SILENT: /* interval chopped */
      break; /* nothing to do */
    default: /* error (slope not monotonically increasing) */
      return NULL;
    }

    /* check number of intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* maximal number of intervals reached */
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"split A stopped, maximal number of intervals reached.");
      break;
    }
  }

  /* pointer to last interval */
  return ((iv->slope < 0) ? iv : iv_last);

} /* end of _unur_tabl_split_a_starting_intervals() */

/*---------------------------------------------------------------------------*/

int
_unur_tabl_run_dars( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling.                         */
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
  struct unur_tabl_interval *iv;
  double Atot, Asqueezetot;    /* total area below hat and squeeze, resp. */
  double Alimit;               /* threshhold value for splitting interval */
  int n_splitted = 1;          /* count splitted intervals */

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* there is no need to run DARS when the DARS factor is INFINITY */
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_arou_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;

  /* now split intervals */
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_ivs < GEN->max_ivs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN->n_ivs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_ivs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 

    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);

      /* do not exceed the maximum number of intervals */
      if (GEN->n_ivs >= GEN->max_ivs)
	break;

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  /* goto next interval */

      switch (_unur_tabl_split_interval( gen, iv, 0., 0., TABL_VARFLAG_SPLIT_ARC )) {
      case UNUR_SUCCESS:  /* splitting succesful */
      case UNUR_ERR_SILENT: /* interval chopped */
	++n_splitted;
	break; /* nothing to do */
      default: /* error (slope not monotonically decreasing) */
	return UNUR_ERR_GEN_DATA;
      }
    }

    if (n_splitted == 0) {
      /* we are not successful in splitting any inteval.
	 abort to avoid endless loop */
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }

  /* ratio between squeeze and hat o.k. ? */
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_ivs >= GEN->max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN->max_ivs = GEN->n_ivs;
  }
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_run_dars() */

/*****************************************************************************/

static int
_unur_tabl_split_interval( struct unur_gen *gen,
			   struct unur_tabl_interval *iv_old, 
			   double x, double fx, 
			   unsigned split_mode )
     /*----------------------------------------------------------------------*/
     /* split interval (replace old one by two new ones in same place)       */
     /* new interval is inserted immedately after old one.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   iv_old     ... pointer to interval that has to be split            */
     /*   x          ... splitting point                                     */
     /*   fx         ... value of PDF at splitting point                     */
     /*   split_mode ... how to split interval                               */
     /*                  TABL_VARFLAG_SPLIT_POINT: split at given point x    */
     /*                  TABL_VARFLAG_SPLIT_MEAN:  at mean point of interval */
     /*                  TABL_VARFLAG_SPLIT_ARC:   at arc mean point         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... splitting successful                           */
     /*   UNUR_ERR_SILENT ... interval chopped off domain                    */
     /*                       (not part of support of PDF)                   */
     /*   others          ... error: PDF not monotone in interval            */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv_new;
  double A_hat_old, A_squ_old;
      
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);     COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_old,UNUR_ERR_NULL);  COOKIE_CHECK(iv_old,CK_TABL_IV,UNUR_ERR_COOKIE);

  /* There are three possibilities for the splitting point:
     (1) use x and avoid computation of PDF(x). 
     (2) use middle of interval. converges faster in many cases.
     (3) use "arc_mean" of interval. 
         converges faster when domain is almost unbounded. */
  switch( split_mode ) {
  case TABL_VARFLAG_SPLIT_POINT:    /* (1) */
    /* nothing to do (default) */
    break;
  case TABL_VARFLAG_SPLIT_MEAN:     /* (2) */
    x = 0.5 * (iv_old->xmin + iv_old->xmax); 
    fx = PDF(x);
    break;
  case TABL_VARFLAG_SPLIT_ARC:      /* (3) */
    x = _unur_arcmean(iv_old->xmin, iv_old->xmax); 
    fx = PDF(x);
    break;
  default: 
    /* this should not happen:
       Invalid variant, use default n*/
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    break;
  }

  /* check for overflow */
  if (_unur_FP_is_infinity(fx)) {
    /* overflow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return UNUR_ERR_GEN_DATA;
  }

  /* store areas of old interval */
  A_hat_old = iv_old->Ahat;
  A_squ_old = iv_old->Asqueeze;

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {
    /* check montonicity */
    if (iv_old->fmin > 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not monotone in slope");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* chop off part out of support */
    iv_old->xmin = x;

    /* compute new area in interval */
    /** TODO: possible overflow/underflow ?? **/
    iv_old->Ahat = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmax;
    /* iv_old->Asqueeze remains 0 */

    /* update total area */
    GEN->Atotal += iv_old->Ahat - A_hat_old;
    /* GEN->Asqueeze remains unchanged */

    /* interval chopped but not split */
    return UNUR_ERR_SILENT;
  }

  /* we need a new interval */
  iv_new = _unur_xmalloc(sizeof(struct unur_tabl_interval));
  ++(GEN->n_ivs);
  COOKIE_SET(iv_new,CK_TABL_IV);

  /* iv_new has the same slope as iv_old */
  iv_new->slope = iv_old->slope;

  /* we have to distinguish between two cases:
     PDF is increasing (slope = +1) or
     PDF is decreasing (slope = -1). */

  switch (iv_old->slope) {
  case -1:
    /* (x) The iv_new inherits the minimum of iv_old.
           iv_old keeps the maximum.
       (x) The splitting point is the maximum of iv_new and
           the minimum of iv_old.
    */
    iv_new->xmin  = iv_old->xmin;  
    iv_new->fmin = iv_old->fmin;
    iv_old->xmin  = iv_new->xmax = x; 
    iv_old->fmin = iv_new->fmax = fx; 
    break;
  case +1: /* the other way round */
    /* (x) The iv_new inherits the maximum of iv_old.
           iv_old keeps the minimum.
       (x) The splitting point is the minimum of iv_new and
           the maximum of iv_old.
    */
    iv_new->xmax  = iv_old->xmax;  
    iv_new->fmax = iv_old->fmax;
    iv_old->xmax  = iv_new->xmin = x; 
    iv_old->fmax = iv_new->fmin = fx; 
    break;
  default: 
    /* this should not happen:
       Invalid slope. Cannot split interval. */
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* compute the areas in both intervals */
  /** TODO: possible overflow/underflow ?? **/
  iv_new->Ahat     = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmax;
  iv_new->Asqueeze = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmin;
  iv_old->Ahat     = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmax;
  iv_old->Asqueeze = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmin;

  /* update total areas */
  GEN->Atotal += iv_old->Ahat + iv_new->Ahat - A_hat_old;
  GEN->Asqueeze += iv_old->Asqueeze + iv_new->Asqueeze - A_squ_old;

  /* insert iv_new into linked list of intervals.
     iv_old is stored on the left hand side of iv_new. */
  iv_new->next = iv_old->next;
  iv_old->next = iv_new;

  /* splitting successful */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_split_interval() */

/*****************************************************************************/

static int
_unur_tabl_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*   1 (--> successful)                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    int max_guide_size = (GEN->guide_factor > 0.) ? (GEN->max_ivs * GEN->guide_factor) : 1;
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tabl_interval*) );
  }

  /* first we need the cumulated areas of rectangles */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
    
  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN->guide_size = GEN->n_ivs;

  /* make table (use variant 2; see dis.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      if( iv->next != NULL )    /* skip to next segment if it exists */
        iv = iv->next;
      else {
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
} /* end of _unur_tabl_make_guide_table() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(log,"%s:\n",gen->genid);

/*---------------------------------------------------------------------------*/

void
_unur_tabl_debug_init( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  empty_line();

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tabl_sample",gen->genid);
  if (par->variant & TABL_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  empty_line();

  fprintf(log,"%s: computation interval = (%g, %g)\n",gen->genid,GEN->bleft,GEN->bright);
  empty_line();

  fprintf(log,"%s: area fraction for equal area rule = %g ",gen->genid,PAR->area_fract);
  _unur_print_if_default(par,TABL_SET_AREAFRACTION);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR->max_ivs);
  _unur_print_if_default(par,TABL_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR->max_ratio*100.);
  _unur_print_if_default(par,TABL_SET_MAX_SQHRATIO);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR->guide_factor);
  _unur_print_if_default(par,TABL_SET_GUIDEFACTOR);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: split intervals at ",gen->genid);
  switch( gen->variant & TABL_VARMASK_SPLIT ) {
  case TABL_VARFLAG_SPLIT_MEAN:
    fprintf(log,"mean point");
    break;
  case TABL_VARFLAG_SPLIT_ARC:
    fprintf(log,"\"arcmean\" point");
    break;
  case TABL_VARFLAG_SPLIT_POINT:
  default: 
    fprintf(log,"sample point");
    break;
  }
  fprintf(log," when using adaptive sampling.\n");
  empty_line();

  if (par->set & TABL_SET_SLOPES) {
    fprintf(log,"%s: slopes = %d\n",gen->genid,PAR->n_slopes);
    for (i=0; i<PAR->n_slopes; i++) {
      if ( PAR->slopes[2*i] > PAR->slopes[2*i+1] )
	fprintf(log,"%s:   (+)  ",gen->genid);
      else
	fprintf(log,"%s:   (-)  ",gen->genid);
      fprintf(log,"< %#g, %#g >\n", PAR->slopes[2*i], PAR->slopes[2*i+1] );
    }
  }
  else
    fprintf(log,"%s: no slopes given. compute from domain and mode.\n",gen->genid);

  if (par->variant & TABL_VARFLAG_STP_A)
    fprintf(log,"%s: split slopes by equal area rule (SPLIT A).\n",gen->genid);

  if (par->variant & TABL_VARFLAG_USEDARS) {
    fprintf(log,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(par,TABL_SET_USE_DARS);
    fprintf(log,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  }
  else {
    fprintf(log,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(par,TABL_SET_USE_DARS);
  }
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting intervals (approx.) = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(par,TABL_SET_N_STP);
  fprintf(log,"\n");
  empty_line();

  if (gen->debug & TABL_DEBUG_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);

} /* end of _unur_tabl_debug_init() */

/*****************************************************************************/

void
_unur_tabl_debug_init_finished( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_tabl_debug_intervals(gen,TRUE);

  empty_line();
  fprintf(log,"%s:  INIT completed **********************\n",gen->genid);
  empty_line();

} /* end of _unur_tabl_debug_init_finished() */

/*****************************************************************************/

void 
_unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print header before runniung DARS into logfile                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: DARS started **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tabl_debug_dars_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tabl_debug_dars( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print infor after generator has run DARS into logfile                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS finished **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  if (gen->debug & TABL_DEBUG_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tabl_debug_dars() */

/*****************************************************************************/

void
_unur_tabl_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  empty_line();
  _unur_tabl_debug_intervals(gen,TRUE);
  empty_line();

  fflush(log);

} /* end of _unur_tabl_debug_free() */

/*****************************************************************************/

void
_unur_tabl_debug_intervals( const struct unur_gen *gen, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tabl_interval *iv;
  double sAsqueeze, Atotal;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: intervals = %d\n",gen->genid,GEN->n_ivs);
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID);
      fprintf(log,"%s:[%3d]: (",gen->genid,i);
      switch (iv->slope) {
      case 1:  fprintf(log,"+"); break;
      case 0:  fprintf(log,"0"); break;
      case -1: fprintf(log,"-"); break;
      default:
	/* this should not happen:
	   invalid value for iv->slope. */
	fprintf(log,"?"); break;
      }
      fprintf(log,")   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }

  if (!print_areas) {
    /* this just happens, when we call _unur_tabl_debug_intervals(), after
       we have finished split A. */
    fprintf(log,"%s: Split A finished.\n",gen->genid);
    empty_line();
    return;
  }

  if (GEN->Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    Atotal = -1.;   /* to avoid floating point exceptions */
  }
  else {
    Atotal = GEN->Atotal;
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID); 
      sAsqueeze += iv->Asqueeze;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
	      iv->Ahat, iv->Ahat * 100. / Atotal, 
	      iv->Acum, iv->Acum * 100. / Atotal);
    }
    fprintf(log,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / Atotal, Atotal);
    empty_line();
  }
    
  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);

  empty_line();

} /* end of _unur_tabl_debug_intervals */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
