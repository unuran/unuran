/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    Polynomial interpolation based INVersion of CDF              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Compute values of CDF incrementally and interpolate resulting points *
 *      by polynomials.                                                      *
 *                                                                           *
 *      Integration:   adaptive Gauss-Lobatto integration with 5 points      *
 *      Interpolation: Newton recursion for interpolating polynomial         *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the PDF, center of distribution                           *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
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
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Method PINV combines numerical integration with interpolation of the    *
 *   inverse CDF.                                                            *
 *                                                                           *
 *   1.  Preprocessing:                                                      *
 *                                                                           *
 *   1a.   Estimate computationally relevant domain (support) of PDF         *
 *         (finite interval where PDF is above some threshold value).        *
 *            _unur_pinv_relevant_support()                                  *
 *            _unur_pinv_searchborder()                                      *
 *                                                                           *
 *   1b.   Compute area below PDF over relevant domain approximately.        *
 *            _unur_pinv_approx_pdfarea()                                    *
 *                                                                           *
 *   1c.   Compute computational domain where inverse CDF is approximated    *
 *         (interval where we safely can compute coefficients of             *
 *         interpolating polynomial).                                        *
 *            _unur_pinv_computational_domain()                              *
 *            _unur_pinv_cut()                                               *
 *            _unur_pinv_tailprob()                                          *
 *                                                                           *
 *   2.  Interpolation:                                                      *
 *                                                                           *
 *   2a.   Compute coefficients for interpolating polynomial for             *
 *         fixed (sub-) interval.                                            *
 *            _unur_pinv_newton_create()                                     *
 *            _unur_pinv_chebyshev_points()                                  *
 *                                                                           *
 *   2b.   Evaluate interpolating polynomial.                                *
 *            _unur_pinv_newton_eval()                                       *
 *                                                                           *
 *   2c.   Estimate approximation error for given interpolating polynomial.  *
 *            _unur_pinv_newton_maxerror()                                   *
 *            _unur_pinv_newton_testpoints()                                 *
 *                                                                           *
 *                                                                           *
 *   Currently the following methods are implemented:                        *
 *                                                                           *
 *   Quadrature (Integration):                                               *
 *     Adaptive Gauss-Lobatto integration with 5 points.                     *
 *        _unur_pinv_lobatto5()                                              *
 *                                                                           *
 *   Interpolation:                                                          *
 *     Newton recursion for coefficients of polynomial                       *
 *     ("Newton interpolation").                                             *
 *                                                                           *
 *****************************************************************************



TODO: Das ganze funktioniert bei exponentialverteilung (und normalverteilung 
auf (-1,1)) nicht!!!!!
Domain der verteilung benutzen!!!!!
vielleicht sollte man auch defaults fuer randsuche davon abhaengig machen,
welchen domain der benutzer angegeben hat (INFINITY?)



Wichtig: die dichte darf nicht zwischen einzelnen modi auf zu langen
intervallen zu klein (fast 0) werden.

ui sind die stuetzstellen der Interpolation von F^{-1} in jedem subinterval.
zi sind die dazugehoerigen polynomkoeffizienten der Newton Interpolation, also
Koeffizienten eines Horner-aehnlichen Schemas.

Die Berechnung der zi folgt der Reukursion der Newton Interpolation. Gerhard hat
das aus einem Deutschen Numerik Lehrbuch, ich hab es mit Wikipedia "Newton
Interpolation" verglichen. Es ist das gleiche.
xi sind die Intervallgrenzen der sub-intervalle.


 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include <tests/unuran_tests.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "pinv.h"
#include "pinv_struct.h"


/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* -- Global parameters                                                      */

#define MAX_ORDER   (19) 
/* maximum order of Newton interpolation polynomial */

#define PINV_UERROR_CORRECTION  (0.9)
/* PINV tries to create an approximation of the inverse CDF where the        */
/* U-error is bounded by the given u-resolution. However, the error caused   */
/* by cutting off the tails introduces some additional errors which must be  */
/* corrected. Thus the upper bound for the pure approximation error of the   */
/* interpolation is set to PINV_UERROR_CORRECTION * u-resolution.            */

#define PINV_MAX_IVS  (10000)
/* maximum number of intervals */

/* -- 1. Preprocessing                                                       */

#define PINV_PDFLLIM    (1.e-13)
/* Threshold value used for finding the boundary of the computational        */
/* domain. When starting the search at some point x0 then the search stops   */
/* when a point x is found where PDF(x) approx. PDF(x0) * PINV_PDFLLIM.      */

#define PINV_TAILCUTOFF_FACTOR(ures) ((ures <= 9.e-13) ? 0.5 : 0.1)
#define PINV_TAILCUTOFF_MAX          (1.e-10) 
/* For unbounded domains the tails has to be cut off. We use the given       */
/* u-resolution for finding the cut points. (The probability for each of the */
/* chopped regions should be less than                                       */
/* HINV_TAILCUTOFF_FACTOR * u-resolution.)                                   */
/* For very small u-resolution we need a higher factor for distributions     */
/* with heavy tails, e.g. the Cauchy distribution, since otherwise we have   */
/* numerical problems to find proper cut-off points .                        */
/*                                                                           */
/* However, the tail probabilities should not be greater than some threshold */
/* value, given by PINV_TAILCUTOFF_MAX which reflects the precision of the   */
/* used stream of uniform pseudo-random numbers (typically about 2^32).      */
/* However, for computational reasons we use a value that is at least twice  */
/* the machine epsilon for the right hand boundary.                          */

/* -- 2. Newton interpolation                                                */

#define PINV_MAX_ITER_IVS    (10 * PINV_MAX_IVS)
/* maximum number of iterations for computing intervals for Newtwon          */
/* interpolation polynomials. Obviously it should be larger than the maximum */
/* number of intervals (or the latter can be reduced).                       */

#define PINV_GUIDE_FACTOR  (1)
/* relative size of guide table for finding the subinterval corresponding    */
/* to the given U-value.                                                     */


/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define PINV_DEBUG_REINIT    0x00000002u   /* print parameters after reinit  */
#define PINV_DEBUG_TABLE     0x00000010u   /* print table                    */
#define PINV_DEBUG_SEARCHBD  0x00010000u   /* trace search boundary          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define PINV_SET_ORDER          0x001u  /* order of polynomial               */
#define PINV_SET_U_RESOLUTION   0x002u  /* maximal error in u                */
#define PINV_SET_BOUNDARY       0x008u  /* boundary of computational region  */
#define PINV_SET_SEARCHBOUNDARY 0x010u  /* search for boundary               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "PINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_init (struct unur_par *par);
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_reinit (struct unur_gen *gen); */
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_create (struct unur_par *par);
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_check_par (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_clone (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_free (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_sample (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_eval_approxinvcdf (const struct unur_gen *gen, double u);
/*---------------------------------------------------------------------------*/
/* evaluate Hermite interpolation of inverse CDF at u.                       */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_lobatto5 (struct unur_gen *gen, double x, double h, double tol);
/*---------------------------------------------------------------------------*/
/* numerical integration of the PDF over the interval (x,x+h)                */
/* using adaptive Gauss-Lobatto integration with 5 points.                   */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_relevant_support (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1a. Estimate computationally relevant domain (support) of PDF             */
/*     (finite interval where PDF is above some threshold value).            */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound, double *dom);
/*---------------------------------------------------------------------------*/
/* [1a.] find left or right hand border of relevant domain.                  */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_approx_pdfarea (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1b. Compute area below PDF over relevant domain approximately.            */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_computational_domain (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1c. Compute computational domain where inverse CDF is approximated        */
/*     (interval where we safely can compute coefficients of                 */
/*     interpolating polynomial).                                            */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_cut (struct unur_gen *gen, double dom, double w, double dw, double crit);
/*---------------------------------------------------------------------------*/
/* [1c.] calculate cut-off points for computational domain of distribution.  */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_tailprob (struct unur_gen *gen, double x, double dx);
/*---------------------------------------------------------------------------*/
/* [1c.] calculate approximate tail probability.                             */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create table for Newton interpolation                                     */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx );
/*---------------------------------------------------------------------------*/
/* make a new interval i with left boundary point x and CDF(x).              */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
				     double *xval, double uerrcrit);
/*---------------------------------------------------------------------------*/
/* 2a. Compute coefficients for Newton interpolation.                        */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_chebyshev_points (int order, double *pt);
/*---------------------------------------------------------------------------*/
/* [2a.] Compute Chebyshev points.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_newton_eval (double q, double ui[], double zi[], int order);
/*---------------------------------------------------------------------------*/
/* 2b. evaluate Newton polynomial.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv, 
					  double xval[], double uerrcrit);
/*---------------------------------------------------------------------------*/
/* 2c. estimate maximal error of Newton interpolation in subinterval         */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_newton_testpoints (int g,double ui[],double utest[]);
/*---------------------------------------------------------------------------*/
/* [2c.] calculate the local maxima of the interpolation polynomial          */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_make_guide_table (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_init_start (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* print before setup starts.                                                */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_searchbd (const struct unur_gen *gen, int aftercut);
/*---------------------------------------------------------------------------*/
/* print computational before or after searching for cut-off points          */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_init (const struct unur_gen *gen, int ok);
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_intervals( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into logfile.               */
/*---------------------------------------------------------------------------*/

#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_pinv_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_pinv_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_pinv_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)  (_unur_cont_PDF((x),(gen->distr)))    /* call to PDF         */
/* #define dPDF(x) (_unur_cont_dPDF((x),(gen->distr)))   /\* call to derivative of PDF *\/    */

/*---------------------------------------------------------------------------*/

#define _unur_pinv_getSAMPLE(gen)  (_unur_pinv_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_pinv_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_pinv_par) );
  COOKIE_SET(par,CK_PINV_PAR);

  /* copy input */
  par->distr   = distr;           /* pointer to distribution object          */

  /* set default values */
  PAR->order = 5;                /* order of polynomial                      */
  PAR->u_resolution = 1.0e-10;   /* maximal error allowed in u-direction     */
  PAR->bleft = -1.e100;          /* left border of the computational domain  */
  PAR->bright = 1.e100;          /* right border of the computational domain */
  PAR->sleft = TRUE;             /* whether to search for left boundary      */
  PAR->sright = TRUE;            /* whether to search for right boundary     */

  par->method   = UNUR_METH_PINV; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_pinv_init;

  return par;

} /* end of unur_pinv_new() */

/*****************************************************************************/

int
unur_pinv_set_order( struct unur_par *par, int order)
     /*----------------------------------------------------------------------*/
     /* Set order of Hermite interpolation.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   order  ... order of interpolation polynome                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (order<2 || order>MAX_ORDER) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"order <2 or >19");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->order = order;

  /* changelog */
  par->set |= PINV_SET_ORDER;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_order() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... maximal error in u                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (u_resolution > 1.e-2) {
    /* this is obviously an error */
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too large --> not changed");
    return UNUR_ERR_PAR_SET;
  }
  if (u_resolution < 5.*DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small --> corrected");
    u_resolution = 5.*DBL_EPSILON;
  }

  /* store date */
  PAR->u_resolution = u_resolution;

  /* changelog */
  par->set |= PINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_boundary( struct unur_par *par, double left, double right )
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
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (!_unur_FP_less(left,right)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (! (_unur_isfinite(left) && _unur_isfinite(right)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->bleft = left;
  PAR->bright = right;

  /* changelog */
  par->set |= PINV_SET_BOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_boundary() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_searchboundary( struct unur_par *par, int left, int right )
     /*----------------------------------------------------------------------*/
     /* set flag for boundary searching algorithm                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   left  ... whether to search for left boundary point                */
     /*   right ... whether to search for right boundary point               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* store date */
  PAR->sleft  = (left)  ? TRUE : FALSE;
  PAR->sright = (right) ? TRUE : FALSE;

  /* changelog */
  par->set |= PINV_SET_SEARCHBOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_searchboundary() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_get_n_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get number of intervals (or more precisely the number of nodes)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   number of intervals ... on success                                 */
     /*   0     ... on error                                                 */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, PINV, 0 );
  return GEN->n_ivs;
} /* end of unur_pinv_get_n_intervals() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_pinv_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_PINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_PINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_pinv_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_pinv_check_par(gen) != UNUR_SUCCESS) {
    _unur_pinv_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_pinv_debug_init_start(gen);
#endif

  /* Preprocessing:                                     */
  /*   find interval for computing Newton interpolation */
  if (_unur_pinv_relevant_support(gen)     != UNUR_SUCCESS ||
      _unur_pinv_approx_pdfarea(gen)       != UNUR_SUCCESS ||
      _unur_pinv_computational_domain(gen) != UNUR_SUCCESS) {

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }

  /* compute table for Newton interpolation */
  if (_unur_pinv_create_table(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_pinv_debug_init(gen,FALSE);
#endif
    _unur_pinv_free(gen); return NULL;
  }

  /* make guide table */
  _unur_pinv_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_pinv_debug_init(gen,TRUE);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_pinv_init() */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_reinit( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* re-initialize (existing) generator.                                  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int rcode; */

/*   /\* check parameters *\/ */
/*   if ( (rcode = _unur_pinv_check_par(gen)) != UNUR_SUCCESS) */
/*     return rcode; */

/*   return UNUR_SUCCESS; */
/* } /\* end of _unur_pinv_reinit() *\/ */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_pinv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_PINV_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_pinv_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_PINV_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_pinv_getSAMPLE(gen);
  gen->destroy = _unur_pinv_free;
  gen->clone = _unur_pinv_clone;
  /* gen->reinit = _unur_pinv_reinit; */

  /* copy parameters into generator object */
  GEN->order = PAR->order;            /* order of polynomial                 */
  GEN->u_resolution = PAR->u_resolution; /* maximal error in u-direction     */
  GEN->bleft_par  = PAR->bleft;          /* border of computational domain   */
  GEN->bright_par = PAR->bright;
  GEN->sleft  = PAR->sleft;              /* whether to search for boundary   */
  GEN->sright = PAR->sright;

  /* initialize variables */
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->dleft = -INFINITY;
  GEN->dright = INFINITY;
  GEN->Umax = 1.;
  GEN->iv = NULL;
  GEN->n_ivs = -1;        /* -1 indicates that there are no intervals at all */
  GEN->guide_size = 0; 
  GEN->guide = NULL;
  GEN->area = 1.;         /* we use 1 as first guess */

  /* allocate maximal array of intervals */
  /* [ Maybe we could move this into _unur_pinv_interval() ] */
  GEN->iv =  _unur_xmalloc(PINV_MAX_IVS * sizeof(struct unur_pinv_interval) );

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_pinv_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_pinv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_check_par( struct unur_gen *gen )
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
  /* points for searching computational domain */
  GEN->bleft = _unur_max(GEN->bleft_par,DISTR.domain[0]);
  GEN->bright = _unur_min(GEN->bright_par,DISTR.domain[1]);

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* domain of distribution (used when x with PDF(x)=0 are found) */
  GEN->dleft =  DISTR.domain[0];
  GEN->dright =  DISTR.domain[1];

  /* center of distribution */
  DISTR.center = unur_distr_cont_get_center(gen->distr);
  if (PDF(DISTR.center)<=0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"PDF(center) <= 0.");
    return UNUR_ERR_GEN_CONDITION;
  }

  /** TODO: check that center is in given domain !! **/

  return UNUR_SUCCESS;
} /* end of _unur_pinv_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_pinv_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_pinv_gen*)clone->datap)

  struct unur_gen *clone;
  int i;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy coefficients for Newton polynomial */
  CLONE->iv =  _unur_xmalloc((GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  memcpy( CLONE->iv, GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );

  for(i=0; i<=GEN->n_ivs; i++) {
    CLONE->iv[i].ui = _unur_xmalloc( GEN->order * sizeof(double) );
    CLONE->iv[i].zi = _unur_xmalloc( GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].ui, GEN->iv[i].ui, GEN->order * sizeof(double) );
    memcpy( CLONE->iv[i].zi, GEN->iv[i].zi, GEN->order * sizeof(double) );
  }

  /* copy guide table */
  CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_pinv_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  int i;

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free guide table */
  if (GEN->guide) free (GEN->guide);

  /* free tables of coefficients of interpolating polynomials */
  if (GEN->iv) {
    for(i=0; i<=GEN->n_ivs; i++){
      free(GEN->iv[i].ui);
      free(GEN->iv[i].zi);
    }
    free (GEN->iv);
  }

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_pinv_free() */

/*****************************************************************************/

double
_unur_pinv_sample( struct unur_gen *gen )
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
  double U,X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  /* sample from U(0,1) */
  U = _unur_call_urng(gen->urng);

  /* compute inverse CDF */
  X = _unur_pinv_eval_approxinvcdf(gen,U);

  if (X<DISTR.trunc[0]) return DISTR.trunc[0];
  if (X>DISTR.trunc[1]) return DISTR.trunc[1];

  return X;

} /* end of _unur_pinv_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
     /* (internal call)                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF (0<=u<=1, no validation!)         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  int i;
  double x,un;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  /* rescale for range (0, Umax) */
  un = u * GEN->Umax;

  /* look up in guide table and search for interval */
  i = GEN->guide[(int)(u * GEN->guide_size)];
  while (GEN->iv[i+1].cdfi < un)
    i++;

  /* rescale for range (0, CDF(right)-CDF(left) for interval */
  un -= GEN->iv[i].cdfi;

  /* evaluate polynomial */
  x = _unur_pinv_newton_eval(un, GEN->iv[i].ui, GEN->iv[i].zi, GEN->order);

  /* return point (add left boundary point to x) */
  return (GEN->iv)[i].xi + x;

} /* end of _unur_pinv_eval_approxinvcdf() */

/*---------------------------------------------------------------------------*/

double
unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
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
  if ( gen->method != UNUR_METH_PINV ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_PINV_GEN,INFINITY);

  if ( u<0. || u>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }
  
  /* validate argument */
  if (u<=0.) return DISTR.domain[0];
  if (u>=1.) return DISTR.domain[1];

  /* compute inverse CDF */
  x = _unur_pinv_eval_approxinvcdf(gen,u);

  /* validate range */
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];

  return x;

} /* end of unur_pinv_eval_approxinvcdf() */

/*****************************************************************************/

int
unur_pinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
     /*----------------------------------------------------------------------*/
     /* Estimate maximal u-error and mean absolute error (MAE) by means of   */
     /* Monte-Carlo simulation.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*   max_error  ... pointer to double for storing maximal u-error       */
     /*   MAE        ... pointer to double for storing MA u-error            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double score;

  /* check arguments */
  _unur_check_NULL(GENTYPE, gen, UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* run test */
  score = unur_test_inverror(gen, max_error, MAE, 1.e-20, samplesize, 
			     FALSE, FALSE, FALSE, NULL);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_pinv_estimate_error() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/** Quadrature                                                              **/

double
_unur_pinv_lobatto5 (struct unur_gen *gen, double x, double h, double tol)
     /*----------------------------------------------------------------------*/
     /* Numerical integration of the PDF over the interval (x,x+h)           */
     /* using adaptive Gauss-Lobatto integration with 5 points.              */
     /*                                                                      */
     /* Halfs intervals recursively if error is too large.                   */
     /*                                                                      */
     /* The recursion stops when the relative error OR the absolute error    */
     /* is less than the respective given tolerances.                        */
     /* (use a negative value if only one of these tolerances has to be      */
     /* checked.)                                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left boundary point of interval                            */
     /*   h   ... length of interval                                         */
     /*   tol ... tolerated error (relative to area below PDF)               */
     /*   reclevel     ... recursion level                                   */   
     /*                                                                      */
     /* return:                                                              */
     /*   integral                                                           */
     /*----------------------------------------------------------------------*/
{ 
  double fl, fr, fc;  /* values of PDF at x-h, x, and x+h */
  double int1, int2;  /* estimated values of integrals */
  double is;          /* estimate for PDFarea rescaled by requested tolerance */

  /* points for Lobatto integration */
#define W1 (0.17267316464601146)   /* = 0.5-sqrt(3/28) */
#define W2 (1.-W1)

  /* check length of interval */
  if (_unur_iszero(h)) return 0.;

  /* compute PDF values */
  fl = PDF(x);
  fr = PDF(x+h);
  fc = PDF(x+h/2.);
 
  /* compute integral on [x,x+h] */
  int1 = (9*(fl+fr)+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64*fc)*h/180.;

  /* compute integral on [x,x+h/2] and on [x+h/2,x+h] */
  int2 = ( (9*(fl+fc)+49.*(PDF(x+h*W1*0.5)+PDF(x+h*W2*0.5))+64*PDF(x+h*0.25))*h/360. +
	   (9*(fc+fr)+49.*(PDF(x+h*(0.5+W1*0.5))+PDF(x+h*(0.5+W2*0.5)))+64*PDF(x+h*0.75))*h/360. );

  /* check whether accuracy goal is reached */
  is = GEN->area * tol / DBL_EPSILON;
  if (is + (int1-int2) == is)
    return int2;
  
  /* else: error above tolerance */

  if (x+h/2. == x) {
    _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,
		  "numeric integration did not reach full accuracy");
    return int2;
    /* Remark: Since we are halving intervals, this comparision */
    /* limits the maximal number of iterations to at most 2048. */
  }

  /* recompute with shorter intervals */
  return ( _unur_pinv_lobatto5(gen,x,    h/2,tol) +
	   _unur_pinv_lobatto5(gen,x+h/2,h/2,tol) );

#undef W1
#undef W2
} /* end of _unur_pinv_lobatto5() */


/*****************************************************************************/
/** Preprocessing                                                           **/

int
_unur_pinv_relevant_support ( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* 1a. Estimate computationally relevant domain (support) of PDF        */
     /*     (finite interval where PDF is above some threshold value).       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  /** FIXME: PDF(boundary of domain) > 0. && < INFINITY --> set sleft=FALSE **/
  /* y = PDF(domain[0]);
     if (_unur_isfinite(y) && y > 1.e-16 ) ...
  */

  /* search for interval of computational relevance (if required) */
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_searchborder(gen,DISTR.center, GEN->bleft, &(GEN->dleft));
    if (!_unur_isfinite(GEN->bleft)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get left boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  
  if(GEN->sright) {
    GEN->bright = _unur_pinv_searchborder(gen,DISTR.center, GEN->bright, &(GEN->dright));
    if (!_unur_isfinite(GEN->bright)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get right boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_relevant_support() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound, double *dom)
     /*----------------------------------------------------------------------*/
     /* [1a.] Find left or right hand border of relevant domain.             */
     /*                                                                      */
     /* Calculate domain of computational relevant region.                   */
     /* Start at 'x0' and search towards 'bound'.                            */
     /* The boundary points of this domain are approximately given as        */
     /*      PDF(x0) * PINV_PDFLLIM                                          */
     /*                                                                      */
     /* As a side effect the support of the distribution is shrinked if      */
     /* points with PDF(x)=0 are found.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   x0    ... starting point for searching boudary                     */
     /*             PDF(x0) must not be too small                            */
     /*   bound ... stop searching at this point                             */
     /*   dom   ... pointer to boundary of domain / support of distribution  */
     /*                                                                      */
     /* return:                                                              */
     /*   boundary point                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double x, xold;         /* current and old previous searching point */
  double xnull;           /* last point where PDF vanishes */
  double fx;              /* PDF at x */
  double fllim;           /* threshold value */
  int i;                  /* aux variable */

  /* threshold value where we stop searching */
  fllim = PDF(x0) * PINV_PDFLLIM;

  /* starting point */
  xold = x0;
  x = _unur_arcmean(x0,bound);
  xnull = INFINITY;

  /* find a point where PDF is less than the threshold value */
  for (i=0; i<100 && PDF(x) > fllim; i++) {
    xold = x;
    x = _unur_arcmean(x,bound);
  }

  /* however: PDF(x) must not be too small */
  do{
    x = (x+xold)*0.5;
    fx = PDF(x);
    if (_unur_iszero(fx)) xnull = x;
    i++;
  } while(i<2048 && fx<fllim);
  /* Remark:
   * The emergency break after 2048 seems very high. On the other hand we can 
   * protect ourselves agains users that provide a PDF like that of the 
   * expontential distribution without providing a domain.
   */

  /* Check whether we have to shrink the support of the distribution */
  if (_unur_isfinite(xnull))
    *dom = (fabs(xnull)<1.e-300) ? 0. : xnull;

  /* return point */
   return x;

} /* end of _unur_pinv_searchborder() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_approx_pdfarea (struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* 1b. Compute area below PDF over relevant domain approximately.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* Remark:                                                              */
     /*   The user can also provide the area below the PDF.                  */
     /*   However, then we probably need not method PINV                     */
     /*----------------------------------------------------------------------*/
{
  double tol;   /* tolerated integration error */
  int i;        /* number of trials            */
  int res = UNUR_SUCCESS; /* return code of computation */

  /* there is no need to be more accurate than the U-resolution.      */
  /* here we assume that the integral is not too much smaller than 1. */
  tol = GEN->u_resolution;


  for (i=1; i<=2; i++) {

    GEN->area  = _unur_pinv_lobatto5( gen, DISTR.center, GEN->bright - DISTR.center, tol);
    GEN->area += _unur_pinv_lobatto5( gen, GEN->bleft,   DISTR.center - GEN->bleft,  tol);

    if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
      res = UNUR_FAILURE;
      break;
    }

    
    /** FIXME: this is only necessary if we use absolute tolerance **/

    if (GEN->area < 1.e-3) {
      /* the area is too small. thus we the relative integration error is too large */
      tol *= GEN->area;
    }
    else {
      break;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_searchbd(gen,FALSE);
#endif

  return res;

} /* end of _unur_pinv_approx_pdfarea() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_computational_domain (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* 1c. Compute computational domain where inverse CDF is approximated   */
     /*     (interval where we safely can compute coefficients of            */
     /*     interpolating polynomial).                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double tailcut_error;    /* threshold values for cut-off points */

  /* parameters for tail cut-off points: maximal area in tails          */
  /* We use the given U-reslution * PINV_TAILCUTOFF_FACTOR * PDFarea.   */
  tailcut_error = GEN->u_resolution * PINV_TAILCUTOFF_FACTOR(GEN->u_resolution);
  tailcut_error = _unur_min( tailcut_error, PINV_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PINV_UERROR_CORRECTION;

  /* compute cut-off points for tails */
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_cut( gen, GEN->dleft, GEN->bleft, (GEN->bleft-GEN->bright)/128, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

  if(GEN->sright) {
    GEN->bright = _unur_pinv_cut( gen, GEN->dright, GEN->bright, (GEN->bright-GEN->bleft)/128, tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_searchbd(gen,TRUE);
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_computational_domain() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_cut( struct unur_gen *gen, double dom, double w, double dw, double crit )
     /*----------------------------------------------------------------------*/
     /* [1c.] Calculate cut-off points for computational domain of           */
     /* distribution.                                                        */
     /* The area outside the cut-off point is given by 'crit'.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   dom  ... boundary of domain / support of distribution              */
     /*   w    ... starting point for searching cut-off point                */
     /*   dw   ... initial step size for searching,                          */
     /*            sign of dw gives searching direction:                     */
     /*               dw < 0 ... left hand side cut-off point                */
     /*               dw > 0 ... right hand side cut-off point               */
     /*   crit ... u-error criterium for tail cut off                        */
     /*                                                                      */
     /* return:                                                              */
     /*   cut-off point                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double u, uplus;     /* tail probabilities */

  int found = FALSE;   /* indicates whether we have found an cut-off point */
  double dx = dw/64.;  /* step size for numeric differentiation */
  int j;               /* aux variable */

  double s = (dw>0) ? 1 : -1; /* searching direction */

  /* search for cut-off point with tail probability less than 'crit'ical value */
  for (j=1; j<1000; j++) {
    
    /* check for boundary */
    if (s*dom < s*w) {
      /* boundary exceeded */
       return dom;
    }

    /* compute approximate tail probability at w */
    u = _unur_pinv_tailprob(gen, w, dw/64.);

    /* below threshold value ? */
    if (u < crit && u >= 0.) {
      found = TRUE;
      break;
    }

    /* else
       the tail probability is too large, or the approximation formula
       is not appropriate for the point.
       Hence we make a step towards + or - infinity ...
    */
    w += dw;

    /* ... and increase stepsize for next interation. */
    if (j>32) dw *= 1.5;
  }

  if(!found){
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cound not find valid cut point");
    return INFINITY;
  }


  if (_unur_iszero(u)) {
    /* the tail probability at w is approximately 0, i.e. _unur_pinv_tailprob() */
    /* assumes it is 0. Thus we only can return this point w.                   */
    return w;
  }

  /* now run secant method to find cut-off point with tail probability approx 'crit' */

  /* step size for numeric differentiation */
  dx=dw/64.;

  for (j=0; j<2048; j++) {

    /* check whether 'u' approx 'crit' */
    if (fabs(crit/u - 1.)<1.e-7) 
      return w;

    /* compute tail probability in the next point towards 'dx' */
    /* (required for numerical derivative) */
    uplus = _unur_pinv_tailprob(gen,w+dx,dx);
    if(uplus<0){
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"negative tail probability");
      return INFINITY;
    }
    if (_unur_iszero(uplus)) 
      return w;

    /* next iteration for secant method */
    w -= dx * (1./u-1./crit)/(1./uplus-1./u);

    /* check for boundary */
    if (s*dom < s*w)
      /* boundary exceeded */
       return dom;

    /* compute tail probability */
    u = _unur_pinv_tailprob(gen,w,dx);
    if(u<0){
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"negative tail probability");
      return INFINITY;
    }
    if (_unur_iszero(u))
      return w;
  }

  /* could not find point till now */
  _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cound not find cut point");
  return INFINITY;

} /* end of _unur_pinv_cut() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_tailprob( struct unur_gen *gen, double x, double dx )
     /*----------------------------------------------------------------------*/
     /* [1c.] calculate approximate tail probability.                        */
     /* use formula                                                          */
     /*   area ~ f(x)^2 / ( (lc_f(x)+1) * abs(f'(x)) )                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... cut-off point                                              */
     /*   dx  ... step size for numeric differentiation                      */
     /*                                                                      */
     /* return:                                                              */
     /*   approximate tail probability                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return -INFINITY                                                   */
     /*----------------------------------------------------------------------*/
{
  double fx, fp, fm; /* value of PDF at x, x+dx, and x-dx */
  double df;         /* derivative of PDF                 */
  double lcplus1;    /* local concavity + 1               */
  double area;       /* tail probability                  */

  /* compute PDF */
  fx = PDF(x);
  fp = PDF(x+dx);
  fm = PDF(x-dx);

  /* We have a serious problem when PDF(X) == 0.                      */
  /* Thus our only chance is to assume that the support of the PDF is */
  /* connected. Hence we assume that we have found the cut point.     */
  /*   if ( _unur_iszero(fp) || _unur_iszero(fm) ) { */
  if ( _unur_iszero(fx) ) {
    /* TODO: we should start simple bisection to find the extremal point */
    /* where PDF(X) is non-zero.                                         */
    return 0.;
  } 

  /* check data */
  if( fm-2.*fx+fp < 0.|| fm<1.e-100 || fx<1.e-100 || fp<1.e-100) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,
		  "numerical problems with tail probability");
  }

  /* approximate local concavity + 1 */
  lcplus1 = fp/(fp-fx) + fm/(fm-fx);

  /* approximate derivative */
  df = (fp-fm)/(2.*dx);
  /** TODO: use dPDF if available **/
 
  /* approximate tail probability */
  area = (fx*fx) / (lcplus1 * fabs(df));

  /* check result */
  if (area < 0.) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"tail probability might be negative");
    /* Remark: 
       It is difficult to distinguish between an invalid PDF given by the user,
       a real numerical (round-off) error or simply a point too far from 
       the region of "computational irrelevant".
       Hence we just print a warning and let the calling routine handle the 
       situation.
     */
  }

  if (_unur_isnan(area)) {
    /* When the PDF is extremely flat than NaN might occur. There are two possibilities: 
     * (1) We are far away from the center. Then the tail probabilty is almost 0.
     * (2) PDF(X) is still quite large. Then we should return INFINITY.
     */
    _unur_warning(gen->genid,UNUR_ERR_NAN,"tail probability gives NaN --> assume 0.");
    return 0.;
  }

  /* return area below PDF in tail ("tail probability") */
  return area;

} /* end of _unur_pinv_tailprob() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Compute coefficients for polynomial                                     **/

int
_unur_pinv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Create table for Newton interpolation.                               */
     /* The computational domain must be already computed / given.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double uerrcrit;           /* threshold for maximum U-error */
  double maxerror;           /* maximum U-error in a particular interval */ 
  double h;                  /* step size / length of intervals */
  int i;                     /* number of interval at work */
  int iter;                  /* number of iterations */
  int cont;                  /* whether we have to continue or loop */
  int k;                     /* auxiliary variable */
  double chebyshev[MAX_ORDER+1]; /* Chebyshev points */
  double xval[MAX_ORDER+1];  /* x-values for construction points for Newton polynomial */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* compute construction points for Chebyshev points */
  _unur_pinv_chebyshev_points(GEN->order,chebyshev);

  /* threshold for tolerated U-error */
  uerrcrit = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION;

  /* initialize step size for subintervals */
  h = (GEN->bright-GEN->bleft)/128.;

  /* initialize array of interval: starting interval */
  if (_unur_pinv_interval( gen, 0, GEN->bleft, 0.) != UNUR_SUCCESS) 
    return UNUR_ERR_GEN_CONDITION;

  /* initialize counter and control variables */
  i = 0;        /* start at interval 0 ;-) */
  cont = TRUE;  /* there is at least one iteration */

  /* compute intervals and coefficients of Newton polynomials */
  for (iter=0; cont ; iter++) {

    /* check number of iterations */
    if (iter > PINV_MAX_ITER_IVS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "maximum number of interations exceeded");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* right boundary reached ? */
    if(GEN->iv[i].xi+h > GEN->bright) {
      h = GEN->bright - GEN->iv[i].xi;
      cont = FALSE;  /* we probably can stop after this iteration */
    }

    /* construction points for Newton interpolation polynomial: */
    /* uses Chebyshev points for x values.                      */
    for(k=0; k<=GEN->order; k++)
      xval[k] = GEN->iv[i].xi + h * chebyshev[k];

    /* compute Newton interpolation polynomial */
    if (_unur_pinv_newton_create(gen,&(GEN->iv[i]),xval,uerrcrit) != UNUR_SUCCESS)
      return UNUR_ERR_GEN_CONDITION;
    /** FIXME: make interval longer ?? **/

    /* estimate error of Newton interpolation */
    maxerror = _unur_pinv_newton_maxerror(gen,&(GEN->iv[i]),xval,uerrcrit);

    if (maxerror > uerrcrit) { 
      /* error too large: reduce step size */
      h *= (maxerror > 4.*uerrcrit) ? 0.81 : 0.9;
      cont = TRUE;  /* we need another iteration */
    }

    else {
      /* create next interval */
      if ( _unur_pinv_interval( gen, i+1, GEN->iv[i].xi+h, 
				GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order-1])
	   /* cdfi holds CDF value at the left border of the interval,                  */
	   /* ui[order-1] holds area below PDF in interval, i.e. CDF(right) - CDF(left) */
      	   != UNUR_SUCCESS )
	return UNUR_ERR_GEN_CONDITION;

      /* increase step size for very small errors */
      if(maxerror < 0.3*uerrcrit) h *= 1.2;
      if(maxerror < 0.1*uerrcrit) h *= 2.;
      
      /* continue with next interval */
      i++;
    }
  }

  /* update size of array */
  GEN->iv = _unur_xrealloc( GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  
  /* set range for uniform random numbers */
  /* Umin = 0, Umax depends on area below PDF, tail cut-off points and round-off errors */
  GEN->Umax = GEN->iv[GEN->n_ivs].cdfi;

  /* o.k. */
  return UNUR_SUCCESS;
}  /* end of _unur_pinv_create_table() */

/*---------------------------------------------------------------------------*/

int 
_unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx )
     /*----------------------------------------------------------------------*/
     /* make a new interval i with left boundary point x and CDF(x).         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   i    ... index (number) of interval                                */
     /*   x    ... left boundary point of new interval                       */
     /*   cdfx ... CDF at x                                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_pinv_interval *iv;

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);

  /* check for free intervalls */
  if (i >= PINV_MAX_IVS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"maximum number of intervals exceeded");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* set values */
  iv = GEN->iv+i;     /* pointer to interval */
  iv->xi = x;         /* left boundary of interval */
  iv->cdfi = cdfx;    /* CDF at left boundary */
  COOKIE_SET(iv,CK_PINV_IV);

  /* allocate space for coefficients for Newton interpolation */
  iv->ui = _unur_xmalloc( GEN->order * sizeof(double) );
  iv->zi = _unur_xmalloc( GEN->order * sizeof(double) );

  /* update size of array (number of intervals) */
  GEN->n_ivs = i;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_interval() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval, double uerrcrit)
     /*----------------------------------------------------------------------*/
     /* 2a. Compute coefficients for Newton interpolation within a           */
     /* subinterval of the domain of the distribution.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv       ... pointer to current interval                           */
     /*   h        ... length of interval                                    */
     /*   xval     ... x-values for constructing polynomial                  */
     /*   uerrcrit ... maximal accepted u-error                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
     /** uses Chebyshev points for x-points **/
     /** not implemented: use Chebyshev points for u-points:
	 setup less stable and more expensive (need explicit inversion),
	 but smaller tables as u-points need not be stored.
	 (also need fewer intervals) (expecially in the tails).
	 need not _unur_pinv_newton_testpoints()
	 (could use maxima of Chebyshev polynomial (closed form))
     **/

{
  double *ui = iv->ui;   /* u-values for Newton interpolation */
  double *zi = iv->zi;   /* coefficients of Newton interpolation */
  
  double xi, dxi;        /* boundary and length of i-th subinterval */

  double area;           /* integral of PDF over subinterval */
  int i,k;

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* compute tuples (ui,zi) for constructing polynomials */
  for(i=0; i<GEN->order; i++) {

    /* left boundary and length of subinterval for integration */
    xi = xval[i];
    dxi = xval[i+1]-xval[i];

    /* compute integral of PDF in interval (xi,xi+dxi) */
    area = _unur_pinv_lobatto5(gen, xi, dxi, uerrcrit*0.1);
    if (area<1.e-50) {
      /** FIXME: see below **/
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"interval too short or PDF 0");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* construction points of interpolation polynomial of CDF^{-1} */
    ui[i] = (i>0) ? (ui[i-1]+area) : area;
    /* rescaled corresponding values of CDF^{-1} */ 
    zi[i] = dxi/area;

    /** FIXME (see below): if (!_unur_isfinite(zi[i]) ... "PDF over interval too close to 0" **/

  }
  /* Remark: ui[GEN->order-1] is the probability of the interval */

  /* compute coefficients of interpolation polynomial */
  for(k=1; k<GEN->order; k++) {
    for(i=GEN->order-1; i>k; i--) {
      zi[i] = (zi[i]-zi[i-1]) / (ui[i]-ui[i-(k+1)]);
    }
    zi[k] = (zi[i]-zi[i-1]) / ui[i];
  }
  
  /** FIXME: if (!_unur_isfinite(zi[i]) ... "PDF over interval too close to 0" **/


  /* ?WH? muesste man da nicht ueberpruefen ob ui[i]-ui[i-k] != 0,
     bzw. of _unur_isfinite(zi[i]) TRUE ist,  ja sollte man machen ? Siehe TODO oben
  */

  /** FIXME: test for inf and NaN **/

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_chebyshev_points (int order, double *pt)
     /*----------------------------------------------------------------------*/
     /* [2a.] Compute Chebyshev points.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   order ... order of polynomial                                      */
     /*   pt    ... pointer to array of size (order+1) for storing points    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;
  double phi = M_PI*0.5/(order+1); /* constant for computing Chebyshev points */
  
  pt[0] = 0.;
  for(i=1; i<order; i++)
    pt[i] = sin(i*phi) * sin((i+1)*phi)/cos(phi);
  pt[order] = 1.;

  return UNUR_SUCCESS;
} /* end of _unur_pinv_chebyshev_points() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_eval ( double q, double ui[], double zi[], int order )
     /*----------------------------------------------------------------------*/
     /* 2b. evaluate Newton interpolation polynomial using Horner scheme.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   q     ... argument                                                 */
     /*   ui    ... coefficients of polynomial (increasing order)            */
     /*   zi    ... coefficients of polynomial (increasing order)            */
     /*   order ... order of polynomial                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   value of interpolating polynomial at u=q                           */
     /*----------------------------------------------------------------------*/
{
  int k;
  double chi;

  chi = zi[order-1];
  for (k=order-2; k>=0; k--)
    chi = chi*(q-ui[k])+zi[k];

  return (chi*q);
} /* end of _unur_pinv_newton_eval() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv,
			    double xval[], double uerrcrit)
     /*----------------------------------------------------------------------*/
     /* 2c. Estimate maximal error of Newton interpolation in subinterval.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   xval ... x-values for constructing polynomial                      */
     /*                                                                      */
     /* return:                                                              */
     /*   estimated maximal u-error                                          */
     /*----------------------------------------------------------------------*/
{
  double x0 = iv->xi;    /* left boundary point of interval */
  double *ui = iv->ui;   /* u-values for Newton interpolation  */
  double *zi = iv->zi;   /* coefficient of Newton interpolation  */

  double maxerror = 0.;  /* maximum error */
  double uerror;         /* error for given U value */
  double x;              /* x = CDF^{-1}(U) */
  double u;              /* u = CDF(x) */

  double testu[MAX_ORDER]; /* array of U values for testing */ 
  int i;                 /* aux variable */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* get U values for test (points with maximal worst case error) */
  _unur_pinv_newton_testpoints(GEN->order,ui,testu);

  /* calculate the max u-error at the test points */
  for(i=0; i<GEN->order; i++){
    
    /* inverse CDF for U test point */
    x = _unur_pinv_newton_eval(testu[i], ui, zi, GEN->order);
    
    /* estimate CDF for interpolated x value */
    if (i==0 || xval==NULL)
      u = _unur_pinv_lobatto5(gen, x0, x, uerrcrit*0.1);
    else
      u = ui[i-1] + _unur_pinv_lobatto5(gen, xval[i], x+x0-xval[i], uerrcrit*0.1);

    /* compute u-error */
    uerror = fabs(u - testu[i]);
    
    /* update maximal error */
    if (uerror>maxerror) maxerror = uerror;
  }

  return maxerror;
} /* end of _unur_pinv_newton_maxerror() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_testpoints (int order, double ui[], double utest[])
     /*----------------------------------------------------------------------*/
     /* [2c.] calcuates the local maxima of the polynomial.                  */
     /* used as control points for error estimate.                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*    order ... order of interpolation polynomial                       */
     /*    ui    ... u-values of interpolation                               */ 
     /*    utest ... pointer to array for storing control points             */
     /*                                                                      */
     /* return:                                                              */
     /*    u-values of control points in the array utest                     */
     /*----------------------------------------------------------------------*/
{
  int k,j,i;
  double sum, qsum,x;
  
  /* compute approximate maxima of error polynomial */
  for(k=0; k<order; k++) {

    /* first approximation: use mean of consecuting construction points */
    x = (k>0) ? 0.5*(ui[k-1]+ui[k]) : 0.5*ui[k];

    /* make two iterations for root finding */
    for(j=1; j<=2; j++) {
      sum = 1./x;
      qsum = sum*sum;
      for(i=0; i<order; i++){
	sum += 1./(x-ui[i]);
	qsum += 1./((x-ui[i])*(x-ui[i]));
      }
      x += sum/qsum;
    }

    /* store result in utest */
    utest[k] = x;
  }
  
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_testpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_make_guide_table (struct unur_gen *gen)
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
  int i,j, imax;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  GEN->guide_size = GEN->n_ivs * PINV_GUIDE_FACTOR;
  if (GEN->guide_size <= 0) GEN->guide_size = 1;
  GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) );

  /* maximum index for array of data */
  imax = GEN->n_ivs;

  /* create guide table */
  i = 0;
  GEN->guide[0] = 0;
  for( j=1; j<GEN->guide_size ;j++ ) {
    while(GEN->iv[i+1].cdfi/GEN->Umax < j/(double)GEN->guide_size && i < imax)
      /* "/GEN->Umax" above is necessary, as we need the guide table for u in (0,umax) */
      /* ?WH? ist das wirklich notwendig ? WH: Hab laenger herumgebastelt und nichts besseres gefunden??? */
      i++;
    if (i >= imax) break;
    GEN->guide[j]=i;
  }

  /* check i */
  i = _unur_min(i,imax);

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = i;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_make_guide_table() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile before setup starts.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = PINV (Polynomial interpolation based INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_pinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: order of polynomial = %d",gen->genid,GEN->order);
  _unur_print_if_default(gen,PINV_SET_ORDER);
  fprintf(log,"\n%s: u-resolution = %g",gen->genid,GEN->u_resolution);
  _unur_print_if_default(gen,PINV_SET_U_RESOLUTION);
  fprintf(log,"\n");

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} /* end of _unur_pinv_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init( const struct unur_gen *gen, int ok )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   ok  ... exitcode of init call                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: domain of computation = [%g,%g]\n",gen->genid, GEN->bleft,GEN->bright);
  fprintf(log,"%s: Umin = 0 [fixed], Umax = %g,  1-Umax = %g\n",gen->genid,
	  GEN->Umax,1.-GEN->Umax);
  fprintf(log,"%s: approx. area below PDF = %g\n",gen->genid, GEN->area);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: # Intervals = %d\n",gen->genid,GEN->n_ivs);
  fprintf(log,"%s:\n",gen->genid);

  _unur_pinv_debug_intervals(gen);

  fprintf(log,"%s: initialization %s\n",gen->genid,((ok)?"successful":"*** FAILED ***")); 
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_pinv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_searchbd( const struct unur_gen *gen, int aftercut )
     /*----------------------------------------------------------------------*/
     /* print computational before or after searching for cut-off points     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   aftercut ... whether we print boundary points after searching      */
     /*                for proper points (before computing cut-off points)   */
     /*                or after computing cut-off points                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (aftercut) {
    fprintf(log,"%s: after cutting-off points: domain = (%g,%g)\n",gen->genid,
	    GEN->bleft,GEN->bright);
  }
  else {
    fprintf(log,"%s: after searching border:   domain = (%g,%g),  area = %g\n",gen->genid,
	    GEN->bleft,GEN->bright,GEN->area);
    fprintf(log,"%s: possible support of distribution = (%g,%g)\n",gen->genid,
	    GEN->dleft,GEN->dright);
  }

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} /* end of _unur_pinv_debug_searchbd() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print table of intervals into logfile                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int n;
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (gen->debug & PINV_DEBUG_TABLE) {
    for (n=0; n<=GEN->n_ivs; n++) {
      fprintf(log,"%s: [%3d] xi = %g, cdfi = %g\n",gen->genid,
	      n, GEN->iv[n].xi, GEN->iv[n].cdfi);
    }
  }

  fflush(log);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_pinv_debug_intervals() */


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_pinv_info( struct unur_gen *gen, int help )
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
  struct unur_distr *distr = gen->distr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.trunc[0],DISTR.trunc[1]);
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }

  if (help) {
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
                          "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
  
  /* method */
  _unur_string_append(info,"method: PINV (Polynomial interpolation based INVerse CDF)\n");
  _unur_string_append(info,"   order of polynomial = %d\n", GEN->order);
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   truncated domain = (%g,%g)\n",GEN->bleft,GEN->bright);
  /*   _unur_string_append(info,"   Prob(X<domain)   = %g\n", _unur_max(0,GEN->tailcutoff_left)); */
  /*   _unur_string_append(info,"   Prob(X>domain)   = %g\n", _unur_max(0,1.-GEN->tailcutoff_right)); */
  {
    double max_error=1.; double MAE=1.;
    unur_pinv_estimate_error( gen, 10000, &max_error, &MAE );
    _unur_string_append(info,"   u-error         <= %g  (mean = %g)\n", max_error, MAE);
  }

  _unur_string_append(info,"   # intervals      = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   order = %d  %s\n", GEN->order,
 			(gen->set & PINV_SET_ORDER) ? "" : "[default]");

    _unur_string_append(info,"   u_resolution = %g  %s\n", GEN->u_resolution,
 			(gen->set & PINV_SET_U_RESOLUTION) ? "" : "[default]");
    
    _unur_string_append(info,"   boundary = (%g,%g)  %s\n", GEN->bleft_par, GEN->bright_par,
			(gen->set & PINV_SET_BOUNDARY) ? "" : "[default]");

    _unur_string_append(info,"   search for boundary: left=%s,  right=%s  %s\n",
			(GEN->sleft ? "TRUE":"FALSE"), (GEN->sright ? "TRUE":"FALSE"), 
			(gen->set & PINV_SET_BOUNDARY) ? "" : "[default]");

    _unur_string_append(info,"\n");
  }


  /* Hints */
  if (help) {
    if ( GEN->order < MAX_ORDER )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can increase \"order\" to decrease #intervals");
    if (! (gen->set & PINV_SET_U_RESOLUTION) )
      _unur_string_append(info,"[ Hint: %s\n\t%s ]\n",
			  "You can decrease the u-error by decreasing \"u_resolution\".",
			  "(it is bounded by the machine epsilon, however.)");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_tdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
