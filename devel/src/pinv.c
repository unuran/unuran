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
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the CDF                                                   *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      pointer to PDF and dPDF                                              *
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
 *   [1] W. Hörmann and J. Leydold:                                          *
 *       Continuous Random Variate Generation by Fast Numerical Inversion,   *
 *       ACM Trans. Model. Comput. Simul. 13(4), pp. 347-362 (2003)          *
 *                                                                           *
 *   [2] W. Hörmann, J. Leydold, and G. Derflinger:                          *
 *       Automatic Nonuniform Random Variate Generation,                     *
 *       Springer-Verlag, Berlin Heidelberg (2004)                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "pinv.h"
#include "pinv_struct.h"


/*---------------------------------------------------------------------------*/
/* Constants                                                                 */



/*   PAR->guide_factor = factor; */
#define PINV_GUIDE_FACTOR  (1)

/*   PAR->max_ivs = max_ivs; */
#define PINV_MAX_IVS  (1000000)



#define CK_PINV_IV     0x00100132u


#define PINVMAXINT  10000

/* #define PINV_MAX_ITER      (300) */
/* Maximal number of iterations for finding the boundary of the              */
/* computational interval, i.e. where CDF(x) is close to 0 and 1, resp.      */

/* #define PINV_MAX_U_LENGTH  (0.05) */
/* Maximal value for |u_i - u_{i-1}|. If for an interval this value is       */
/* larger then it is splitted (independently of its u-error).                */

/* #define PINV_TAILCUTOFF    (1.e-10)  */
/* For unbounded domains the tails has to be cut of. We use the given        */
/* u-resolution for finding the cut points. (The probability of the chop     */
/* regions should be less than the 1 fifth of the u-resolution.)             */
/* However, it should not be greater than some threshold value, given by     */
/* PINV_TAILCUTOFF which reflects the precision of the used stream of        */
/* uniform pseudo-random numbers (typically about 2^32).                     */
/* However, for computational reasons we use a value that is at least twice  */
/* the machine epsilon for the right hand boundary.                          */

/* #define PINV_XDEVIATION    (0.05) */
/* Used for splitting intervals. When the u-error is estimated for an        */
/* interval then the CDF is evaluated in the approximate center of the       */
/* u-interval. This could be used as splitting point of the interval.        */
/* However, this might result in slow convergence. A much more stable        */
/* point is the center of the x-interval. However, this requires an          */
/* additional evalution of the CDF.                                          */
/* Thus we use the following rule: If CDF(approx. center of u-int) is        */
/* close to the center of the x-interval use the first, otherwise use the    */
/* latter. PINV_XDEVIATION is the threshold value for relative distance      */
/* between these two points.                                                 */
/* As a rule-of-thumb larger values of PINV_XDEVIATION result in more        */
/* intervals but less evaluations of the CDF (until there are too many       */
/* intervals).                                                               */

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
#define PINV_DEBUG_CHG       0x00001000u   /* print changed parameters       */
#define PINV_DEBUG_SAMPLE    0x01000000u   /* trace sampling                 */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define PINV_SET_ORDER          0x001u  /* order of polynomial               */
#define PINV_SET_U_RESOLUTION   0x002u  /* maximal error in u                */
#define PINV_SET_STP            0x004u  /* starting design points            */
#define PINV_SET_BOUNDARY       0x008u  /* boundary of computational region  */
#define PINV_SET_SEARCHBOUNDARY 0x010u  /* search for boundary               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "PINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_reinit( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_eval_approxinvcdf( const struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* evaluate Hermite interpolation of inverse CDF at u.                       */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_find_boundary( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* find boundary of computational interval                                   */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_create_table( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* create the table with splines                                             */
/*---------------------------------------------------------------------------*/

/* static struct unur_pinv_interval *_unur_pinv_interval_new( struct unur_gen *gen, double p, double u ); */
/*---------------------------------------------------------------------------*/
/* make a new interval with node (u=F(p),p).                                 */
/*---------------------------------------------------------------------------*/

/* static struct unur_pinv_interval *_unur_pinv_interval_adapt( struct unur_gen *gen,  */
/* 							     struct unur_pinv_interval *iv,  */
/* 							     int *error_count_shortinterval ); */
/*---------------------------------------------------------------------------*/
/* check parameters in interval and split or truncate where necessary.       */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_interval_is_monotone( struct unur_gen *gen, struct unur_pinv_interval *iv ); */
/*---------------------------------------------------------------------------*/
/* check whether the given interval is monotone.                             */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_interval_parameter( struct unur_gen *gen, struct unur_pinv_interval *iv ); */
/*---------------------------------------------------------------------------*/
/* compute all parameter for interval (spline coefficients).                 */
/*---------------------------------------------------------------------------*/

/* static double _unur_pinv_eval_polynomial( double x, double *coeff, int order ); */
/*---------------------------------------------------------------------------*/
/* evaluate polynomial.                                                      */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_list_to_array( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* copy list of intervals into double array.                                 */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_make_guide_table( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

/* static double _unur_pinv_CDF( const struct unur_gen *gen, double x ); */
/*---------------------------------------------------------------------------*/
/* compute CDF of truncated distribution.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_init( const struct unur_gen *gen, int ok);
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

/* static void _unur_pinv_debug_intervals( const struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into logfile.               */
/*---------------------------------------------------------------------------*/

/* static void _unur_pinv_debug_chg_truncated( const struct unur_gen *gen); */
/*---------------------------------------------------------------------------*/
/* trace changes of the truncated domain.                                    */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
/* static void _unur_pinv_info( struct unur_gen *gen, int help ); */
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif







int  pinvsetup ( struct unur_gen *gen );
                             
int setup (struct unur_gen *gen, double a, double b, double hh, double uerror);

static double lobato5 (struct unur_gen *gen, double x, double h, double fx, double *fxph);

static double nint_12 (struct unur_gen *gen, double a,double b,double *res_relerror);

static double nint_monoton_dens (struct unur_gen *gen, double a,double b,double step,double crit);

static double evalnewtoninterpol (double u,int g,double ui[],double zi[]);

static int newtoninterpol (struct unur_gen *gen, double x0, double h,double ui[],double zi[],double *x);

static int tstpt (int g,double ui[],double utest[]);

static double maxerrornewton (struct unur_gen *gen,double ui[],double zi[],double x0,double xval[]);

static double searchborder (struct unur_gen *gen, double x0, double step,double border);

static double tail (struct unur_gen *gen, double x, double d);

static double cut (struct unur_gen *gen, double w,double dw, double crit);














/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_pinv_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_pinv_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

/* CDF, PDF, and dPDF are rescaled such that the CDF is a "real" CDF with    */
/* u (range) in (0,1) on the interval (DISTR.domain[0], DISTR.domain[1]).    */
/* call to CDF: */
/* #define CDF(x)  (_unur_pinv_CDF((gen),(x))) */
/* --> ((_unur_cont_CDF((x),(gen->distr))-GEN->CDFmin)/(GEN->CDFmax-GEN->CDFmin)) */

/* call to PDF: */
/* #define PDF(x)  (_unur_cont_PDF((x),(gen->distr))/(GEN->CDFmax-GEN->CDFmin))  */
#define PDF(x)  (_unur_cont_PDF((x),(gen->distr)))    /* call to PDF         */

/* call to derivative of PDF: */   
/* #define dPDF(x) (_unur_cont_dPDF((x),(gen->distr))/(GEN->CDFmax-GEN->CDFmin)) */

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
/*   PAR->stp = NULL;                /\* starting nodes                          *\/ */
/*   PAR->n_stp = 0;                 /\* number of starting nodes                *\/ */

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
  if (order<2 || order>19) {
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
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution");
    return UNUR_ERR_PAR_SET;
  }
  if (u_resolution < 5.*DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution");
    u_resolution = 5.*DBL_EPSILON;
  }

  /* store date */
  PAR->u_resolution = u_resolution;

  /* changelog */
  par->set |= PINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

/* int */
/* unur_pinv_set_cpoints( struct unur_par *par, const double *stp, int n_stp ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* set starting construction points (nodes) for Hermite interpolation.  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   par    ... pointer to parameter for building generator object      *\/ */
/*      /\*   stp    ... pointer to array of starting points                     *\/ */
/*      /\*   n_stp  ... number of starting points                               *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i; */

/*   /\* check arguments *\/ */
/*   _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL ); */
/*   _unur_check_par_object( par, PINV ); */

/*   /\* check starting construction points *\/ */
/*   if (n_stp < 1 || stp==NULL) { */
/*     _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 1"); */
/*     return UNUR_ERR_PAR_SET; */
/*   } */

/*   /\* starting points must be strictly monontonically increasing *\/ */
/*   for( i=1; i<n_stp; i++ ) */
/*     if (stp[i] <= stp[i-1]) { */
/*       _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing"); */
/*       return UNUR_ERR_PAR_SET; */
/*     } */

/*   /\* store date *\/ */
/*   PAR->stp = stp; */
/*   PAR->n_stp = n_stp; */

/*   /\* changelog *\/ */
/*   par->set |= PINV_SET_STP; */

/*   return UNUR_SUCCESS; */

/* } /\* end of unur_pinv_set_cpoints() *\/ */

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
  if (left <= -INFINITY || right >= INFINITY) {
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

/* int */
/* unur_pinv_get_n_intervals( const struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* get number of intervals (or more precisely the number of nodes)      *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen  ... pointer to generator object                               *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   number of intervals ... on success                                 *\/ */
/*      /\*   0     ... on error                                                 *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   /\* check input *\/ */
/*   _unur_check_NULL( GENTYPE, gen, 0 ); */
/*   _unur_check_gen_object( gen, PINV, 0 ); */
/*   return GEN->N; */
/* } /\* end of unur_pinv_get_n_intervals() *\/ */

/*---------------------------------------------------------------------------*/

/* int  */
/* unur_pinv_chg_truncated( struct unur_gen *gen, double left, double right ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* change the left and right borders of the domain of the distribution  *\/ */
/*      /\* the new domain should not exceed the original domain given by        *\/ */
/*      /\* unur_distr_cont_set_domain(). Otherwise it is truncated.             *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen      ... pointer to generator object                           *\/ */
/*      /\*   left  ... left boundary point                                      *\/ */
/*      /\*   right ... right boundary point                                     *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*                                                                      *\/ */
/*      /\* comment:                                                             *\/ */
/*      /\*   the new boundary points may be +/- INFINITY                        *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   double Umin, Umax, Uminbound, Umaxbound; */

/*   /\* check arguments *\/ */
/*   _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL ); */
/*   _unur_check_gen_object(gen, PINV, UNUR_ERR_GEN_INVALID); */

/*   /\* the truncated domain must be a subset of (computational) domain *\/ */
/*   if (left < GEN->bleft) { */
/*     _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, increase left boundary"); */
/*     left = GEN->bleft; */
/*   } */
/*   if (right > GEN->bright) { */
/*     _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, decrease right boundary"); */
/*     right = GEN->bright; */
/*   } */

/*   /\* the truncated domain must have non-empty intersection *\/ */
/*   if (!_unur_FP_less(left,right)) { */
/*     _unur_error(gen->genid,UNUR_ERR_DISTR_SET,"domain, left >= right"); */
/*     return UNUR_ERR_DISTR_SET; */
/*   } */

/*   /\* compute Uminbound and Umaxbound using the u-value of the first and  */
/*      the last design point. *\/ */
/*   /\* this setting of Uminbound and Umaxbound guarantees that in the  */
/*      sampling algorithm U is always in a range where a table */
/*      is available for the inverse CDF. */
/*      So this is a safe guard against segfault for U=0. or U=1. *\/  */
/*   Uminbound = _unur_max(0.,GEN->intervals[0]); */
/*   Umaxbound = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]); */

/*   /\* set bounds of U -- in respect to given bounds *\/ */
/*   Umin = (left > -INFINITY) ? CDF(left)  : 0.; */
/*   Umax = (right < INFINITY) ? CDF(right) : 1.; */

/*   /\* check result *\/ */
/*   if (Umin > Umax) { */
/*     /\* this is a serios error that should not happen *\/ */
/*     _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); */
/*     return UNUR_ERR_SHOULD_NOT_HAPPEN; */
/*   } */

/*   if (_unur_FP_equal(Umin,Umax)) { */
/*     /\* CDF values very close *\/ */
/*     _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close"); */
/*     if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) { */
/*       /\* this is very bad *\/ */
/*       _unur_error(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close"); */
/*       return UNUR_ERR_DISTR_SET; */
/*     } */
/*   } */

/*   /\* copy new boundaries into generator object *\/ */
/*   DISTR.trunc[0] = left; */
/*   DISTR.trunc[1] = right; */
/*   GEN->Umin = _unur_max(Umin, Uminbound); */
/*   GEN->Umax = _unur_min(Umax, Umaxbound); */

/*   /\* changelog *\/ */
/*   gen->distr->set |= UNUR_DISTR_SET_TRUNCATED; */

/* #ifdef UNUR_ENABLE_LOGGING */
/*   /\* write info into log file *\/ */
/*   if (gen->debug & PINV_DEBUG_CHG)  */
/*     _unur_pinv_debug_chg_truncated( gen ); */
/* #endif */
  
/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */
  
/* } /\* end of unur_pinv_chg_truncated() *\/ */


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



  pinvsetup(gen);


  /* compute splines */
/*   if (_unur_pinv_create_table(gen)!=UNUR_SUCCESS) { */
/*     /\* make entry in log file *\/ */
/* #ifdef UNUR_ENABLE_LOGGING */
/*     if (gen->debug) { */
/*       _unur_pinv_list_to_array( gen ); */
/*       _unur_pinv_debug_init(gen,FALSE); */
/*     } */
/* #endif */
/*     _unur_pinv_free(gen); return NULL; */
/*   } */

  /* copy linked list into array */
/*   _unur_pinv_list_to_array( gen ); */

  /* adjust minimal and maximal U value */
/*   GEN->Umin = _unur_max(0.,GEN->intervals[0]); */
/*   GEN->Umax = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]); */

  /* this setting of Umin and Umax guarantees that in the
     sampling algorithm U is always in a range where a table
     is available for the inverse CDF.
     So this is a safe guard against segfault for U=0. or U=1. */

  /* These values for Umin and Umax are only changed in
     unur_pinv_chg_truncated(). */

  /* make guide table */
/*   _unur_pinv_make_guide_table(gen); */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_pinv_debug_init(gen,TRUE);
#endif

  /* we do not use 'stp' any more. since it contains a pointer to
     an array outside the generator object. thus we set the pointer
     to NULL
  */
/*   GEN->stp = NULL; */
/*   GEN->n_stp = 0; */

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

/*   /\* compute splines *\/ */
/*   if ( (rcode = _unur_pinv_create_table(gen)) != UNUR_SUCCESS) */
/*     return rcode; */

/*   /\* copy linked list into array *\/ */
/*   _unur_pinv_list_to_array( gen ); */

/*   /\* adjust minimal and maximal U value *\/ */
/*   GEN->Umin = _unur_max(0.,GEN->intervals[0]); */
/*   GEN->Umax = _unur_min(1.,GEN->intervals[(GEN->N-1)*(GEN->order+2)]); */

/*   /\* (re)set sampling routine *\/ */
/*   SAMPLE = _unur_pinv_getSAMPLE(gen); */

/* #ifdef UNUR_ENABLE_LOGGING */
/*   /\* write info into log file *\/ */
/*   if (gen->debug & PINV_DEBUG_REINIT) _unur_pinv_debug_init(gen,TRUE); */
/* #endif */

/*   /\* make guide table *\/ */
/*   _unur_pinv_make_guide_table(gen); */

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
/*   gen->reinit = _unur_pinv_reinit; */

  /* copy parameters into generator object */
  GEN->order = PAR->order;            /* order of polynomial                 */
  GEN->u_resolution = PAR->u_resolution; /* maximal error in u-direction     */
  GEN->bleft_par  = PAR->bleft;          /* border of computational domain   */
  GEN->bright_par = PAR->bright;
  GEN->sleft  = PAR->sleft;              /* whether to search for boundary   */
  GEN->sright = PAR->sright;

/*   GEN->stp = PAR->stp;               /\* pointer to array of starting points  *\/ */
/*   GEN->n_stp = PAR->n_stp;           /\* number of construction points        *\/ */

  /* default values */
/*   GEN->tailcutoff_left  = -1.;       /\* no cut-off by default                *\/ */
/*   GEN->tailcutoff_right = 10.; */

  /* initialize variables */
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;
  GEN->Umin = 0.;
  GEN->Umax = 1.;
  GEN->umax = 1.;
  /*   GEN->N = 0; */
  GEN->iv = NULL;
  GEN->ni = 0;
  /*   GEN->intervals = NULL; */
  GEN->guide_size = 0; 
  GEN->guide = NULL;



  GEN->iv =  _unur_xmalloc(sizeof(struct siv)*PINVMAXINT);


#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
/*   gen->info = _unur_pinv_info; */
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
  /* computational domain as given by user */
  GEN->bleft = GEN->bleft_par;
  GEN->bright = GEN->bright_par;

  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

/*   /\* set bounds of U -- in respect to given bounds                          *\/ */
/*   GEN->CDFmin = (DISTR.domain[0] > -INFINITY) ? _unur_cont_CDF((DISTR.domain[0]),(gen->distr)) : 0.; */
/*   GEN->CDFmax = (DISTR.domain[1] < INFINITY)  ? _unur_cont_CDF((DISTR.domain[1]),(gen->distr)) : 1.; */

/*   if (!_unur_FP_less(GEN->CDFmin,GEN->CDFmax)) { */
/*     _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing"); */
/*     return UNUR_ERR_GEN_DATA; */
/*   } */

/*   /\* cut points for tails *\/ */
/*   if (DISTR.domain[0] <= -INFINITY ||  */
/*       (DISTR.pdf!=NULL && _unur_cont_PDF((DISTR.domain[0]),(gen->distr))<=0.) ) { */
/*     GEN->tailcutoff_left = _unur_min(PINV_TAILCUTOFF, 0.1*GEN->u_resolution); */
/*     GEN->tailcutoff_left = _unur_max(GEN->tailcutoff_left,2*DBL_EPSILON); */
/*   } */
/*   if (DISTR.domain[1] >= INFINITY ||  */
/*       (DISTR.pdf!=NULL && _unur_cont_PDF((DISTR.domain[1]),(gen->distr))<=0.) ) { */
/*     GEN->tailcutoff_right = _unur_min(PINV_TAILCUTOFF, 0.1*GEN->u_resolution); */
/*     GEN->tailcutoff_right = _unur_max(GEN->tailcutoff_right,2*DBL_EPSILON); */
/*     GEN->tailcutoff_right = 1. - GEN->tailcutoff_right; */
/*   } */

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

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

/*   /\* copy tables for generator object *\/ */
/*   CLONE->intervals = _unur_xmalloc( GEN->N*(GEN->order+2) * sizeof(double) ); */
/*   memcpy( CLONE->intervals, GEN->intervals, GEN->N*(GEN->order+2) * sizeof(double) ); */
/*   CLONE->guide = _unur_xmalloc( GEN->guide_size * sizeof(int) ); */
/*   memcpy( CLONE->guide, GEN->guide, GEN->guide_size * sizeof(int) ); */

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

  /* free tables */
/*   if (GEN->intervals) free (GEN->intervals); */
  if (GEN->guide)     free (GEN->guide);

  if (GEN->iv) {
    for(i=0;i<=GEN->ni;i++){
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

  /* sample from U( Umin, Umax ) */
  U = GEN->Umin + _unur_call_urng(gen->urng) * (GEN->Umax - GEN->Umin);

  /* compute inverse CDF */
  X = U;
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

/*   /\* look up in guide table and search for interval *\/ */
/*   i =  GEN->guide[(int) (GEN->guide_size*u)]; */
/*   while (u > GEN->intervals[i+GEN->order+2]) */
/*     i += GEN->order+2; */

/*   /\* rescale uniform random number *\/ */
/*   u = (u-GEN->intervals[i])/(GEN->intervals[i+GEN->order+2] - GEN->intervals[i]); */

/*   /\* evaluate polynome *\/ */
/*   return _unur_pinv_eval_polynomial( u, GEN->intervals+i+1, GEN->order ); */
  

  un = u * GEN->umax;
  for(i= GEN->guide[(int)(u * GEN->guide_size)]; GEN->iv[i+1].cdfi < un ; i++)
    ;

  x = evalnewtoninterpol(un- GEN->iv[i].cdfi,GEN->order,GEN->iv[i].ui,GEN->iv[i].zi);

  return (GEN->iv)[i].xi+x;

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

/* int */
/* unur_pinv_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* Estimate maximal u-error and mean absolute error (MAE) by means of   *\/ */
/*      /\* Monte-Carlo simulation.                                              *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen        ... pointer to generator object                         *\/ */
/*      /\*   samplesize ... sample size for Monte Carlo simulation              *\/ */
/*      /\*   max_error  ... pointer to double for storing maximal u-error       *\/ */
/*      /\*   MEA        ... pointer to double for storing MA u-error            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* {  */
/*   double U, ualt, X; */
/*   double max=0., average=0., uerror, errorat=0.; */
/*   int j, outside_interval=0; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE); */

/*   for(j=0;j<samplesize;j++) {   */
/*     /\* sample from U( Umin, Umax ) *\/ */
/*     U = GEN->Umin + _unur_call_urng(gen->urng) * (GEN->Umax - GEN->Umin); */
/*     ualt=U; */

/*     /\* compute inverse CDF *\/ */
/*     X = _unur_pinv_eval_approxinvcdf(gen,U); */

/*     if (X<DISTR.trunc[0]) { X = DISTR.trunc[0]; outside_interval++; } */
/*     if (X>DISTR.trunc[1]) { X = DISTR.trunc[1]; outside_interval++; } */

/*     uerror = fabs(ualt-CDF(X)); */

/*     average += uerror; */
/*     if(uerror>max) { */
/*       max = uerror; */
/*       errorat = X; */
/*     } */
/*     /\* printf("j %d uerror %e maxerror %e average %e\n",j,uerror,max,average/(j+1)); *\/ */
/*   } */
/*   /\* */
/*   printf("maximal error occured at x= %.16e percentage outside interval %e \n", */
/* 	 errorat,(double)outside_interval/(double)samplesize);  */
/*   *\/ */

/*   *max_error = max; */
/*   *MAE = average/samplesize; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */

/* } /\* end of unur_pinv_estimate_error() *\/ */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/* int */
/* _unur_pinv_find_boundary( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* find boundary of computational interval                              *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer generator object                                   *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   double x,u; */
/*   int i; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE); */

/*   /\* reset counter for intervals *\/ */
/*   GEN->N = 0; */

/*   /\* boundary of the computational domain must not exceed domain of distribution *\/ */
/*   if (GEN->bleft  < DISTR.domain[0]) GEN->bleft  = DISTR.domain[0]; */
/*   if (GEN->bright > DISTR.domain[1]) GEN->bright = DISTR.domain[1]; */

/*   /\* find left boundary point *\/ */
/*   for (x = GEN->bleft, i=0; i<PINV_MAX_ITER; i++) { */
    
/*     /\* next trial *\/ */
/*     GEN->bleft = x; */
/*     u = CDF(GEN->bleft); */

/*     /\* everything fine ? *\/ */
/*     if (u <= GEN->tailcutoff_left || GEN->tailcutoff_left < 0.)  */
/*       break; */

/*     /\* otherwise ... *\/ */
/*     if (DISTR.domain[0] <= -INFINITY) { */
/*       /\* domain not bounded from below *\/ */
/*       x = (GEN->bleft > -1.) ? -1. : 10.*GEN->bleft; */
/*       if (! _unur_isfinite(x) )   */
/* 	i = PINV_MAX_ITER; */
/*     } */
    
/*     else { */
/*       /\* domain bounded from below *\/ */
/*       x = _unur_arcmean(GEN->bleft, DISTR.domain[0]); */
/*       if (_unur_FP_equal(x,DISTR.domain[0]))  */
/* 	i = PINV_MAX_ITER; */
/*     } */
/*   } */
  
/*   /\* computation successful ? *\/ */
/*   if (i >= PINV_MAX_ITER) */
/*     _unur_warning(gen->genid,UNUR_ERR_DISTR_PROP,"cannot find l.h.s. of domain"); */

/*   /\* make l.h.s. starting interval *\/ */
/*   GEN->iv = _unur_pinv_interval_new(gen,GEN->bleft,u); */
/*   if (GEN->iv == NULL) return UNUR_ERR_GEN_DATA; */


/*   /\* find right boundary point *\/ */
/*   for (x = GEN->bright, i=0; i<PINV_MAX_ITER; i++) { */
    
/*     /\* next trial *\/ */
/*     GEN->bright = x; */
/*     u = CDF(GEN->bright); */

/*     /\* everything fine ? *\/ */
/*     if (u >= GEN->tailcutoff_right || GEN->tailcutoff_right > 1.1)  */
/*       break; */

/*     /\* otherwise ... *\/ */
/*     if (DISTR.domain[1] >= INFINITY) { */
/*       /\* domain not bounded from above *\/ */
/*       x = (GEN->bright < 1.) ? 1. : 10.*GEN->bright; */
/*       if (! _unur_isfinite(x) )   */
/* 	i = PINV_MAX_ITER; */
/*     } */
    
/*     else { */
/*       /\* domain bounded from below *\/ */
/*       x = _unur_arcmean(GEN->bright, DISTR.domain[1]); */
/*       if (_unur_FP_equal(x,DISTR.domain[1]))  */
/* 	i = PINV_MAX_ITER; */
/*     } */
/*   } */
  
/*   /\* computation successful ? *\/ */
/*   if (i >= PINV_MAX_ITER) */
/*     _unur_warning(gen->genid,UNUR_ERR_DISTR_PROP,"cannot find r.h.s. of domain"); */

/*   /\* make r.h.s. starting intervals *\/ */
/*   GEN->iv->next = _unur_pinv_interval_new(gen,GEN->bright,u); */
/*   if (GEN->iv->next == NULL) return UNUR_ERR_GEN_DATA; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */

/* } /\* end of _unur_pinv_find_boundary() *\/ */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_create_table( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* create a table of splines                                            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer generator object                                   *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   struct unur_pinv_interval *iv, *iv_new; */
/*   int i, error_count_shortinterval=0; */
/*   double Fx; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE); */

/*   /\* find boundary point of computational interval *\/ */
/*   if (_unur_pinv_find_boundary(gen) != UNUR_SUCCESS) */
/*     return UNUR_ERR_GEN_DATA; */

/*   /\* use starting design points of given *\/ */
/*   if (GEN->stp) { */
/*     iv = GEN->iv; */
/*     for (i=0; i<GEN->n_stp; i++) { */
/*       if (!_unur_FP_greater(GEN->stp[i],GEN->bleft)) continue; /\* skip *\/ */
/*       if (!_unur_FP_less(GEN->stp[i],GEN->bright))   break;    /\* no more points *\/ */
 
/*       Fx = CDF(GEN->stp[i]); */
/*       iv_new = _unur_pinv_interval_new(gen,GEN->stp[i],Fx); */
/*       if (iv_new == NULL) return UNUR_ERR_GEN_DATA; */
/*       iv_new->next = iv->next; */
/*       iv->next = iv_new; */
/*       iv = iv_new; */

/*       if (Fx > GEN->tailcutoff_right) */
/* 	/\* there is no need to add another starting point in the r.h. tail *\/ */
/* 	break; */
/*     } */
/*   } */

/*   else /\* mode - if known - is inserted as "default design point" *\/ */
/*     if( (gen->distr->set & UNUR_DISTR_SET_MODE) && */
/*         _unur_FP_greater(DISTR.mode, GEN->bleft) && */
/*         _unur_FP_less(DISTR.mode, GEN->bright) ) { */
/*       iv = GEN->iv; */
/*       iv_new = _unur_pinv_interval_new(gen,DISTR.mode,CDF(DISTR.mode)); */
/*       if (iv_new == NULL) return UNUR_ERR_GEN_DATA; */
/*       iv_new->next = iv->next; */
/*       iv->next = iv_new; */
/*     } */

/*   /\* now split intervals where approximation error is too large *\/ */
/*   for (iv=GEN->iv; iv->next!=NULL; ) { */
/*     COOKIE_CHECK(iv,CK_PINV_IV,UNUR_ERR_COOKIE); */
/*     if (GEN->N >= PINV_MAX_IVS) { */
/*       /\* emergency break *\/ */
/*       _unur_error(GENTYPE,UNUR_ERR_GEN_CONDITION,"too many intervals"); */
/*       return UNUR_ERR_GEN_CONDITION; */
/*     } */
/*     iv = _unur_pinv_interval_adapt(gen,iv, &error_count_shortinterval); */
/*     if (iv == NULL) return UNUR_ERR_GEN_DATA; */
/*   } */

/*   /\* last interval is only used to store right boundary *\/ */
/*   iv->spline[0] = iv->p; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */
/* }  /\* end of _unur_pinv_create_table() *\/ */

/*---------------------------------------------------------------------------*/

/* struct unur_pinv_interval * */
/* _unur_pinv_interval_new( struct unur_gen *gen, double p, double u ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* make a new interval with node (u=F(p),p).                            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*   p   ... left design point of new interval                          *\/ */
/*      /\*   u   ... value of CDF at p, u=CDF(p)                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   pointer to new interval                                            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* error:                                                               *\/ */
/*      /\*   return NULL                                                        *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   struct unur_pinv_interval *iv; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,NULL); */

/*   /\* first check u *\/ */
/*   if (u<0.) { */
/*     if (u < -UNUR_SQRT_DBL_EPSILON) { */
/*       _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF(x) < 0."); */
/*       return NULL; */
/*     } */
/*     else { /\* round off error *\/ */
/*       u = 0.; */
/*     } */
/*   } */
/*   if (u>1.) { */
/*     _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF(x) > 1."); */
/*     return NULL; */
/*   } */

/*   /\* we need new interval *\/ */
/*   iv = _unur_xmalloc( sizeof(struct unur_pinv_interval) ); */
/*   COOKIE_SET(iv,CK_PINV_IV); */

/*   /\* compute and store data *\/ */
/*   switch (GEN->order) { */
/*   case 5: */
/*     iv->df = dPDF(p); */
/*   case 3: */
/*     iv->f = PDF(p); */
/*   case 1: */
/*     iv->p = p; */
/*     iv->u = u; */
/*     break; */
/*   default: */
/*     _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); */
/*     free(iv); */
/*     return NULL; */
/*   } */

/*   iv->next = NULL;  /\* add eol marker *\/ */
/*   ++(GEN->N);   /\* increment counter for intervals *\/ */

/*   /\* o.k. *\/ */
/*   return iv; */

/* } /\* end of _unur_pinv_interval_new() *\/ */

/*---------------------------------------------------------------------------*/

/* struct unur_pinv_interval * */
/* _unur_pinv_interval_adapt( struct unur_gen *gen, struct unur_pinv_interval *iv, */
/*                            int *error_count_shortinterval ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* check parameters in interval and split or truncate where necessary.  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*   iv  ... pointer to interval                                        *\/ */
/*      /\*   error_count_shortinterval ... pointer to errorcount to supress too *\/ */
/*      /\*                                 many error messages                  *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   iv       ... if splitted                                           *\/ */
/*      /\*   iv->next ... if interval was o.k.                                  *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   double p_new;   /\* new design point *\/ */
/*   struct unur_pinv_interval *iv_new, *iv_tmp; */
/*   double x, Fx; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,NULL);       COOKIE_CHECK(gen,CK_PINV_GEN,NULL); */
/*   CHECK_NULL(iv,NULL);        COOKIE_CHECK(iv,CK_PINV_IV,NULL); */
/*   CHECK_NULL(iv->next,NULL);  COOKIE_CHECK(iv->next,CK_PINV_IV,NULL); */

/*   /\* 1st check: right most interval (of at least 2) */
/*      with CDF greater than GEN->tailcutoff_right *\/ */

/*   iv_tmp = iv->next->next; */
/*   if(iv_tmp && iv->next->u > GEN->tailcutoff_right) { */
/*     /\* chop off right hand tail *\/ */
/*     free (iv_tmp); */
/*     iv->next->next = NULL; */
/*     GEN->N--; */
/*     /\* update right boundary *\/ */
/*     GEN->bright = iv->next->p; */
/*     return iv; */
/*   } */

/*   /\* 2nd check: is the left most interval (of at least 2)  */
/*      with CDF less than GEN->tailcutoff_left *\/ */

/*   if (iv==GEN->iv && iv->next->next && iv->next->u < GEN->tailcutoff_left) { */
/*     /\* chop off left hand tail *\/ */
/*     iv_tmp = GEN->iv; */
/*     GEN->iv = iv->next; */
/*     free (iv_tmp); */
/*     GEN->N--; */
/*     /\* update left boundary *\/ */
/*     GEN->bleft = GEN->iv->p; */
/*     return GEN->iv; */
/*   } */

/*   /\* center of x-interval as splitting point *\/ */
/*   p_new = 0.5 * (iv->next->p + iv->p); */

/*   /\* we do not split an interval if is too close *\/ */
/*   /\*  changing the below FP_equal to FP_same can strongly increase the number of */
/*       intervals needed and may slightly decrease the MAError. In both cases the */
/*       required u-precision is not reached due to numerical problems with very steep CDF*\/ */
/*   if (_unur_FP_equal(p_new,iv->p) || _unur_FP_equal(p_new,iv->next->p)) { */
/*     if(!(*error_count_shortinterval)){  */
/*       _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF, */
/* 		     "one or more intervals very short; possibly due to numerical problems with a pole or very flat tail"); */
/*       (*error_count_shortinterval)++; */
/*     }  */
/*     /\* skip to next interval *\/ */
/*     _unur_pinv_interval_parameter(gen,iv); */
/*     return iv->next; */
/*   } */

/*   /\* 3rd check: |u_i - u_{i-1}| must not exceed threshold value *\/ */
/*   /\* 4th check: monotonicity                                    *\/ */

/*   if ( (iv->next->u - iv->u > PINV_MAX_U_LENGTH) || */
/*        (! _unur_pinv_interval_is_monotone(gen,iv)) ) { */
/*     /\* insert new interval into linked list *\/ */
/*     iv_new = _unur_pinv_interval_new(gen,p_new,CDF(p_new)); */
/*     if (iv_new == NULL) return NULL; */
/*     iv_new->next = iv->next; */
/*     iv->next = iv_new; */
/*     return iv; */
/*   } */

/*   /\* compute coefficients for spline (only necessary if monotone) *\/ */
/*   _unur_pinv_interval_parameter(gen,iv); */

/*   /\* 5th check: error in u-direction *\/ */

/*   /\* compute approximate value for inverse CDF in center of interval *\/ */
/*   x = _unur_pinv_eval_polynomial( 0.5, iv->spline, GEN->order ); */
/*   Fx = CDF(x); */

/*   /\* check for FP errors *\/ */
/*   if (_unur_isnan(x)) {  */
/*     _unur_error(gen->genid,UNUR_ERR_ROUNDOFF, */
/*  		"NaN occured; possibly due to numerical problems with a pole or very flat tail"); */
/*     return NULL; */
/*    } */

/*   /\* check error *\/ */
/*   if (!(fabs(Fx - 0.5*(iv->next->u + iv->u)) < GEN->u_resolution)) { */
/*     /\* error in u-direction too large *\/ */
/*     /\* if possible we use the point x instead of p_new *\/ */
/*     if(fabs(p_new-x)< PINV_XDEVIATION * (iv->next->p - iv->p)) */
/*       iv_new = _unur_pinv_interval_new(gen,x,Fx); */
/*     else */
/*       iv_new = _unur_pinv_interval_new(gen,p_new,CDF(p_new)); */
/*     if (iv_new == NULL) return NULL; */
/*     iv_new->next = iv->next; */
/*     iv->next = iv_new; */
/*     return iv; */
/*   } */

/*   /\* interval o.k. *\/ */
/*   return iv->next; */

/* } /\* end of _unur_pinv_interval_adapt() *\/ */

/*---------------------------------------------------------------------------*/

/* int  */
/* _unur_pinv_interval_is_monotone( struct unur_gen *gen, struct unur_pinv_interval *iv ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* check whether the given interval is monotone.                        *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*   iv  ... pointer to interval                                        *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   TRUE  ... if monotone                                              *\/ */
/*      /\*   FALSE ... otherwise                                                *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   double bound; */

/*   switch (GEN->order) { */
/*   case 5: */
/*     /\* the monotone check is in the moment only implemented for order 3 */
/*        as approximation we use the same check for order 5 *\/ */
/*   case 3: */
/*     /\* we skip the test if computing the bound has too many round-off errors *\/ */
/*     if (_unur_iszero(iv->u) || _unur_FP_approx(iv->u,iv->next->u)) */
/*       return TRUE; */
/*     /\* difference quotient *\/ */
/*     bound = 3.*(iv->next->p - iv->p)/(iv->next->u - iv->u); */
/*     return (1./iv->next->f > bound || 1./iv->f > bound) ? FALSE : TRUE; */
/*   case 1: */
/*     /\* linear interpolation is always monotone *\/ */
/*   default:  /\* we assume that we have checked GEN->order very often till now *\/ */
/*     return TRUE; */
/*   } */

/* } /\* end of _unur_pinv_interval_is_monotone() *\/ */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_interval_parameter( struct unur_gen *gen, struct unur_pinv_interval *iv ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* compute all parameter for interval (spline coefficients).            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*   iv  ... pointer to interval                                        *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   double delta_u, delta_p; */
/*   double f1, fs0, fs1, fss0, fss1; */

/*   delta_u = iv->next->u - iv->u; */
/*   delta_p = iv->next->p - iv->p; */

/*   switch (GEN->order) { */

/*   case 5:    /\* quintic Hermite interpolation *\/ */
/*     if (iv->f > 0. && iv->next->f > 0. && */
/* 	iv->df < INFINITY && iv->df > -INFINITY &&  */
/* 	iv->next->df < INFINITY && iv->next->df > -INFINITY ) { */
/*       f1   = delta_p; */
/*       fs0  = delta_u / iv->f;       */
/*       fs1  = delta_u / iv->next->f; */
/*       fss0 = -delta_u * delta_u * iv->df / (iv->f * iv->f * iv->f); */
/*       fss1 = -delta_u * delta_u * iv->next->df / (iv->next->f * iv->next->f * iv->next->f); */
      
/*       iv->spline[0] = iv->p; */
/*       iv->spline[1] = fs0; */
/*       iv->spline[2] = 0.5*fss0; */
/*       iv->spline[3] = 10.*f1 - 6.*fs0 - 4.*fs1 - 1.5*fss0 + 0.5*fss1; */
/*       iv->spline[4] = -15.*f1 + 8.*fs0 + 7.*fs1 + 1.5*fss0 - fss1; */
/*       iv->spline[5] = 6.*f1 - 3.*fs0 - 3.*fs1 - 0.5*fss0 + 0.5*fss1; */
/*       return UNUR_SUCCESS; */
/*     } */
/*     else { */
/*       /\* cannot use quintic interpolation in interval; use cubic instead *\/ */
/*       iv->spline[4] = 0.; */
/*       iv->spline[5] = 0.; */
/*     } */

/*   case 3:    /\* cubic Hermite interpolation *\/ */
/*     if (iv->f > 0. && iv->next->f > 0.) { */
/*       iv->spline[0] = iv->p; */
/*       iv->spline[1] = delta_u / iv->f; */
/*       iv->spline[2] = 3.* delta_p - delta_u * (2./iv->f + 1./iv->next->f); */
/*       iv->spline[3] = -2.* delta_p + delta_u * (1./iv->f + 1./iv->next->f); */
/*       return UNUR_SUCCESS; */
/*     } */
/*     else { */
/*       /\* cannot use cubic interpolation in interval; use linear instead *\/ */
/*       iv->spline[2] = 0.; */
/*       iv->spline[3] = 0.; */
/*     } */

/*   case 1:    /\* linear interpolation *\/ */
/*     iv->spline[0] = iv->p; */
/*     iv->spline[1] = delta_p; */
/*     return UNUR_SUCCESS; */

/*   default: */
/*     _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,""); */
/*     return UNUR_ERR_SHOULD_NOT_HAPPEN; */
/*   } */

/* } /\* end of _unur_pinv_interval_parameter() *\/ */

/*---------------------------------------------------------------------------*/

/* double */
/* _unur_pinv_eval_polynomial( double x, double *coeff, int order ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* evaluate polynomial using Horner scheme.                             *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   x     ... argument                                                 *\/ */
/*      /\*   coeff ... coefficients of polynomial (increasing order)            *\/ */
/*      /\*   order ... order of polynomial                                      *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   value of spline at u                                               *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i; */
/*   double poly; */

/*   poly = coeff[order]; */
/*   for (i=order-1; i>=0; i--) */
/*     poly = x*poly + coeff[i]; */

/*   return poly; */
/* } /\* end of _unur_pinv_eval_polynomial() *\/ */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_list_to_array( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* copy list of intervals into double array.                            *\/ */
/*      /\* the linked list is freed.                                            *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i;  */
/*   struct unur_pinv_interval *iv, *next; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE); */

/*   /\* allocate memory *\/ */
/*   GEN->intervals =  */
/*     _unur_xrealloc( GEN->intervals, GEN->N*(GEN->order+2)*sizeof(double) ); */

/*   i = 0; */
/*   for (iv=GEN->iv; iv!=NULL; iv=next) { */
/*     /\* copy *\/ */
/*     GEN->intervals[i] = iv->u; */
/*     memcpy( GEN->intervals+(i+1), &(iv->spline), (GEN->order+1)*sizeof(double) ); */
/*     i += GEN->order+2; */
/*     /\* and free linked list *\/ */
/*     next = iv->next; */
/*     free(iv); */
/*   } */

/*   /\* linked list is now empty *\/ */
/*   GEN->iv = NULL; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */
/* } /\* end of _unur_pinv_list_to_array() *\/ */

/*---------------------------------------------------------------------------*/

/* int */
/* _unur_pinv_make_guide_table( struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* make a guide table for indexed search                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* return:                                                              *\/ */
/*      /\*   UNUR_SUCCESS ... on success                                        *\/ */
/*      /\*   error code   ... on error                                          *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i,j, imax; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE); */

/*   /\* allocate blocks for guide table (if necessary). */
/*      (we allocate blocks for maximal guide table.) *\/ */
/*   GEN->guide_size = GEN->N * PINV_GUIDE_FACTOR; */
/*   if (GEN->guide_size <= 0) GEN->guide_size = 1;  */
/*   GEN->guide = _unur_xrealloc( GEN->guide, GEN->guide_size * sizeof(int) ); */

/*   imax = (GEN->N-2) * (GEN->order+2); */

/*   /\* u value at end of interval *\/ */
/* # define u(i)  (GEN->intervals[(i)+GEN->order+2]) */

/*   i = 0; */
/*   GEN->guide[0] = 0; */
/*   for( j=1; j<GEN->guide_size ;j++ ) { */
/*     while( u(i) < (j/(double)GEN->guide_size) && i <= imax) */
/*       i += GEN->order+2; */
/*     if (i > imax) break; */
/*     GEN->guide[j]=i; */
/*   } */

/* # undef u */

/*   /\* check i *\/ */
/*   i = _unur_min(i,imax); */

/*   /\* if there has been an round off error, we have to complete the guide table *\/ */
/*   for( ; j<GEN->guide_size ;j++ ) */
/*     GEN->guide[j] = i; */

/*   /\* o.k. *\/ */
/*   return UNUR_SUCCESS; */
/* } /\* end of _unur_pinv_make_guide_table() *\/ */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
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
/*   int i; */

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = PINV (Hermite approximation of INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_pinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

/*   fprintf(log,"%s: order of polynomial = %d",gen->genid,GEN->order); */
/*   _unur_print_if_default(gen,PINV_SET_ORDER); */

/*   fprintf(log,"\n%s: u-resolution = %g",gen->genid,GEN->u_resolution); */
/*   _unur_print_if_default(gen,PINV_SET_U_RESOLUTION); */
/*   fprintf(log,"\n%s: tail cut-off points = ",gen->genid); */
/*   if (GEN->tailcutoff_left < 0.)  fprintf(log,"none, "); */
/*   else                            fprintf(log,"%g, ",GEN->tailcutoff_left); */
/*   if (GEN->tailcutoff_right > 1.) fprintf(log,"none\n"); */
/*   else                            fprintf(log,"1.-%g\n",1.-GEN->tailcutoff_right); */
  
/*   fprintf(log,"%s: domain of computation = [%g,%g]\n",gen->genid,GEN->bleft,GEN->bright); */
/*   fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax); */
/*   fprintf(log,"%s:\n",gen->genid); */

/*   if (GEN->stp && gen->set & PINV_SET_STP) { */
/*     fprintf(log,"%s: starting points: (%d)",gen->genid,GEN->n_stp); */
/*     for (i=0; i<GEN->n_stp; i++) { */
/*       if (i%5==0) fprintf(log,"\n%s:\t",gen->genid); */
/*       fprintf(log,"   %#g,",GEN->stp[i]); */
/*     } */
/*   fprintf(log,"\n%s:\n",gen->genid); */
/*   } */
 
/*   fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid); */
/*   fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PINV_GUIDE_FACTOR); */
/*   _unur_print_if_default(gen,PINV_SET_GUIDEFACTOR); */
/*   fprintf(log,"\n%s:\n",gen->genid); */

/*   _unur_pinv_debug_intervals(gen); */

  fprintf(log,"%s: initialization %s\n",gen->genid,((ok)?"successful":"failed")); 
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_pinv_debug_init() */

/*---------------------------------------------------------------------------*/

/* void */
/* _unur_pinv_debug_intervals( const struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* print table of intervals into logfile                                *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   int i,n; */
/*   FILE *log; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID); */

/*   log = unur_get_stream(); */

/*   fprintf(log,"%s: Intervals: %d\n",gen->genid,GEN->N-1); */

/*   if (gen->debug & PINV_DEBUG_TABLE) { */
/*     fprintf(log,"%s:   Nr.      u=CDF(p)     p=spline[0]   spline[1]    ...\n",gen->genid); */
/*     for (n=0; n<GEN->N-1; n++) { */
/*       i = n*(GEN->order+2); */
/*       fprintf(log,"%s:[%4d]: %#12.6g  %#12.6g  %#12.6g", gen->genid, n, */
/* 	      GEN->intervals[i], GEN->intervals[i+1], GEN->intervals[i+2]); */
/*       if (GEN->order>1) */
/* 	fprintf(log,"  %#12.6g  %#12.6g", GEN->intervals[i+3], GEN->intervals[i+4]); */
/*       if (GEN->order>3) */
/* 	fprintf(log,"  %#12.6g  %#12.6g", GEN->intervals[i+5], GEN->intervals[i+6]); */
/*       fprintf(log,"\n"); */
/*     } */
/*     /\* the following might cause troubles when creating the tables fails. *\/ */
/*     /\* so we remove it.                                                   *\/ */
/*     /\*     i = n*(GEN->order+2); *\/ */
/*     /\*     fprintf(log,"%s:[%4d]: %#12.6g  %#12.6g  (right boundary)\n", gen->genid, n, *\/ */
/*     /\* 	    GEN->intervals[i], GEN->intervals[i+1] ); *\/ */
/*   } */

/*   fprintf(log,"%s:\n",gen->genid); */

/* } /\* end of _unur_pinv_debug_intervals() *\/ */

/*---------------------------------------------------------------------------*/

/* void  */
/* _unur_pinv_debug_chg_truncated( const struct unur_gen *gen ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* print new (changed) domain of (truncated) distribution               *\/ */
/*      /\*                                                                      *\/ */
/*      /\* parameters:                                                          *\/ */
/*      /\*   gen ... pointer to generator object                                *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   FILE *log; */

/*   /\* check arguments *\/ */
/*   CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID); */

/*   log = unur_get_stream(); */

/*   fprintf(log,"%s: domain of (truncated) distribution changed:\n",gen->genid); */
/*   fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]); */
/*   fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax); */

/* } /\* end of _unur_pinv_debug_chg_truncated() *\/ */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

/* void */
/* _unur_pinv_info( struct unur_gen *gen, int help ) */
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

/*   /\* generator ID *\/ */
/*   _unur_string_append(info,"generator ID: %s\n\n", gen->genid); */
  
/*   /\* distribution *\/ */
/*   _unur_string_append(info,"distribution:\n"); */
/*   _unur_distr_info_typename(gen); */
/*   _unur_string_append(info,"   functions = CDF"); */
/*   if (GEN->order > 1) */
/*     _unur_string_append(info," PDF"); */
/*   if (GEN->order > 3) */
/*     _unur_string_append(info," dPDF"); */
/*   _unur_string_append(info,"\n"); */
/*   _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]); */
/*   if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) { */
/*     _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]); */
/*   } */
/*   _unur_string_append(info,"\n"); */

/*   if (distr->set & UNUR_DISTR_SET_MODE) { */
/*     _unur_string_append(info,"   mode      = %g\n", DISTR.mode); */
/*   } */

/*   if (help) */
/*     if (! (distr->set & UNUR_DISTR_SET_MODE) ) */
/*       _unur_string_append(info,"\n[ Hint: %s ]\n", */
/* 			  "You may set the \"mode\" of the distribution in case of a high peak"); */

/*   _unur_string_append(info,"\n"); */
      
/*   /\* method *\/ */
/*   _unur_string_append(info,"method: PINV (Hermite approximation of INVerse CDF)\n"); */
/*   _unur_string_append(info,"   order of polynomial = %d\n", GEN->order); */
/*   _unur_string_append(info,"\n"); */

/*   /\* performance *\/ */
/*   _unur_string_append(info,"performance characteristics:\n"); */
/*   _unur_string_append(info,"   truncated domain = (%g,%g)\n",GEN->bleft,GEN->bright); */
/*   _unur_string_append(info,"   Prob(X<domain)   = %g\n", _unur_max(0,GEN->tailcutoff_left)); */
/*   _unur_string_append(info,"   Prob(X>domain)   = %g\n", _unur_max(0,1.-GEN->tailcutoff_right)); */
/*   { */
/*     double max_error=1.; double MAE=1.; */
/*     unur_pinv_estimate_error( gen, 10000, &max_error, &MAE ); */
/*     _unur_string_append(info,"   u-error         <= %g  (mean = %g)\n", max_error, MAE); */
/*   } */

/*   _unur_string_append(info,"   # intervals      = %d\n", GEN->N-1); */
/*   _unur_string_append(info,"\n"); */
  

/*   /\* parameters *\/ */
/*   if (help) { */
/*     _unur_string_append(info,"parameters:\n"); */
/*     _unur_string_append(info,"   order = %d  %s\n", GEN->order, */
/*  			(gen->set & PINV_SET_ORDER) ? "" : "[default]"); */

/*     _unur_string_append(info,"   u_resolution = %g  %s\n", GEN->u_resolution, */
/*  			(gen->set & PINV_SET_U_RESOLUTION) ? "" : "[default]"); */
    
/*     _unur_string_append(info,"   boundary = (%g,%g)  %s\n", GEN->bleft, GEN->bright, */
/* 			(gen->set & PINV_SET_BOUNDARY) ? "" : "[computed]"); */

/*     _unur_string_append(info,"\n"); */
/*     /\* Not displayed: */
/*        int unur_pinv_set_cpoints( UNUR_PAR *parameters, const double *stp, int n_stp ); */
/*        int unur_pinv_set_guidefactor( UNUR_PAR *parameters, double factor ); */
/*     *\/ */
/*   } */


/*   /\* Hints *\/ */
/*   if (help) { */
/*     if ( GEN->order < 5 )  */
/*       _unur_string_append(info,"[ Hint: %s ]\n", */
/* 			  "You can set \"order=5\" to decrease #intervals"); */
/*     if (! (gen->set & PINV_SET_U_RESOLUTION) ) */
/*       _unur_string_append(info,"[ Hint: %s\n\t%s ]\n", */
/* 			  "You can decrease the u-error by decreasing \"u_resolution\".", */
/* 			  "(it is bounded by the machine epsilon, however.)"); */
/*     _unur_string_append(info,"\n"); */
/*   } */

/* } /\* end of _unur_tdr_info() *\/ */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/



int check_inversion_unuran(struct unur_gen *gen,double uerror,double (*cdf)(double x),int printyn){
  /* checks the inversion for a generator object by calculating the u-error using the CDF
     checks many u-values especially in the etreme tails
  */

  double u,x,uerr,maxerror=0.,maxu;

  uerror*=1.1;
/* "*1.1" is necessary for this controll as the left cut-off tail has about 10% of the uerror */
  maxu=1.;
  for(u=1.e-10;u<maxu;u+=1./33211){
    x = _unur_pinv_eval_approxinvcdf(gen,u);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror) printf("%g %g %g\n",u,x,uerr);
  }
  for(u=1.e-13;u<1.e-5;u+=1.e-8){//check left tail
    x = _unur_pinv_eval_approxinvcdf(gen,u);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror)printf("%g %g %g\n",u,x,uerr);
  }
  for(u=maxu-1.e-13;u>maxu-1.e-5;u-=1.e-8){//check right tail
    x = _unur_pinv_eval_approxinvcdf(gen,u);
    uerr=(*cdf)(x)-u;
    if(fabs(uerr)>maxerror) maxerror=uerr;
    if(printyn&&fabs(uerr)>uerror)printf("%g %g %g\n",u,x,uerr);
  }
  printf("aimed uerror: %g ;  observed maxerror : %g\n",uerror,maxerror);
  return 0;
}

/**********************************************/


















double lobato5 (struct unur_gen *gen, double x, double h, double fx, double *fxph)
/************************************
 * Numerical Integration of the interval (x,x+h)
 * using Gauss-Lobato integration with 5 points.
 * fx ... f(x) to save calls to f()
 * *fxph ... f(x+h) 
 ************************************/
#define W1 0.17267316464601146  //= 0.5-sqrt(3/28)
#define W2 (1.-W1)
{ double ifx,ifxph;

 ifxph = PDF(x+h);
 if(fxph!=NULL){ 
   ifx = fx;
   *fxph = ifxph;
 }
 else {
   ifx = PDF(x);
 }
/*  return (9*(ifx+(ifxph))+49.*((*f)(x+h*W1)+(*f)(x+h*W2))+64.*(*f)(x+h/2.))*h/180.; */
 return (9*(ifx+(ifxph))+49.*(PDF(x+h*W1)+PDF(x+h*W2))+64.*PDF(x+h/2.))*h/180.;
}
#undef W1
#undef W2


double nint_12 (struct unur_gen *gen, double a,double b,double *res_relerror){
  //with 1 and two intervals reports the resulting (relative-error) in res_relerror 
  int i;
  double res,reso;

  reso = lobato5(gen,a,(b-a),0,NULL);
  i=1;
  res = lobato5(gen,a,(b-a)*0.5,0,NULL)+lobato5(gen,(b+a)*0.5,(b-a)*0.5,0,NULL);
  *res_relerror = fabs(res-reso)/res; 
  return res;
}

double nint_monoton_dens(struct unur_gen *gen, double a,double b,double step,double crit){
  // step <= b-1 !!!
// numerical integration with a variable step-size for a monoton density, starting with "step"
// crit ... maximal accepted relative error

//TODO: Kann man verbessern, wenn man in den Tails den Error relativ zum Gesamtintegral rechnet

  double sum=0.,sumi,x,error;
  int i=1;

  x=a;
  while(x<b){
    sumi = nint_12(gen,x,x+step,&error);
    //printf("%d: x%g step%g sumi%g error%g\n",i,x,step,sumi,error);
    if(error > crit){ 
      step = step*pow(0.5*crit/error,1./9.);
    }else{
      i++;
      sum += sumi;
      x = x+step;
      //      step = step*1.2*pow(crit/error,1./9.);
      step = step*pow(crit/error,1./9.);
      if(x+step >b){
        step = (b-x)*(1+5.e-16);
      }
    }    
  }
  return sum;
}


////////////////////////////////////////////////////////////////////////////////////

double evalnewtoninterpol(double u,int g,double ui[],double zi[]){
  /* zi = pi;  ui=yi */
  double q,chi;
  int k;
  
  q = u;
  chi = zi[g];
  for(k=g-1;k>=1;k--) chi = chi*(q-ui[k])+zi[k];
  return( chi*q);
}


int newtoninterpol (struct unur_gen *gen, double x0, double h,double ui[],double zi[],double *x)
{
  /*calculates ui and zi arrays, xi pointer may be NULL */
  double xi, dxi, temp,// zi[20], ui[20]={0.},
    phi;/*20 statt g+1 */
  int i,k;

  phi = M_PI*0.5/(GEN->order+1);
  ui[0]=0.;
  if(x!=NULL)x[0]=xi;
  for(i=1; i<=GEN->order; i++){
    xi = x0 + h*sin((i-1)*phi)*sin(i*phi)/cos(phi);
    dxi = h*sin(2*i*phi)*tan(phi);
    if(x!=NULL)x[i]=xi+dxi;
   temp = lobato5(gen,xi, dxi,0,NULL);
   if(temp<1.e-50){
     printf("ERROR!! Newtoninterpolation interval too short. or density 0 EXITINGNG\n");
     exit(1);
   }
    zi[i]=dxi/temp;
    ui[i]=ui[i-1]+temp; /* ui[0] initialisation??? */
  } // ui[g] is the probability of the interval
  for(k=2; k<=GEN->order; k++)
    for(i=GEN->order; i>=k; i--)
      zi[i]=(zi[i]-zi[i-1])/(ui[i]-ui[i-k]);

  //  printf("x:");
  //    printvec(g+1,x);
  //printf("u:");
  //  printvec(g+1,ui);

  return 1;
}


int tstpt(int g,double ui[],double utest[]){
  int k,j,i;
  double sum, qsum,x;
  utest[0]=0.;
  for(k=1; k<=g; k++){
    x = 0.5*(ui[k-1]+ui[k]);
    for(j=1; j<=2; j++){
      sum = 1./x;
      qsum = sum*sum;
      for(i=1; i<=g; i++){
	sum = sum + 1./(x-ui[i]);
	qsum = qsum + 1./((x-ui[i])*(x-ui[i]));
      }
      x +=sum/qsum;
    }
    utest[k] = x;
  }
  return 1;
}

double maxerrornewton (struct unur_gen *gen, double ui[],double zi[],double x0,double xval[]){
  double maxerror=0.,uerror,*testu,uarr[21],x;
  int i,n;
    n=GEN->order;
    testu=uarr;
    tstpt(GEN->order,ui,testu);
    //printvec(n+1,testu);
    //printvec(n+1,ui);
    //if(xval!=NULL)printvec(3,xval);
  for(i=0;i<n;i++){
    //Naechste Zeile: TODO Verbesserung moeglich, wenn man die xi verwendet und so kuerzere
    // intervalle fuer LObato integration bekommt
    x=evalnewtoninterpol(testu[i+1],GEN->order,ui,zi);
    if(i==0||xval==NULL) uerror=fabs(lobato5(gen,x0,x ,0.,NULL)-testu[i+1]);
    else uerror=fabs(ui[i]+lobato5(gen,xval[i],x+x0-xval[i],0.,NULL)-testu[i+1]);
    if(uerror>maxerror) maxerror=uerror;
    //    printf("%d:u max %g %g\n",i,uerror,maxerror);
  }
  //printf("maxerror %g\n",maxerror);
  return maxerror;
}


int  setup(struct unur_gen *gen, double a, double b, double hh, double uerror){
  /*
    f ... PDF
    g ... order of polynomial
    a ...
    b ...
    hh ...
    uerror ... u-error
  */
  double maxerror,h=hh,*xval;
  int i,j,cont;
  int countextracalc=0;

  xval=malloc(sizeof(double)*(GEN->order+1));

  GEN->iv[0].ui = malloc(sizeof(double)*(GEN->order+1));
  GEN->iv[0].zi = malloc(sizeof(double)*(GEN->order+1));
  GEN->iv[0].xi = a;
  GEN->iv[0].cdfi = 0.;//cdfi holds cdf value at the left border of the interval
  cont=1;
  i=0;
  while(cont){
    if(GEN->iv[i].xi+h >b){
      h = b - GEN->iv[i].xi;
      cont=0;
    }
    newtoninterpol(gen,GEN->iv[i].xi,h,GEN->iv[i].ui,GEN->iv[i].zi,xval);
    maxerror = maxerrornewton(gen,GEN->iv[i].ui,GEN->iv[i].zi,GEN->iv[i].xi,xval);
    if(maxerror > uerror){ 
      countextracalc++;
      h*= 0.9;
     if(maxerror>4.*uerror) h*=0.9;
 
    }else{
	GEN->iv[i+1].ui = malloc(sizeof(double)*(GEN->order+1));
  	GEN->iv[i+1].zi = malloc(sizeof(double)*(GEN->order+1));
  	GEN->iv[i+1].xi = GEN->iv[i].xi+h;
  	GEN->iv[i+1].cdfi = GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order];//cdfi holds cdf value at the left border of the interval
        if(maxerror < 0.3*uerror) h*=1.2;
        if(maxerror < 0.1*uerror) h*=2.;
        i++;
   }
   if(i>PINVMAXINT){
     printf("error setup(); i>maxint; EXITING\n");
     exit(1);
   }
 }
 GEN->ni = i;
 GEN->iv = realloc(GEN->iv,sizeof(struct siv)*(GEN->ni+1));

 free(xval);

 GEN->umax = GEN->iv[GEN->ni].cdfi;


 GEN->guide_size = GEN->ni;//size of guide-table
 GEN->guide = malloc(sizeof(int)*GEN->guide_size);  
 GEN->guide[0] = 0;
 i=0;
 for(j=1; j<GEN->guide_size; j++){
   while(j/(double)GEN->guide_size > GEN->iv[i+1].cdfi/GEN->umax) i++;
   /* "/geno->umax" above is necessary, as we need the guide table for u in (0,umax) */
   GEN->guide[j] = i;
 }
 printf("Set-up finished: g=%d,  Number of intervals = %d,\n         additional calculated interpolations=%d\n",
	GEN->order,GEN->ni,countextracalc);
 printf("u in (0,%.18g)   1-umax%g\n",GEN->umax,1-GEN->umax);

 return UNUR_SUCCESS;
} 



double searchborder (struct unur_gen *gen, double x0, double step,double border){
  // x0 starting point with f(x0) not small; may but need not be the mode
  // step first step size, includes direction
  double fx0=PDF(x0),fx=fx0,x=x0,xa;
  int i;

  for(i=0;i<100 && PDF(x)>fx0*1.e-13;i++){
    xa=x;
    if(i>10) step*=2;
    if((x-border)*(x+step-border)>0.) x+= step;
    else x=(x+border)*0.5;
    fx=PDF(x);
    //printf("%d:x%g fx%g\n",i,x,fx);
  }

  do{
    x=(x+xa)*0.5;
    fx=PDF(x);
    //printf("%d:x%g fx%g\n",i,x,fx);
  }while(fx<fx0*1.e-13);
 


  return x;
}

/************************************/
double tail (struct unur_gen *gen, double x, double d){
  /********************************
     calculates approximate tail area = f(x)^2/((lc_f(x)+1)*abs(f'(x)) 
     x...cut off point
     d... step length for numeric differentiation
  **********************************/
  double ff,fp,fm,cplus1;
  ff = PDF(x);
  fp = PDF(x+d);
  fm = PDF(x-d);
  if(fm-2.*ff+fp < 0.||fm<1.e-100||ff<1.e-100||fp<1.e-100){
    printf("warning possible problem in function tail() !!!\n");
  }
  cplus1 = fp/(fp-ff) + fm/(fm-ff);
  return ff*ff/(cplus1*fabs((fp-fm)/(2.*d)));
}




double cut (struct unur_gen *gen, double w,double dw, double crit){
  /**********************
    calculates starting from w the left (dw<0) or right (dw>0) cut off point.
    crit... u-error criterium for tail cut off 
    dw ... initial step-length
    the are outside is approximately crit/10
  ****************/

  double H,rezH,y,d,rezy,yplus,rezys,corr;
  int j,k,cont;

  H = crit;
  rezH=1./H;
  cont=1;
  for(j=1;j<1000&&cont;j++){
    y=tail(gen,w,dw/64.);
    if(y<H) cont=0;
    else{
      w+=dw;
      if(j>32) dw*=1.5;
    }
  }
  if(cont){
    printf("error in cut(), first loop; exiting!!!\n");
    exit(1);
  }

  d=dw/64.;
  for(k=0;k<50;k++){
    rezy=1./y;
    yplus = tail(gen,w+d,d);
    if(yplus<0){
      printf("error in cut(), yplus negative; exiting!!!\n");
      exit(1);
    }
    rezys = (1./yplus-rezy)/d;
    corr = -(rezy-rezH)/rezys;
    if(fabs(rezy/rezH-1.)<1.e-7) return w;
    w+=corr;
    y=tail(gen,w,d);
    if(y<0){
      printf("error in cut(), y negative; exiting!!!\n");
      exit(1);
    }
  }
    printf("error in cut(), second loop; exiting!!!\n");
    exit(1);


}

/**********************************************************************/




/******************************************************/

int pinvsetup( struct unur_gen *gen)
/***********************************************
  starts the set-up and returns a pointer to the generator object

  f      ... pointer to the pdf)
  g      ... degree of polynomial
  uerror ... maximal accepted u-error
  x0     ... value with f(x0) not small
  asearch... 1 ... search for a cut off value >= a
  bsearch... 1 ... search for a cut off value <= b
  a   ...left domain border, can be set very small eg. -1.e100 for asearch==1
  b   ...right domain border, can be set very large eg. 1.e100 for bsearch==1

  asearch==0 means that the direct a-value is left unchanged
  bsearch==0 means that the direct a-value is left unchanged
  note that this may lead to numericl problems in the set-up
  if f(x) is very small close to the border.

 So the algorithm should run more stable for  asearch and bsearch =1.
 
 asearch=0 is useful eg. for the Gamma(2) distribution where the left border 0
 is fixed and should not be searched for, as f(x) has no left tail

 ***********************************************/
{
/*   double a = GEN->bleft; */
/*   double b = GEN->bright; */

  double area,areal,tailcutfact;

  if(GEN->sleft) GEN->bleft = searchborder(gen,DISTR.center, -1, GEN->bleft);
  if(GEN->sright) GEN->bright = searchborder(gen,DISTR.center, 1, GEN->bright);
  area = nint_monoton_dens(gen,DISTR.center,GEN->bright,1.,1.e-8);
  areal = nint_monoton_dens(gen,GEN->bleft,DISTR.center,1.,1.e-8);
  area+=areal;
  printf("after searchborder: a=%g  b=%g area=%g\n",GEN->bleft,GEN->bright,area);
//  printf("a=%g  b=%g area=%g integralerror%g\n",a,b,area,area/(cdf(b)-cdf(a))-1.);

 
  tailcutfact= 0.1;
  if(GEN->u_resolution<=9.e-13) tailcutfact=0.5;
/* above command necessary for Cauchy distribution where cut has problems with very small values*/
 
  if(GEN->sleft) GEN->bleft = cut(gen,GEN->bleft,(GEN->bleft-GEN->bright)/128,GEN->u_resolution*area*tailcutfact);
  if(GEN->sright) GEN->bright = cut(gen,GEN->bright,(GEN->bright-GEN->bleft)/128,GEN->u_resolution*area*tailcutfact);
  printf("after cut: a=%g b=%g\n",GEN->bleft,GEN->bright);

  return setup(gen,GEN->bleft,GEN->bright,(GEN->bright-GEN->bleft)/128,GEN->u_resolution*area);
}

/****************************************/
