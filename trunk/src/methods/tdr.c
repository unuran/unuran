/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and .... of a T-concave distribution                     *
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
 *   REFERENCES:                                                             *
 *   [1] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *                                                                           *
 *   [2] Chen, H. C. and Asau, Y. (1974): On generating random variates      *
 *       from an empirical distribution, AIIE Trans. 6, pp. 163-166          *
 *                                                                           *
 *   [3] Gilks, W. R. and Wild,  P. (1992): Adaptive rejection sampling      *
 *       for Gibbs sampling, Applied Statistics 41, pp. 337-348              *
 *                                                                           *
 *   [4] Evans, M. and Swartz, T. (1998): Random variable generation         *
 *       using concavity properties of transformed densities,                *
 *       Journal of Computational and Graphical Statistics 7(4), pp. 514-528 *
 *                                                                           *
 *   [5] Derflinger, G. and Hoermann, W. (1999): The optimal selection of    *
 *       hat functions for rejection algorithms, preprint                    *
 *                                                                           *
 *   [6] Leydold J. (1999): Automatic Sampling with the ratio-of-uniforms    *
 *       method, preprint, 15pp.                                             *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 * Transformed density rejection (see [1,3,4]) is an acceptance/rejection    *
 * technique that uses the fact, that the probability density function f(x)  *
 * for many distribution is T-concave, to construct a hat function.          *
 * That is, there exists a transformation T(x), such that T(f(x)) is         *
 * concave. Then it is easy to construct a majorizing or hat function Th(x)  *
 * for the transformed density T(F(x)) by the pointwise minima of several    *
 * tangents. Transforming this back into the original scale gives the        *
 * hat h(x) = T^(-1)(Th(x)). Squeezes can be constructed by secants.         *
 *                                                                           *
 * However, transformation T(x) has to satisfy (see [1]):                    *
 *                                                                           *
 *  (1) T(f(x)) is concave.                                                  *
 *  (2) T is differentiable.                                                 *
 *  (3) T is strictly monotonically increasing (i.e., T'(x) > 0),            *
 *      which implies that T^(-1) exists.                                    *
 *  (4) the area below the hat is finite.                                    *
 *  (5) it is easy to generate from the hat distribution.                    *
 *                                                                           *
 * We use a family T_c(x) of distribution defined by (see [1]):              *
 *                                                                           *
 *           c    |   T_c(x)                                                 *
 *       ---------+---------------                                           *
 *           0    |   log(x)                                                 *
 *         -1/2   |   -1/sqrt(x)                                             *
 *        -1<c<0  |   -x^c                                                   *
 *                                                                           *
 * Figures 1 and 2 show the situation for the standard normal distribution   *
 * with T(x) = log(x) and two construction points at -1 and 0.3 for the      *
 * tangents.                                                                 *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *           figure 1: transformed density with hat and squeeze              *
 *                     T(x)   = log(x)                                       *
 *                     pdf(x) = exp(-x^2)                                    *
 *                                                                           *
 *                                .                                          *
 *     -3        -2        -1     ....0         1         2         3        *
 *    --+---------+---------X----.--*****-X-----+---------+---------+--->    *
 *                              .**     / **.                                *
 *                             *     /       *...                            *
 *                            *   /           *  ....                        *
 *             construction  * / squeeze       *     ...                     *
 *                   point  X                   *       ...                  *
 *                         *                     *         ....              *
 *                        *                       *            ...           *
 *                       *                         *              ..         *
 *                      .                                           .        *
 *                     .*                           *                        *
 *                    .*                             *                       *
 *                   .                                                       *
 *                  . *                               *                      *
 *                 . *                                 *                     *
 *                .                                                          *
 *               .  * transformed p.d.f.                *                    *
 *  transformed .                                                            *
 *         hat .   *                                     *                   *
 *            .                                                              *
 *           .    *                                       *                  *
 *          .                                                                *
 *         .     *                                         *                 *
 *        .                                                                  *
 *       .      *                                           *                *
 *      .                                                                    *
 *             *                                             *               *
 *                                                                           *
 *                                                                           *
 *            *                                               *              *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *                   figure 2: density with hat and squeeze                  *
 *                             T(x)   = log(x)                               *
 *                             pdf(x) = exp(-x^2)                            *
 *                                                                           *
 *                                ..                                         *
 *                                . .                                        *
 *                                .  .                                       *
 *                               .    .                                      *
 *                               .   ***.                                    *
 *                               . **   **                                   *
 *                              . *       X                                  *
 *                              .*       /|*                                 *
 *                              .*      / |* .                               *
 *                             .*      /  | * .                              *
 *                             .*      /  | *  .                             *
 *                             *      /   |  *  ..                           *
 *                             *     /    |  *    .                          *
 *                            *    /      |   *    ..                        *
 *                           *   /squeeze |    *     .                       *
 *         construction      * /          |    *      ..                     *
 *                point --- X/            |     *       ..                   *
 *                         *|             |      *        ..                 *
 *                         *|             |      *          ..               *
 *                        * |             |       *           ..             *
 *                        * |             |       *             ...          *
 *                       *  |             |        *               ..        *
 *                      .*  |             |        *                         *
 *                     .*   |             |         *                        *
 *                    .*    |             |          *                       *
 *                   . *    |             |           *                      *
 *                  . *     |             |           *                      *
 *           hat  .. *      |             |            *                     *
 *               .  *       |             |             *                    *
 *             ..  * p.d.f. |             |              *                   *
 *          ...   *         |             |               *                  *
 *       ...    **          |             |                **                *
 *      .     **            |             |                  **              *
 *        ****              |             |                    ****          *
 *    --**--------+---------X---------+---X-----+---------+--------**--->    *
 *     -3        -2        -1         0         1         2         3        *
 *                                                                           *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * To generate from the hat distribution we have to partition the domain     *
 * of the p.d.f. by means of the construction points of the tangents into    *
 * several intervals. Each interval has to be divided into two parts         *
 * according to the tangent line that defines the hat function.              *
 * Although the interval can be divided at any point in the interval,        *
 * the best possible splitting point obviously is the intersection point     *
 * of respective the tangents at the left and right construction point of    *
 * the interval. However there are two cases when this choice is not         *
 * appropriate: (1) The computation of this point is unstable. Then we use   *
 * the mean point of the interval. (2) One of the two construction points    *
 * is very much greater (or smaller) than the other one. Since we generate   *
 * points by adding or subtracting a random number to the construction       *
 * point, it then happens, that we have the difference of two large numbers, *
 * which results in serious roundoff error. Then it is better to choose the  *
 * intersection point such that we have only one interval with the point     *
 * nearer to the origin as the only construction point.                      *
 * (This procedure is a little bit different from [1] or [3], but the author *
 * of this program finds the implemented version more convenient.) Then we   *
 * have to compute the area below the hat for each of the two parts of the   *
 * interval i (A_l^i for the left hand part, A_r^i for the right hand part)  *
 * for each interval. Then we have to sample from a discrete random variate  *
 * with probability vector proportional to (A_l^0,A_r^0; A_l^1,A_r^1; ...;   *
 * A_l^n,A_r^n) to get one of the intervals. We use indexed search (or       *
 * guide tables) to perform this task ([2], see also description of DIS).    *
 * Then we can sample from the hat distribution in this part of the          *
 * interval by inversion. Notice that we can use reuse the uniform random    *
 * number from the first step (i.e., sampling an interval) for this second   *
 * step without danger. The whole generation can be seen as inversion from   *
 * the hat distribution.                                                     *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Algorithm TDR                                                             *
 *                                                                           *
 * [Required]                                                                *
 * p.d.f. f(x), transformation T(x), construction points c_1,...,c_n         *
 *                                                                           *
 * [Setup]                                                                   *
 *  1: Construct hat h(x) and squeeze s(x).                                  *
 *  2: Foreach Interval (c_i,c_{i+1}) Do                                     *
 *  3:    Compute areas A_l^i, A_r^i, A_squeeze.                             *
 *                                                                           *
 * [Generate]                                                                *
 *  4: Generate I proportional to (A_l^1,A_r^1; ...).                        *
 *  5: Generate X with p.d.f. proportional to h|I.                           *
 *  6: Generate U ~ U(0,1).                                                  *
 *  7: If U * h(x) <= min(f(c_i),f(c_{i+1})) Return X.                       *
 *  8: If U * h(x) <= s(X) Return X.                                         *
 *  9: If U * h(x) <= f(X) Return X.                                         *
 * 10: Goto 4.                                                               *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * There are several alternatives for the choice of the construction points  *
 *                                                                           *
 * (1) Adaptive rejection sampling:                                          *
 *     Use two points on both sides of the mode of the p.d.f. (provided      *
 *     that the hat has finite area).                                        *
 *     Whenever the p.d.f. has to be evaluated in step 8, add a new          *
 *     construction point at X.                                              *
 *                                                                           *
 * (2) Optimal construction points:                                          *
 *     [5] have show an asymptotic formula for optimal construction points,  *
 *     i.e., with minimiza the ration Ahat/Asqueeze for a given number of    *
 *     points. Although this formula is only valid for infinitely many       *
 *     construction points, it works very well even for 3 or 4 points.       *
 *                                                                           *
 *     [1] gives a formula for three optimal construction points.            *
 *                                                                           *
 * (3) Empirical formulas:                                                   *
 *     [6] uses points c_i, such that arctan(c_i) are equidistributed.       *
 *     ("equiangular rule").                                                 *
 *     In most cases it has turned out as an acceptable good choice.         *
 *     (We made one change to [6]: If the mode is known, all points are      *
 *     moved, such the mode becomes the center of these points.)             *
 *                                                                           *
 * This implementation uses (3) to get good points for the start and adds    *
 * additional construction points by means of method (1).                    *
 * (The number of starting points and the maximum number of points can       *
 * be set for each generator.)                                               *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Implementation details:                                                   *
 *                                                                           *
 * (1) If possible the boundary points of the domain are used as             *
 *     construction points. (So do not include these in a list of starting   *
 *     points.) As a consequence 0 starting points are allowed.              *
 *                                                                           *
 * (2) To make use of the automatic generation of starting points call       *
 *     `unur_set_cpoints()� with the NULL pointer as the last argument       *
 *     and the number of construction points (besides the boundary points    *
 *     and the mode) as its third argument.                                  *
 *                                                                           *
 * (3) If the mode is given, we use a tangent with slope 0 at this point;    *
 *     even when this does not result in the best possible hat (e.g. for     *
 *     the exponential distribution).                                        *
 *                                                                           *
 * (4) Starting points and the mode outside the domain are ignored.          *
 *                                                                           *
 * (5) If the given starting points does not result in a hat with finite     *
 *     area, the program tries to find some proper additional construction   *
 *     points by splitting interval with infinite area. (Here we use the     *
 *     "equiangular rule" again.)                                            *
 *                                                                           *
 * (6) The use of the additional "fast" squeeze in step 7 is due to a        *
 *     suggestion of G. Derflinger. It speeds up the generation about 20 %   *
 *     when many construction points are used.                               *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TDR_VARFLAG_VERIFY     0x100u   /* flag for verifying mode          */
#define TDR_VARFLAG_USECENTER  0x200u   /* flag whether center is used as cpoint or not */
#define TDR_VARFLAG_USEMODE    0x400u   /* flag whether mode is used as cpoint or not */

#define TDR_VARMASK_T          0x003u   /* indicates transformation         */
#define TDR_VAR_T_SQRT         0x001u   /* T(x) = -1/sqrt(x)                */
#define TDR_VAR_T_LOG          0x002u   /* T(x) = log(x)                    */
#define TDR_VAR_T_POW          0x003u   /* T(x) = -x^c                      */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TDR_DEBUG_IV           0x00000010u
#define TDR_DEBUG_SPLIT        0x00010000u
#define TDR_DEBUG_SAMPLE       0x01000000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TDR_SET_CENTER         0x002u
#define TDR_SET_STP            0x004u
#define TDR_SET_N_STP          0x008u
#define TDR_SET_GUIDEFACTOR    0x010u
#define TDR_SET_C              0x020u
#define TDR_SET_MAX_SQHRATIO   0x040u
#define TDR_SET_MAX_IVS        0x080u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TDR"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

/* No reinit() cal                                                           */
/*  static int _unur_tdr_reinit( struct unur_gen *gen );                     */
/*---------------------------------------------------------------------------*/
/* Re-initialize (existing) generator.                                       */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_sample_log( struct unur_gen *generator );
static double _unur_tdr_sample_sqrt( struct unur_gen *generator );
static double _unur_tdr_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_get_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute intervals from given starting construction points.                */
/*---------------------------------------------------------------------------*/

static struct unur_tdr_interval *_unur_tdr_interval_new( struct unur_gen *gen, 
							 double x, double fx, int is_mode );
/*---------------------------------------------------------------------------*/
/* make a new segment with left construction point x.                        */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for interval.                                  */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_interval_split( struct unur_gen *gen, 
				      struct unur_tdr_interval *iv_old, double x, double fx );
/*---------------------------------------------------------------------------*/
/* split am interval point x. return 0 if not successful.                    */                                           
/*---------------------------------------------------------------------------*/

static int _unur_tdr_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_interval_division_point( struct unur_gen *gen,
					      struct unur_tdr_interval *iv, double *ipt );
/*---------------------------------------------------------------------------*/
/* compute cutting point of interval into left and right part.               */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv,
				       double slope, double x );
/*---------------------------------------------------------------------------*/
/* compute area below piece of hat or slope in                               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_intervals( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for intervals                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_sample( struct unur_gen *gen, struct unur_tdr_interval *iv, struct unur_tdr_interval *pt, 
				    double x, double fx, double hx, double sqx );
/*---------------------------------------------------------------------------*/
/* print data while sampling from generators.                                */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_split_start( struct unur_gen *gen, 
					 struct unur_tdr_interval *iv, double x, double fx );
static void _unur_tdr_debug_split_stop( struct unur_gen *gen, 
					struct unur_tdr_interval *iv_left, struct unur_tdr_interval *iv_right );
/*---------------------------------------------------------------------------*/
/* print before and after an interval has been split (not / successfully).   */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.tdr         /* data for parameter object         */
#define GEN       gen->data.tdr         /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),&(gen->distr))   /* call to p.d.f.      */
#define dPDF(x)   _unur_cont_dPDF((x),&(gen->distr))  /* call to derivative of p.d.f. */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_tdr_new( struct unur_distr* distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.d.f."); return NULL; }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of p.d.f."); return NULL; }

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_TDR_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR.guide_factor        = 3.;     /* size of guide table / number of intervals */

  PAR.c_T                 = -0.5;   /* parameter for transformation (-1. <= c < 0.) */

  PAR.starting_cpoints    = NULL;   /* pointer to array of starting points   */
  PAR.n_starting_cpoints  = 10;     /* number of starting points             */
  PAR.max_ivs             = 50;     /* maximum number of intervals           */
  PAR.max_ratio           = 0.95;   /* bound for ratio  Atotal / Asqueeze    */
  PAR.bound_for_adding    = 0.5;    /* do not add a new construction point in an interval,
				       where ambigous region is too small, i.e. if 
				       area / ((A_hat - A_squeeze)/number of segments) < bound_for_adding */
 
  par->method   = UNUR_METH_TDR;                 /* method                   */
  par->variant  = ( TDR_VARFLAG_USECENTER |      /* default variant          */
		    TDR_VARFLAG_USEMODE );

  par->set      = 0u;               /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* we use the mode (if known) as center of the distribution */
  if (distr->set & UNUR_DISTR_SET_MODE) {
    PAR.center = DISTR_IN.mode;
    par->set |= TDR_SET_CENTER;
  }
  else
    PAR.center = 0.;        /* the default */

  /* routine for starting generator */
  par->init = _unur_tdr_init;

  return par;

} /* end of unur_tdr_new() */

/*****************************************************************************/

int
unur_tdr_set_cpoints( struct unur_par *par, int n_stp, double *stp )
     /*----------------------------------------------------------------------*/
     /* set construction points for envelope                                 */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*   stp    ... pointer to array of starting points                     */
     /*              (NULL for changing only the number of default points)   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return 0;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return 0;
      }

  /* store date */
  PAR.starting_cpoints = stp;
  PAR.n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TDR_SET_N_STP | ((stp) ? TDR_SET_STP : 0);

  return 1;

} /* end of unur_tdr_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= TDR_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_tdr_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_max_sqhratio( struct unur_par *par, double max_ratio )
     /*----------------------------------------------------------------------*/
     /* set bound for ratio A(squeeze) / A(hat)                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ratio ... upper bound for ratio to add a new construction point*/
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1.+DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return 0;
  }

  /* store date */
  PAR.max_ratio = max_ratio;

  /* changelog */
  par->set |= TDR_SET_MAX_SQHRATIO;

  return 1;

} /* end of unur_tdr_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_max_intervals( struct unur_par *par, int max_ivs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ivs   ... maximum number of intervals                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return 0;
  }

  /* store date */
  PAR.max_ivs = max_ivs;

  /* changelog */
  par->set |= TDR_SET_MAX_IVS;

  return 1;

} /* end of unur_tdr_set_max_intervals() */

/*---------------------------------------------------------------------------*/


int
unur_tdr_set_center( struct unur_par *par, double center )
     /*----------------------------------------------------------------------*/
     /* set center (approximate mode) of p.d.f.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   center ... center of p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* store data */
  PAR.center = center;

  /* changelog */
  par->set |= TDR_SET_CENTER;

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_center() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_usecenter( struct unur_par *par, int usecenter )
     /*----------------------------------------------------------------------*/
     /* set flag for using center as construction point                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usecenter ... 0 = do not use,  !0 = use                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   using center as construction point is the default                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* we use a bit in variant */
  par->variant = (usecenter) ? (par->variant | TDR_VARFLAG_USECENTER) : (par->variant & (~TDR_VARFLAG_USECENTER));

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_usecenter() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_usemode( struct unur_par *par, int usemode )
     /*----------------------------------------------------------------------*/
     /* set flag for using (exact) mode as construction point                */
     /* (this overwrites "use_center"!)                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usemode   ... 0 = do not use,  !0 = use                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   using mode as construction point is the default                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* we use a bit in variant */
  par->variant = (usemode) ? (par->variant | TDR_VARFLAG_USEMODE) : (par->variant & (~TDR_VARFLAG_USEMODE));

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_usemode() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_c( struct unur_par *par, double c )
     /*----------------------------------------------------------------------*/
     /* set parameter c for transformation T_c                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   c    ... parameter c                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* check new parameter for generator */
  if (c > 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c > 0");
    return 0;
  }
  /** TODO: ... **/
/*    if (c <= -1.) { */
/*      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c <= -1 only if domain is bounded. Use `TABL' method then."); */
/*      return 0; */
/*    } */
  /** TODO: ... **/
  if (c < -0.5) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"c < -0.5 not implemented yet");
    return 0;
  }
  if (c != 0 && c > -0.5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
    c = -0.5;
  }
    
  /* store date */
  PAR.c_T = c;

  /* changelog */
  par->set |= TDR_SET_C;

  return 1;

} /* end of unur_tdr_set_c() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TDR_VARFLAG_VERIFY) : (par->variant & (~TDR_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_verify() */

/*****************************************************************************/

struct unur_gen *
_unur_tdr_init( struct unur_par *par )
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

  /* check input */
  if ( par->method != UNUR_METH_TDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tdr_create(par);
  if (!gen) { free(par); return NULL; }

  /* get starting points */
  if (!_unur_tdr_get_starting_cpoints(par,gen) ) {
    free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* compute intervals for given starting points */
  if ( !_unur_tdr_get_starting_intervals(par,gen) ) {
    free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* we have to update the maximal number of intervals,
     if the user wants more starting points. */
  if (GEN.n_ivs > GEN.max_ivs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"maximal number of intervals too small. increase.");
    GEN.max_ivs = GEN.n_ivs;
  }

  /* make initial guide table */
  _unur_tdr_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdr_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* is there any hat at all ? */
  if (GEN.Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdr_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of _unur_tdr_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_tdr_create( struct unur_par *par )
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
  unsigned variant;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* which transformation */
  if      (PAR.c_T == 0.)    variant = TDR_VAR_T_LOG;
  else if (PAR.c_T == -0.5)  variant = TDR_VAR_T_SQRT;
  else                       variant = TDR_VAR_T_POW;
  par->variant = (par->variant & (~TDR_VARMASK_T)) | variant;

  /** TODO: remove this **/
  if ((par->variant & TDR_VARMASK_T) == TDR_VAR_T_POW) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"c != 0. and c != -0.5 not implemented!");
    return NULL;
  }

  /* routines for sampling and destroying generator */
  gen->destroy = _unur_tdr_free;
  gen->reinit = _unur_reinit_error;

  if (par->variant & TDR_VARFLAG_VERIFY)
    SAMPLE = _unur_tdr_sample_check;
  else
    switch( par->variant & TDR_VARMASK_T ) {
    case TDR_VAR_T_LOG:
      SAMPLE = _unur_tdr_sample_log;
      break;
    case TDR_VAR_T_SQRT:
      SAMPLE = _unur_tdr_sample_sqrt;
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
      SAMPLE = NULL;
      return NULL;
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return NULL;
    }

  /* set all pointers to NULL */
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.iv          = NULL;
  GEN.n_ivs       = 0;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */
  GEN.c_T = PAR.c_T;                /* parameter for transformation          */

  /* bounds for adding construction points  */
  GEN.max_ivs = PAR.max_ivs;        /* maximum number of segments            */
  GEN.max_ratio = PAR.max_ratio;    /* bound for ratio  Atotal / Asqueeze    */
  GEN.bound_for_adding = PAR.bound_for_adding;

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* mode known and in given domain ?? */
  if ( !(par->distr->set & UNUR_DISTR_SET_MODE)
       || (DISTR.mode < DISTR.BD_LEFT)
       || (DISTR.mode > DISTR.BD_RIGHT))
    /* we cannot use the mode as construction point */
    par->variant = par->variant & (~TDR_VARFLAG_USEMODE);

  /* center known ?? */
  if (!(par->set & TDR_SET_CENTER))
    /* we cannot use the center as construction point */
    par->variant = par->variant & (~TDR_VARFLAG_USECENTER);
  else {
    /* center must be in domain */
    PAR.center = max(PAR.center,DISTR.BD_LEFT);
    PAR.center = min(PAR.center,DISTR.BD_RIGHT);
  }

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_tdr_create() */

/*****************************************************************************/

double
_unur_tdr_sample_log( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator; T(x) = log(x)                                 */
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
  struct unur_tdr_interval *iv, *pt;
  double u,v,x;
  double sqx, hx, fx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    u -= iv->Acum;    /* result: u in (-A_hat, 0) */

    /* left or right side of hat */
    if (-u < iv->Ahatr) { /* right */
      pt = iv->next;
      /* u unchanged */
    }
    else {                /* left */
      pt = iv;
      u += iv->Ahatl + iv->Ahatr;
    }

    /* random variate */
    if (pt->dTfx == 0.)
      x = pt->x + u/pt->fx;
    else
      x = pt->x + log(pt->dTfx / pt->fx * u + 1.) / pt->dTfx;     /** TODO: possible over/underflow **/

    /** TODO: possible over/underflow in the following calculations (?) **/

    /* accept or reject */
    hx = pt->fx * exp(pt->dTfx*(x - pt->x));    /* value of hat at x */   
    v = _unur_call_urng(gen) * hx;  /* a random point between 0 and hat at x */

    /* below mininum of density in interval ? */
    if (v <= iv->fx && v <= iv->next->fx)
      return x;

    /* below squeeze ? */
    sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(x - iv->x)) : 0.;     /* value of squeeze at x */
    if (v <= sqx)
      return x;

    /* value of p.d.f. at x */
    fx = PDF(x);

    /* being above squeeze is bad. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze)
      if ( !_unur_tdr_interval_split(gen,iv,x,fx) ) {
	/* condition for pdf is violated! */
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	/* replace sampling routine by dummy routine that just returns INFINITY */
	SAMPLE = _unur_sample_cont_error;
	return INFINITY;
      }

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject and try again */

  }
} /* end of _unur_tdr_sample_log() */

/*****************************************************************************/

double
_unur_tdr_sample_sqrt( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator; T(x) = -1./sqrt(x)                            */
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
  struct unur_tdr_interval *iv, *pt;
  double u,v,x;
  double sqx, hx,fx;
  double Tsqx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    u -= iv->Acum;    /* result: u in (-A_hat, 0) */

    /* left or right side of hat */
    if (-u < iv->Ahatr) { /* right */
      pt = iv->next;
      /* u unchanged */
    }
    else {                /* left */
      pt = iv;
      u += iv->Ahatl + iv->Ahatr;
    }

    /* random variate */
    if (pt->dTfx == 0.)
      x = pt->x + u/pt->fx;
    else {
      /* it would be less expensive to use:
	 x = pt->x + pt->Tfx/pt->dTfx * (1. - 1./(1. + pt->dTfx * pt->Tfx * u) )
	 however, this is unstable for small pt->dTfx */
      x = pt->x + (pt->Tfx*pt->Tfx*u) / (1.-pt->Tfx*pt->dTfx*u);  
      /* It cannot happen, that the denominator becomes 0 ! */
      /** TODO: underflow possible ?? **/
    }

    /** TODO: possible over/underflow in the following calculations (?) **/

    /* accept or reject */
    Thx = pt->Tfx + pt->dTfx * (x - pt->x);     /* transformed hat at x */ 
    hx = 1./(Thx*Thx);
    v = _unur_call_urng(gen) * hx;  /* a random point between 0 and hat at x */

    /* below mininum of density in interval ? */
    if (v <= iv->fx && v <= iv->next->fx)
      return x;

    /* below squeeze ? */
    Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (x - iv->x)) : -INFINITY; /* transformed squeeze at x */ 
    sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
    if (v <= sqx)
      return x;

    /* value of p.d.f. at x */
    fx = PDF(x);

    /* being above squeeze is bad. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze)
      if ( !_unur_tdr_interval_split(gen,iv,x,fx) ) {
	/* condition for pdf is violated! */
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	/* replace sampling routine by dummy routine that just returns INFINITY */
	SAMPLE = _unur_sample_cont_error;
	return INFINITY;
      }

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject and try again */

  }
} /* end of _unur_tdr_sample_sqrt() */

/*****************************************************************************/

double
_unur_tdr_sample_check( struct unur_gen *gen )
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
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... p.d.f.                                                     */
     /*   Tf  ... transformed p.d.f.                                         */
     /*   dTf ... derivative of transformed p.d.f.                           */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
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
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   squeeze(x) = 1 / (Tf(x0) + sq * (x-x0))^2                          */
     /*                                                                      */
     /*   left hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                  */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = 1 / (Tf(x1) + (Tf)'(x1) * (x-x1))^2                 */
     /*   generation:                                                        */
     /*      X = x1 + (Tf(x1)^2 * U) / (1 - Tf(x1) * (Tf)'(x1) * U)          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_tdr_interval *iv, *pt;
  double u,v,x;
  double fx, sqx, hx;
  double Tfx, Tsqx, Thx;
  int error = 0;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    u -= iv->Acum;    /* result: u in (-A_hat, 0) */

    /* left or right side of hat */
    if (-u < iv->Ahatr) { /* right */
      pt = iv->next;
      /* u unchanged */
    }
    else {                /* left */
      pt = iv;
      u += iv->Ahatl + iv->Ahatr;
    }

    /* random variate */
    if (pt->dTfx == 0.)
      x = pt->x + u/pt->fx;
    else
      switch( gen->variant & TDR_VARMASK_T ) {
      case TDR_VAR_T_LOG:
	x = pt->x + log(pt->dTfx / pt->fx * u + 1.) / pt->dTfx;     /** TODO: possible over/underflow **/
	break;
      case TDR_VAR_T_SQRT:
	x = pt->x + (pt->Tfx*pt->Tfx*u) / (1.-pt->Tfx*pt->dTfx*u);  
	/* It cannot happen, that the denominator becomes 0 (in theory!) */
	/** TODO: underflow possible ?? **/
	break;
      case TDR_VAR_T_POW:
	/** TODO **/
	x = 0.;
	break;
      default:  /* this should not happen */
	x = 0.;
      }

    /** TODO: possible over/underflow in the following calculations (?) **/
    fx = PDF(x);                                /* value of p.d.f. at x */
    Thx = pt->Tfx + pt->dTfx * (x - pt->x);     /* transformed hat at x */ 
    Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (x - iv->x)) : -INFINITY; /* transformed squeeze at x */ 

    /** TODO: possible over/underflow in the following calculations (?) **/
    switch( gen->variant & TDR_VARMASK_T ) {
    case TDR_VAR_T_LOG:
      Tfx = (fx>0.) ? log(fx) : -INFINITY;
      hx = pt->fx * exp(pt->dTfx*(x - pt->x));    /* value of hat at x */   
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(x - iv->x)) : 0.;     /* value of squeeze at x */
      break;
    case TDR_VAR_T_SQRT:
      Tfx = (fx>0.) ? -1./sqrt(fx) : -INFINITY;
      hx = 1./(Thx*Thx);
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
      Tfx = 0.;
      hx = 0.;
      sqx = 0.;
      break;
    default:  /* this should not happen */
      Tfx = 0.; hx = 0.; sqx = 0.;
    }

    /* check result */
    if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (Tfx > Thx * (1.+DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"pdf > hat. Not T-concave!");
      error = 1;
    }
    if (Tsqx > Tfx * (1.+DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"pdf < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_debug_sample( gen, iv, pt, x, fx, hx, sqx ); 
#endif

    /* accept or reject */
    v = _unur_call_urng(gen) * hx;  /* a random point between 0 and hat at x */

    /* below mininum of density in interval ? */
    if (v <= iv->fx && v <= iv->next->fx)
      return x;

    /* below squeeze ? */
    if (v <= sqx)
      return x;

    /* being above squeeze is bad. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze)
      if ( !_unur_tdr_interval_split(gen,iv,x,fx) ) {
	/* condition for pdf is violated! */
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	/* replace sampling routine by dummy routine that just returns INFINITY */
	SAMPLE = _unur_sample_cont_error;
	return INFINITY;
      }

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject and try again */

  }
} /* end of _unur_tdr_sample_check() */

/*****************************************************************************/

void
_unur_tdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */


#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdr_debug_free(gen);
#endif

  /* free linked list of intervals and others */
  //   _unur_free_mblocks(GEN.mblocks);

  /* free other memory not stored in list */
  _unur_free_genid(gen);
  free(GEN.guide);
  free(gen);

} /* end of _unur_tdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static int
_unur_tdr_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
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
     /*   1 ... if successful                                                */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int use_center, use_mode, is_mode, was_mode;
  int i, is_increasing;
  double extra_cpoint;

  /* check arguments */
  CHECK_NULL(par,0);  COOKIE_CHECK(par,CK_TDR_PAR,0);
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* use mode as construction point ? */
  use_mode = (par->variant & TDR_VARFLAG_USEMODE) ? TRUE : FALSE;

  /* use center as construction point ? */
  use_center = (!use_mode && (par->variant & TDR_VARFLAG_USECENTER)) ? TRUE : FALSE;

  /* add extra construction point        */
  /* (use either mode or center or none) */
  extra_cpoint = use_mode ? DISTR.mode : (use_center ? PAR.center : 0. );

  /* reset counter of intervals */
  GEN.n_ivs = 0;

  /* prepare for computing construction points */
  if (!PAR.starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  ( DISTR.BD_LEFT  <= -INFINITY ) ? -M_PI/2. : atan(DISTR.BD_LEFT  - PAR.center);  
    right_angle = ( DISTR.BD_RIGHT >= INFINITY )  ? M_PI/2.  : atan(DISTR.BD_RIGHT - PAR.center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR.n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = DISTR.BD_LEFT;
  if (use_mode && DISTR.mode <= x) {
    /* this is the mode of the distribution */
    is_mode = TRUE;
    use_mode = FALSE;  /* do not use the mode again */
    is_increasing = FALSE;
  }
  else if (use_center && PAR.center <= x) {
    is_mode = FALSE;
    use_center = FALSE;     /* do not use the center again */
    is_increasing = TRUE;   /* the center may be left of (unknown) mode */
  }
  else {
    is_mode = FALSE;
    is_increasing = TRUE;
  }
    
  fx = fx_last = (x <= -INFINITY) ? 0. : PDF(x);
  iv = GEN.iv = _unur_tdr_interval_new( gen, x, fx, is_mode );
  if (iv == NULL) return 0;  /* pdf(x) < 0 !! */


  /* now all the other points */
  for( i=0; i<=PAR.n_starting_cpoints; i++ ) {
    was_mode = is_mode;

    /* construction point */
    if (i < PAR.n_starting_cpoints) {
      if (PAR.starting_cpoints) {   
	/* construction points provided by user */
	x = PAR.starting_cpoints[i];
	/* check starting point */
	if (x <= DISTR.BD_LEFT || x >= DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle ) + PAR.center;
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store 
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /* insert mode or center ? */
    if ((use_mode || use_center) && x >= extra_cpoint) {
      is_mode = use_mode;              /* the next construction point is the mode */
      use_center = use_mode = FALSE;   /* we use the mode only once (of course) */
      if (x>extra_cpoint) {
	x = extra_cpoint;     /* use the mode now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!PAR.starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x == extra_cpoint --> nothing to do */
    }
    else
      is_mode = FALSE;

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of p.d.f. at starting point */
    fx = (x >= INFINITY) ? 0. : PDF(x);

    /* check value of p.d.f. at starting point */
    if (!is_increasing && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"p.d.f. not unimodal!");
      return 0;
    }
    if (is_mode && fx < fx_last * (1.-DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"mode -> ignore");
      continue;
    }
    if (was_mode && fx > fx_last * (1+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"mode");
      return 0;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such point */
      if (is_increasing) {
	/* p.d.f. is still increasing, i.e., constant 0 til now */
	if (i<PAR.n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the p.d.f. is constant 0 on all construction points.
	     then we need both boundary points. */
	  iv->x = x;  /* we only have to change x, everything else remains unchanged */
	  x_last = x;
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with pdf(x) > 0 */
	break;
    }
    
    /* need a new interval */
    iv = iv->next = _unur_tdr_interval_new( gen, x, fx, is_mode );
    if (iv == NULL) return 0;  /* pdf(x) < 0 !! */

    /* p.d.f. still increasing ? */
    if (is_increasing && fx < fx_last)
      is_increasing = 0;

    /* store last computed values */
    x_last = x;
    fx_last = fx;

  }

  /* we have left the loop with the right boundary of the support of p.d.f.
     make shure that we will never use iv for sampling. */
  iv->Asqueeze = iv->Ahatl = iv->Ahatr = iv->sq = 0.;
  iv->Acum = INFINITY;
  iv->next = NULL;         /* terminate list */
  --(GEN.n_ivs);           /* we do not count this interval */

  /* o.k. */
  return 1;

} /* end of _unur_tdr_get_starting_cpoints() */

/*****************************************************************************/

static int
_unur_tdr_get_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              /* construction point, value of p.d.f. at x */

  /* check arguments */
  CHECK_NULL(par,0);     COOKIE_CHECK(par,CK_TDR_PAR,0);
  CHECK_NULL(gen,0);     COOKIE_CHECK(gen,CK_TDR_GEN,0);
  CHECK_NULL(GEN.iv,0);  COOKIE_CHECK(GEN.iv,CK_TDR_IV,0); 

  /* compute paramters for all intervals */
  for( iv=GEN.iv; iv->next != NULL; ) {

    /* compute parameters for interval */
    switch (_unur_tdr_interval_parameter(gen, iv)) {
    case 0:     /* p.d.f. not T-concave */
      return 0;
    case 1:     /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case -1:    /* interval unbounded */
      /* split interval */
      break;
    case -2:    /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN.n_ivs);

      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make shure that we will never use this segment */
	iv->Asqueeze = iv->Ahatl = iv->Ahatr = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      continue;
    }

    /* area below hat infinite.
       insert new construction point. */
    x = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */

    /* value of p.d.f. at x */
    fx = PDF(x);

    /* add a new interval, but check if we had to used too many intervals */
    if (GEN.n_ivs >= GEN.max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return 0;
    }
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) return 0; /* pdf(x) < 0. */

    /* insert into linked list */
    iv_new->next = iv->next;
    iv->next = iv_new;
  }

  /* o.k. */
  return 1;

} /* end of _unur_tdr_get_starting_intervals() */

/*****************************************************************************/

static struct unur_tdr_interval *
_unur_tdr_interval_new( struct unur_gen *gen, double x, double fx, int is_mode )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x       ... left point of new interval                             */
     /*   fx      ... value of p.d.f. at x                                   */
     /*   is_mode ... if TRUE, x is a mode of the p.d.f.                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"pdf(x) < 0.!");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_malloc( sizeof(struct unur_tdr_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN.n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_TDR_IV);


  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->fx = fx;            /* value of p.d.f. at x */

  if (fx<=0.) {           /* --> -INFINITY */
    iv->Tfx = -INFINITY;  /* transformed density */
    iv->dTfx = INFINITY;  /* derivative of transformed density */
    return iv;
  }

  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    iv->Tfx = log(fx);
    iv->dTfx = (is_mode) ? 0. : (1./fx * dPDF(x));
    /* we can set dPDF(x) = 0. for the mode */
    break;
  case TDR_VAR_T_SQRT:
    {
      double tmp = pow(fx,1.5);
      if (tmp == 0.) {
	iv->Tfx = -INFINITY;
	iv->dTfx = INFINITY;
      }
      else {
	iv->Tfx = -1./sqrt(fx);
	/* we can set dPDF(x) = 0. for the mode */
	iv->dTfx = (is_mode) ? 0. : (0.5/tmp * dPDF(x));
      }	
    }
    break;
  case TDR_VAR_T_POW:
    /** TODO **/
    iv->Tfx = -pow(fx,GEN.c_T);
    iv->dTfx = (is_mode) ? 0. : (GEN.c_T*pow(fx,GEN.c_T-1.) * dPDF(x)); /** TODO: possible over/underflow **/
    /* we can set dPDF(x) = 0. for the mode */
  }
  
  /* the program requires dTfx > -INFINITY */
  if (iv->dTfx <= -INFINITY)
    iv->dTfx = INFINITY;

  return iv;

} /* end of _unur_tdr_interval_new() */

/*****************************************************************************/

static int
_unur_tdr_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*  -1 ... area = INFINITY                                              */
     /*  -2 ... do not add this construction point                           */
     /*   0 ... error (not p.d.f. T-concave)                                 */
     /*----------------------------------------------------------------------*/
{
  double ipt;   /* point at which the interval iv is divided into two parts */

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  CHECK_NULL(iv,0);   COOKIE_CHECK(iv,CK_TDR_IV,0); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,0);  COOKIE_CHECK(iv->next,CK_TDR_IV,0); 

  /* get division point of interval 
     (= intersection point of tangents in almost all cases) */
  if ( !_unur_tdr_interval_division_point(gen,iv,&ipt) )
    return 0;

  /* squeeze and area below squeeze */
  if (iv->Tfx > -INFINITY && iv->next->Tfx > -INFINITY) {

    /* we do not compute the slope when the construction points
       are too close. at least 8 significant digits should remain. */
    if ( (iv->x <= 0. && iv->x * (1.-FLT_EPSILON) > iv->next->x) ||
    	 (iv->x > 0.  && iv->x > iv->next->x * (1.-FLT_EPSILON)) )
      return -2;   /* construction points too close */

    /* slope of transformed squeeze */
    iv->sq = (iv->next->Tfx - iv->Tfx) / (iv->next->x - iv->x);

    /* check squeeze */
    /* iv->sq might become 0 because of round-off errors. we can ignore this.
       for simplicity we accept all squeezes with zero slope. */
    if ( iv->sq != 0. && 
	 ( (_FP_greater( iv->sq, iv->dTfx )) || 
	   (_FP_less( iv->sq, iv->next->dTfx ) && iv->next->dTfx < INFINITY) ) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. p.d.f. not T-concave!");
      return 0;
    }

    /* volume below squeeze */
    /* always integrate from point with greater value of transformed density
       to the other point */
    iv->Asqueeze = (iv->Tfx > iv->next->Tfx) ?
      _unur_tdr_interval_area( gen, iv, iv->sq, iv->next->x)
      : _unur_tdr_interval_area( gen, iv->next, iv->sq, iv->x);
  }
  else {  /* no squeeze */
    iv->sq = 0.;
    iv->Asqueeze = 0.;
  }

  /* volume below hat */
  /** TODO: it is not always a good idea to integrate from construction point
      to division point when the hat is increasing in this direction. **/
  iv->Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, ipt);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, ipt);

  /* areas below head unbounded ? */
  if (iv->Ahatl >= INFINITY || iv->Ahatr >= INFINITY)
    return -1;

  /* check area */
  if ( (iv->Asqueeze - iv->Ahatl - iv->Ahatr)/(iv->Ahatl + iv->Ahatr) > DBL_EPSILON) {
    /** TODO: is this o.k. to distiguish between roundoff errors and violated condition ?? **/
    /** TODO: possible over/underflow ( ?? ) **/
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"A(squeeze) > A(hat). p.d.f. not T-concave!");
    return 0; 
  }

  /* o.k. */
  return 1;

} /* end of _unur_tdr_interval_parameter() */

/*---------------------------------------------------------------------------*/

static int
_unur_tdr_interval_division_point( struct unur_gen *gen, struct unur_tdr_interval *iv, double *ipt )
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
     /*   1  ... if successful                                               */
     /*   0  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  double delta;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  CHECK_NULL(iv,0);   COOKIE_CHECK(iv,CK_TDR_IV,0); 

  /* case: there is no tangent at one of the boundary points of the interval */
  if (iv->dTfx >= INFINITY) { 
    *ipt = iv->x;        /* intersection point = left boundary of interval */
    return 1; 
  }
  if (iv->next->dTfx >= INFINITY) {
    *ipt = iv->next->x;   /* intersection point = right boundary of interval */
    return 1; 
  }

  /* test for T-concavity */
  if ( _FP_less( iv->dTfx, iv->next->dTfx ) ) {

    /* it might happen because of round-off errors 
       that iv->next->dTfx is almost zero although it should be large.
       thus we ignore this case. */
    if ( fabs(iv->dTfx) < DBL_EPSILON * fabs(iv->next->dTfx) ) {
      *ipt = iv->x;        /* intersection point = left boundary of interval */
      iv->dTfx = INFINITY;
      return 1; 
    }
    else if ( fabs(iv->next->dTfx) < DBL_EPSILON * fabs(iv->dTfx) ) {
      *ipt = iv->next->x;   /* intersection point = right boundary of interval */
      iv->next->dTfx = INFINITY;
      return 1; 
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). p.d.f. not T-concave!");
      return 0;
    }
  }
  /** TODO: the following test is too sensitve to roundoff errors **/
  /*    if (iv->next->Tfx > iv->x + iv->dTfx*(iv->next->x - iv->x)) { */
  /*      _unur_warning(gen->genid,UNUR_ERR_INIT,"tangent below p.d.f.  not T-concave!"); */
  /*      return 0; */
  /*    } */
  
  /* case (3) */
  /** TODO: very very important to avoid some numerical errors **/

  delta = iv->dTfx - iv->next->dTfx;

  /* case (2): computing intersection of tangents is unstable */
  if (delta < TOLERANCE) {          /** TODO: define "very small" instead of TOLERANCE **/
    /* use mean point */
    *ipt = 0.5 * (iv->x + iv->next->x);
    return 1;
  }

  /* case (1): compute intersection point of tangents (regular case) */
  *ipt = ((iv->next->Tfx - iv->Tfx) - (iv->next->dTfx * iv->next->x - iv->dTfx * iv->x)) / delta;
  /* check position of intersection point */
  if (*ipt < iv->x || *ipt > iv->next->x) {
    /** TODO: skip over simple round off error. use mean in this case */
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"intersection point of tangents not in interval. p.d.f. not T-concave!");
    return 0;
  }

  /* o.k. */
  return 1;

} /* end of _unur_tdr_interval_division_point() */

/*---------------------------------------------------------------------------*/

static double
_unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute area below piece of hat or slope in                               */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent of secant of transformed p.d.f.              */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   area                                                                    */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /* log(x)                                                                    */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(x-x0)) dx |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /* 1/sqrt(x)                                                                 */
     /*   area = | \int_{x0}^x 1/(Tf(x0) + slope*(x-x0))^2 dx |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = infinity                                      if T(f(x)) >= 0    */
     /*        = | (x-x0) / (Tf(x0)*(Tf(x0)+slope*(x-x0))) |   otherwise          */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double area = 0.;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  CHECK_NULL(iv,0);   COOKIE_CHECK(iv,CK_TDR_IV,0); 

  /* length of interval > 0 ? */
  if (x == iv->x)
    return 0.;

  /* unbounded? */
  if ( (slope >= INFINITY)         ||
       (x<=-INFINITY && slope<=0.) ||
       (x>= INFINITY && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:
    /* T(x) = log(x) */
    if (slope != 0.) {                         
      /** TODO: cannot use this case if slope is `very small' **/
      if (x<=-INFINITY || x>= INFINITY)
	area = iv->fx / slope;                                     /** TODO: possible over/underflow **/
      else
	area = iv->fx / slope * ( exp(slope*(x - iv->x)) - 1. );   /** TODO: possible over/underflow **/
    }
    else { /* hat/squeeze almost constant */
      if (x<=-INFINITY || x>= INFINITY)
	return INFINITY;
      area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_SQRT:
    /* T(x) = -1./sqrt(x) */
    if (slope != 0.) {
      /** TODO: cannot use this case if slope is `very small' **/
      double z;
      if (x<=-INFINITY || x>= INFINITY)
	area = 1. / ( iv->Tfx * slope );
      else {
	z = iv->Tfx + slope*(x - iv->x);
	if (z>=0.)      /* the hat must not cut the x-axis */
	  /** TODO: this is very sensitive to roundoff errors, when Tfx und slope are 
	      large numbers. see get_division_point, case (3) **/
	  return INFINITY; 
	else
	  area = (x - iv->x) / ( iv->Tfx * z );
      }
    }
    else { /* hat/squeeze almost constant */
      if (x<=-INFINITY || x>= INFINITY)
	return INFINITY;
      area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_POW:
    /* T(x) = -1./x^c */
    /** TODO **/
    break;
  }

  return ( (area<0.) ? -area : area );

} /* end of _unur_tdr_interval_area() */

/*****************************************************************************/

static int
_unur_tdr_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv_oldl, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* split interval iv_oldl into two intervals at point x                 */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv_oldl ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   fx      ... value of p.d.f. at x                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1  ... if successful                                               */
     /*   0  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv_newr;
  int success;

  /* check arguments */
  CHECK_NULL(gen,0);      COOKIE_CHECK(gen,CK_TDR_GEN,0);
  CHECK_NULL(iv_oldl,0);  COOKIE_CHECK(iv_oldl,CK_TDR_IV,0);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_debug_split_start( gen,iv_oldl,x,fx );
#endif

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN.n_ivs * (iv_oldl->Ahatl + iv_oldl->Ahatr - iv_oldl->Asqueeze) / (GEN.Atotal - GEN.Asqueeze))
       < GEN.bound_for_adding)
    return 1;

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"pdf(x) < 0.!");
    return 0;
  }

  /* check if the new interval is completely outside the support of p.d.f. */
  if (fx <= 0.) {

    /* one of the two boundary points must be 0, too! */
    if (iv_oldl->fx <= 0.) {
      /* chop off left part (it's out of support) */
      iv_oldl->x = x;
    }
    else if (iv_oldl->next->fx <= 0.) {
      /* chop off right part (it's out of support) */
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"p.d.f. not T-concave");
      return 0;
    }

    /* compute parameters for chopped interval */
    if ( !_unur_tdr_interval_parameter(gen,iv_oldl) ) {
      /** TODO **/
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"p.d.f. not T-concave");
      return 0;
    }

    /* for _unur_tdr_debug_split_stop only */
    iv_newr = iv_oldl;

  }
  else {

    /* we need a new interval */
    iv_newr = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_newr == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 0;
    }

    /* link into list */
    iv_newr->next = iv_oldl->next;
    iv_oldl->next = iv_newr;

    /* compute parameters for interval */
    if ( (success = _unur_tdr_interval_parameter(gen, iv_oldl)) <= 0 ||
	 (success = _unur_tdr_interval_parameter(gen, iv_newr)) <= 0 ) {

      /* p.d.f. not T-concave, or new interval to narrow, 
	 or area below hat not bounded */

      /* new construction point not suitable --> do not add */
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
#ifdef UNUR_ENABLE_LOGGING
      /* write info into log file */
      if (gen->debug & TDR_DEBUG_SPLIT) 
	_unur_tdr_debug_split_stop( gen,iv_oldl,iv_newr );
#endif

      /* remove from linked list */
      iv_oldl->next = iv_newr->next;
      --(GEN.n_ivs);   /* decrement counter for intervals */
      free( iv_newr );

      /* we have to restore the old interval.
	 (this case should not happen, so it is faster not to make a 
	 backup of the old interval) */
      if ( !_unur_tdr_interval_parameter(gen, iv_oldl) ) {
	/* this should not happen:
	   Cannot restore interval. */
	_unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
	exit (EXIT_FAILURE);
      }

      if ( !success ) {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"p.d.f. not T-concave");
	return 0;
      }
      else
	return 1;
    }
  }

  /* update guide table */ 
  _unur_tdr_make_guide_table(gen);
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_debug_split_stop( gen,iv_oldl,iv_newr );
#endif

  /* o.k. */
  return 1;

} /* end of _unur_tdr_interval_split() */

/*****************************************************************************/

static int
_unur_tdr_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 (--> successful)                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN.guide) {
    int max_guide_size = (GEN.guide_factor > 0.) ? (GEN.max_ivs * GEN.guide_factor) : 1;
    GEN.guide = _unur_malloc( max_guide_size * sizeof(struct unur_tdr_interval*) );
  }

  /* first we need cumulated areas in segments */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,0);
    Acum += iv->Ahatl + iv->Ahatr;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN.Atotal = Acum;
  GEN.Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN.guide_size = (int)(GEN.n_ivs * GEN.guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN.Atotal / GEN.guide_size;
  Acum=0.;
  for( j=0, iv=GEN.iv; j < GEN.guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDR_IV,0);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   /* this is the last virtual intervall --> do not use */
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN.guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide[j] = iv;

  return 1;
} /* end of _unur_tdr_make_guide_table() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_tdr_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_TDR_PAR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = transformed density rejection\n",gen->genid);
  fprintf(log,"%s: transformation T_c(x) = ",gen->genid);
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    fprintf(log,"log(x)  ... c = 0");                   break;
  case TDR_VAR_T_SQRT:
    fprintf(log,"-1/sqrt(x)  ... c = -1/2");            break;
  case TDR_VAR_T_POW:
    fprintf(log,"-x^(%g)  ... c = %g",PAR.c_T,PAR.c_T); break;
  }
  _unur_print_if_default(par,TDR_SET_C);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tdr_sample_",gen->genid);
  if (par->variant & TDR_VARFLAG_VERIFY)
    fprintf(log,"check()\n");
  else
    switch( gen->variant & TDR_VARMASK_T ) {
    case TDR_VAR_T_LOG:   
      fprintf(log,"log()\n");   break;
    case TDR_VAR_T_SQRT:
      fprintf(log,"sqrt()\n");  break;
    case TDR_VAR_T_POW:
      fprintf(log,"pow()\n");   break;
    }
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: center = %g",gen->genid,PAR.center);
  _unur_print_if_default(par,TDR_SET_CENTER);
  if (par->variant & TDR_VARFLAG_USEMODE)
    fprintf(log,"\n%s: use mode as construction point",gen->genid);
  else if (par->variant & TDR_VARFLAG_USECENTER)
    fprintf(log,"\n%s: use center as construction point",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR.max_ivs);
  _unur_print_if_default(par,TDR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,TDR_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
  _unur_print_if_default(par,TDR_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,TDR_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & TDR_SET_STP)
    for (i=0; i<PAR.n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",PAR.starting_cpoints[i]);
    }
  else
    fprintf(log," use \"equdistribution\" rule [default]");
  fprintf(log,"\n%s:\n",gen->genid);
  
  _unur_tdr_debug_intervals(gen);

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_init() */

/*****************************************************************************/

static void
_unur_tdr_debug_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_tdr_debug_intervals(gen);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_free() */

/*****************************************************************************/

static void
_unur_tdr_debug_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tdr_interval *iv;
  double sAsqueeze, sAhatl, sAhatr;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:Intervals: %d\n",gen->genid,GEN.n_ivs);
  if (GEN.iv) {
    if (gen->debug & TDR_DEBUG_IV) {
      fprintf(log,"%s: Nr.            tp          f(tp)        T(f(tp))    d(T(f(tp)))        squeeze\n",gen->genid);
      for (iv = GEN.iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
	fprintf(log,"%s:[%3d]: %#12.6g   %12.6g   %#12.6g   %#12.6g   %12.6g\n", gen->genid, i,
		iv->x, iv->fx, iv->Tfx, iv->dTfx, iv->sq);
      }
      COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
      fprintf(log,"%s:[...]: %#12.6g   %12.6g   %#12.6g   %#12.6g\n", gen->genid,
	      iv->x, iv->fx, iv->Tfx, iv->dTfx);
    }
    fprintf(log,"%s:\n",gen->genid);
  }
  else
    fprintf(log,"%s: No intervals !\n",gen->genid);

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TDR_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\tbelow squeeze\t\t  below hat (left and right)\t\t  cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
    if (GEN.iv) {
      for (iv = GEN.iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
	sAsqueeze += iv->Asqueeze;
	sAhatl += iv->Ahatl;
	sAhatr += iv->Ahatr;
	fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		iv->Asqueeze, iv->Asqueeze * 100. / GEN.Atotal,
		iv->Ahatl, iv->Ahatr, (iv->Ahatl+iv->Ahatr) * 100. / GEN.Atotal, 
		iv->Acum, iv->Acum * 100. / GEN.Atotal);
      }
      fprintf(log,"%s:       ----------  ---------  |  ------------------------  ---------  +\n",gen->genid);
      fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)            %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAsqueeze, sAsqueeze * 100. / GEN.Atotal,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / GEN.Atotal);
      fprintf(log,"%s:\n",gen->genid);
    }
  }

  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_tdr_debug_intervals() */

/*****************************************************************************/

static void
_unur_tdr_debug_sample( struct unur_gen *gen, 
			struct unur_tdr_interval *iv, 
			struct unur_tdr_interval *pt, 
			double x, double fx, double hx, double sqx )
     /*----------------------------------------------------------------------*/
     /* write info about generated point                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   pt  ... pointer to interval that stores construction point         */
     /*   x   ... generated point                                            */
     /*   fx  ... value of p.d.f. at x                                       */
     /*   hx  ... value of hat at x                                          */
     /*   sqx ... value of squeeze at x                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv,/*void*/);   COOKIE_CHECK(iv,CK_TDR_IV,/*void*/);
  CHECK_NULL(pt,/*void*/);   COOKIE_CHECK(pt,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (iv == pt)
    fprintf(log,"%s: point generated in left part:\n",gen->genid);
  else
    fprintf(log,"%s: point generated in right part:\n",gen->genid);

  fprintf(log,"%s: construction point: x0 = %g\n",gen->genid,pt->x);
  fprintf(log,"%s: transformed hat     Th(x) = %g + %g * (x - %g)\n",gen->genid,pt->Tfx,pt->dTfx,pt->x);
  fprintf(log,"%s: transformed squeeze Ts(x) = %g + %g * (x - %g)\n",gen->genid,iv->Tfx,iv->sq,iv->x);
  fprintf(log,"%s: generated point: x = %g\n",gen->genid,x);
  fprintf(log,"%s:  h(x) = %.20g\n",gen->genid,hx);
  fprintf(log,"%s:  f(x) = %.20g\n",gen->genid,fx);
  fprintf(log,"%s:  s(x) = %.20g\n",gen->genid,sqx);
  fprintf(log,"%s:    hat: x - x0 = %g",gen->genid,x-pt->x);
  if (x < pt->x && iv == pt) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    h(x) - f(x) = %g",gen->genid,hx-fx);
  if (hx<fx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    squeeze: x - x0 = %g",gen->genid,x-iv->x);
  if (x > pt->x && iv != pt) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    f(x) - s(x) = %g",gen->genid,fx-sqx);
  if (fx<sqx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_sample() */

/*****************************************************************************/

static void
_unur_tdr_debug_split_start( struct unur_gen *gen, struct unur_tdr_interval *iv, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting interval                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   x   ... split at this point                                        */
     /*   fx  ... value of p.d.f. at x                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv,/*void*/);   COOKIE_CHECK(iv,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: split interval at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(log,"%s: old interval:\n",gen->genid);
  fprintf(log,"%s:   left  construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->x,iv->fx);
  fprintf(log,"%s:   right construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->fx);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv->Asqueeze,iv->Asqueeze*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv->Ahatl+iv->Ahatr-iv->Asqueeze),(iv->Ahatl+iv->Ahatr-iv->Asqueeze)*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv->Ahatl,iv->Ahatr,(iv->Ahatl+iv->Ahatr)*100./GEN.Atotal);

  fflush(log);

} /* end of _unur_tdr_debug_split_start() */

/*****************************************************************************/

static void
_unur_tdr_debug_split_stop( struct unur_gen *gen, 
			    struct unur_tdr_interval *iv_left, 
			    struct unur_tdr_interval *iv_right )
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
  CHECK_NULL(gen,/*void*/);       COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv_left,/*void*/);   COOKIE_CHECK(iv_left,CK_TDR_IV,/*void*/);
  CHECK_NULL(iv_right,/*void*/);  COOKIE_CHECK(iv_right,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: inserted point:\n",gen->genid);
  fprintf(log,"%s: x = %g, f(x) = %g, Tf(x)=%g, dTf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->fx, iv_right->Tfx, iv_right->dTfx, iv_right->sq);
  fprintf(log,"%s: new intervals:\n",gen->genid);
  fprintf(log,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(log,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(log,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);

  fprintf(log,"%s: left interval:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv_left->Asqueeze,
	  iv_left->Asqueeze*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv_left->Ahatl + iv_left->Ahatr - iv_left->Asqueeze),
	  (iv_left->Ahatl + iv_left->Ahatr - iv_left->Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv_left->Ahatl,
	  iv_left->Ahatr,
	  (iv_left->Ahatl + iv_left->Ahatr) * 100./GEN.Atotal);

  if (iv_left == iv_right)
    fprintf(log,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(log,"%s: right interval:\n",gen->genid);
    fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,
	    iv_right->Asqueeze*100./GEN.Atotal);
    fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    (iv_right->Ahatl + iv_right->Ahatr - iv_right->Asqueeze),
	    (iv_right->Ahatl + iv_right->Ahatr - iv_right->Asqueeze) * 100./GEN.Atotal);
    fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	    iv_right->Ahatl,
	    iv_right->Ahatr,
	    (iv_right->Ahatl + iv_right->Ahatr) * 100./GEN.Atotal);
  }

  fprintf(log,"%s: total areas:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_split_stop() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
