/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr.h                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and .... of a T-concave distribution                     *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:  pointer to the density, ....                                 *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *      double *pdf_param    ... parameters of p.d.f.                        *
 *                               (default: NULL)                             *
 *      int     n_pdf_param  ... number of parameters of p.d.f.              *
 *                               (default: 0)                                *
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
 *     `unur_set_cpoints()´ with the NULL pointer as the second argument     *
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
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Thu Oct 14 15:15:26 CEST 1999                        *
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

#define GENTYPE "TDR"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
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
/* return 0 if p.d.f. is not T-concave.                                      */
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

static struct unur_tdr_interval *_unur_tdr_iv_stack_pop( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* pop an interval from the stack of free intervals.                         */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_iv_stack_push( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* push the last popped interval back onto the stack.                        */
/*---------------------------------------------------------------------------*/


#if UNUR_DEBUG & UNUR_DB_INFO
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

#define PAR     par->data.tdr
#define GEN     gen->data.tdr
#define SAMPLE  gen->sample.cont

#define PDF(x) ((*(GEN.pdf))((x),GEN.pdf_param,GEN.n_pdf_param))
#define dPDF(x) ((*(GEN.dpdf))((x),GEN.pdf_param,GEN.n_pdf_param))

/*---------------------------------------------------------------------------*/

/* indicate variant of method */
#define TDR_MASK_T      0x03UL          /* indicates transformation */

#define TDR_METH_SQRT   0x01UL
#define TDR_METH_LOG    0x02UL
#define TDR_METH_POW    0x03UL

/* Special debugging flags (do not use the first 3 bits) */
#define TDR_DB_SPLIT    0x010UL
#define TDR_DB_SAMPLE   0x020UL
#define TDR_DB_IV       0x040UL

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_tdr_new( double (*pdf)(double x,double *pdf_param, int n_pdf_params), 
	      double (*dpdf)(double x,double *pdf_param, int n_pdf_params) )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdf  ... probability density function of the desired distribution  */
     /*   dpdf ... derivative of p.d.f.                                      */
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
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_TDR_PAR);

  /* copy input */
  PAR.pdf                 = pdf;    /* pointer to p.d.f.                     */
  PAR.dpdf                = dpdf;   /* pointer to derivative of p.d.f.       */

  /* set default values */
  PAR.pdf_param           = NULL;   /* no parameters for pdf                 */
  PAR.n_pdf_param         = 0;      /* number of parameters                  */
  PAR.bleft               = - INFINITY; /* left boundary of domain           */
  PAR.bright              = INFINITY;   /* right boundary of domain          */
  PAR.mode                = 0.;     /* (exact!) location of mode             */

  PAR.guide_factor        = 3.;     /* size of guide table / number of intervals */

  PAR.c_T                 = -0.5;   /* parameter for transformation (-1. <= c < 0.) */

  PAR.starting_cpoints    = NULL;   /* pointer to array of starting points   */
  PAR.n_starting_cpoints  = 10;     /* number of starting points             */
  PAR.max_ivs             = 50;     /* maximum number of intervals           */
  PAR.max_ratio           = 0.95;   /* bound for ratio  Atotal / Asqueeze    */
  PAR.bound_for_adding    = 0.5;    /* do not add a new construction point in an interval,
				       where ambigous region is too small, i.e. if 
				       area / ((A_hat - A_squeeze)/number of segments) < bound_for_adding */
 
  par->method             = UNUR_METH_TDR;  /* method and default variant    */
  par->set                = 0UL;    /* inidicate default parameters          */    
  par->urng               = unur_get_default_urng(); /* use default urng     */

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  _unur_set_genid(par,GENTYPE);     /* set generator identifier              */

  /* routine for starting generator */
  par->init = unur_tdr_init;

  return par;

} /* end of unur_tdr_new() */

/*****************************************************************************/

struct unur_gen *
unur_tdr_init( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tdr_create(par);
  if (!gen) { free(par); return NULL; }

  /* get starting points */
  if (!_unur_tdr_get_starting_cpoints(par,gen) ) {
    free(par); unur_tdr_free(gen);
    return NULL;
  }

  /* compute intervals for given starting points */
  if ( !_unur_tdr_get_starting_intervals(par,gen) ) {
    free(par); unur_tdr_free(gen);
    return NULL;
  }

  /* we have to update the maximal number of intervals,
     if the user wants more starting points. */
  if (GEN.n_ivs > GEN.max_ivs) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"maximal number of intervals too small. increase.");
    GEN.max_ivs = GEN.n_ivs;
  }

  /* make initial guide table */
  _unur_tdr_make_guide_table(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_tdr_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* is there any hat at all ? */
  if (GEN.Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_INIT,"cannot construct hat function. bad construction points.");
    unur_tdr_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of unur_tdr_init() */

/*****************************************************************************/

double
unur_tdr_sample_log( struct unur_gen *gen )
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
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    COOKIE_CHECK(iv,CK_TDR_IV,0.);

    u *= GEN.Atotal;
    while (iv->Acum < u) {
      iv = iv->next;
      COOKIE_CHECK(iv,CK_TDR_IV,0.);
    }

    /* reuse of uniform random number */
    u = iv->Acum - u;

    /* left or right side of hat */
    if (u < iv->Ahatl) {   /* left */
      pt = iv;
      /* u unchanged */
    }
    else {                 /* right */
      pt = iv->next;
      u = u - iv->Ahatl - iv->Ahatr;
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
      _unur_tdr_interval_split(gen,iv,x,fx);

    /** TODO: test fx >= sqx ?? (use transformed denisty) **/

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject */

  }
} /* end of unur_tdr_sample_log() */

/*****************************************************************************/

double
unur_tdr_sample_sqrt( struct unur_gen *gen )
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
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    COOKIE_CHECK(iv,CK_TDR_IV,0.);

    u *= GEN.Atotal;
    while (iv->Acum < u) {
      iv = iv->next;
      COOKIE_CHECK(iv,CK_TDR_IV,0.);
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
      _unur_tdr_interval_split(gen,iv,x,fx);

    /** TODO: test fx >= sqx ?? (use transformed denisty) **/

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject */

  }
} /* end of unur_tdr_sample_sqrt() */

/*****************************************************************************/

double
unur_tdr_sample_check( struct unur_gen *gen )
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
     /*   squeeze(x) = f(x0) * \exp(sq * (x-x0))                             */
     /*                                                                      */
     /*   left hat:                                                          */
     /*   X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )             */
     /*   U ~ U(0,area below left hat)                                       */
     /*                                                                      */
     /*   right hat:                                                         */
     /*   X = x1 + 1/(Tf)'(x1) * \log( (Tf)'(x1)/f(x1) * U + 1 )             */
     /*   U ~ U(- area below right hat,0)                                    */
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
     /*   right hat(x) = 1 / (Tf(x0) + (Tf)'(x1) * (x-x1))^2                 */
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
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

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
      switch( gen->method & TDR_MASK_T ) {
      case TDR_METH_LOG:
	x = pt->x + log(pt->dTfx / pt->fx * u + 1.) / pt->dTfx;     /** TODO: possible over/underflow **/
	break;
      case TDR_METH_SQRT:
	x = pt->x + (pt->Tfx*pt->Tfx*u) / (1.-pt->Tfx*pt->dTfx*u);  
	/* It cannot happen, that the denominator becomes 0 (in theory!) */
	/** TODO: underflow possible ?? **/
	break;
      case TDR_METH_POW:
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
    switch( gen->method & TDR_MASK_T ) {
    case TDR_METH_LOG:
      Tfx = (fx>0.) ? log(fx) : -INFINITY;
      hx = pt->fx * exp(pt->dTfx*(x - pt->x));    /* value of hat at x */   
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(x - iv->x)) : 0.;     /* value of squeeze at x */
      break;
    case TDR_METH_SQRT:
      Tfx = (fx>0.) ? -1./sqrt(fx) : -INFINITY;
      hx = 1./(Thx*Thx);
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      break;
    case TDR_METH_POW:
      /** TODO **/
      Tfx = 0.;
      hx = 0.;
      sqx = 0.;
      break;
    default:  /* this should not happen */
      Tfx = 0.; hx = 0.; sqx = 0.;
    }

    /* check result */
    if (x<GEN.bleft || x>GEN.bright) {
      _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"generated point out of domain");
      error = 1;
    }
    if (Tfx > Thx && 
	fabs((Tfx-Thx)/Tfx) > 2. * DBL_EPSILON ) {   /* this construct should skip over simple roundoff errors */
      /** TODO: is factor 2 a good choice ?? **/
      _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf > hat. Not T-concave!");
      error = 1;
    }
    if (Tsqx > Tfx &&
	fabs((Tsqx-Tfx)/Tsqx) > 2. * DBL_EPSILON ) {   /* this construct should skip over simple roundoff errors */
      _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf < squeeze. Not T-concave!");
      error = 1;
    }

#if UNUR_DEBUG & UNUR_DB_INFO
    /* write info into log file (in case error) */
    if (error && (gen->debug & TDR_DB_SAMPLE)) 
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
      _unur_tdr_interval_split(gen,iv,x,fx);

    if (v <= fx)
      /* between p.d.f. and squeeze */
      return x;

    /* else reject */

  }
} /* end of unur_tdr_sample_check() */

/*****************************************************************************/

void
unur_tdr_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  if (gen->debug) _unur_tdr_debug_free(gen);
#endif

  /* free linked list of intervals and others */
  _unur_free_mblocks(GEN.mblocks);

  /* free other memory not stored in list */
  _unur_free_genid(gen);
  free(GEN.guide);
  free(gen);

} /* end of unur_tdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

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
  unsigned long variant;
  int i;

  /* check arguments */
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDR_GEN);

  /* which transformation */
  if      (PAR.c_T == 0.)    variant = TDR_METH_LOG;
  else if (PAR.c_T == -0.5)  variant = TDR_METH_SQRT;
  else                       variant = TDR_METH_POW;
  par->method = (par->method & (~TDR_MASK_T)) | variant;

  /** TODO: remove this **/
  if ((par->method & TDR_MASK_T) == TDR_METH_POW) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"c != 0. and c != -0.5 not implemented!");
    return NULL;
  }

  /* routines for sampling and destroying generator */
  gen->destroy = unur_tdr_free;
  if (par->method & UNUR_MASK_SCHECK)
    SAMPLE = unur_tdr_sample_check;
  else
    switch( par->method & TDR_MASK_T ) {
    case TDR_METH_LOG:
      SAMPLE = unur_tdr_sample_log;
      break;
    case TDR_METH_SQRT:
      SAMPLE = unur_tdr_sample_sqrt;
      break;
    case TDR_METH_POW:
      /** TODO **/
      SAMPLE = NULL;
      break;
    default:
      _unur_warning(gen->genid,UNUR_ERR_INIT,"internal error.");
      return NULL;
    }

  /* set all pointers to NULL */
  GEN.pdf_param   = NULL;
  GEN.n_pdf_param = 0;
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.iv          = NULL;
  GEN.n_ivs       = 0;
  GEN.iv_stack    = NULL;
  GEN.iv_free     = 0;
  GEN.mblocks     = NULL;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN.pdf = PAR.pdf;                /* p.d.f. of distribution                */
  GEN.dpdf = PAR.dpdf;              /* derivative of p.d.f.                  */
  GEN.bleft = PAR.bleft;            /* left boundary of domain               */
  GEN.bright = PAR.bright;          /* right boundary of domain              */

  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */

  GEN.c_T = PAR.c_T;                /* parameter for transformation          */

  /* bounds for adding construction points  */
  GEN.max_ivs = PAR.max_ivs;        /* maximum number of segments            */
  GEN.max_ratio = PAR.max_ratio;    /* bound for ratio  Atotal / Asqueeze    */
  GEN.bound_for_adding = PAR.bound_for_adding;

  gen->method = par->method;        /* indicates method and variant          */
  _unur_copy_urng_pointer(par,gen); /* pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);    /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);        /* copy generator identifier             */

  /* copy parameters of distribution */
  GEN.n_pdf_param = PAR.n_pdf_param;
  if( PAR.n_pdf_param > 0 ) {
    GEN.pdf_param = _unur_malloc( PAR.n_pdf_param * sizeof(double) );
    _unur_add_mblocks( &(GEN.mblocks), GEN.pdf_param);
    for (i=0; i<PAR.n_pdf_param; i++)
      GEN.pdf_param[i] = PAR.pdf_param[i];
  }
    
  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_tdr_create() */

/*****************************************************************************/

static int
_unur_tdr_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting intervals.                  */
     /* if not provided as arguments compute these                           */
     /* by means of the "equiangular rule" from AROU.                        */
     /** TODO (?) **/
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*   use_boundary ... indicates whether boundary  points are used as    */
     /*                    construction points.                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int i, use_mode, is_mode, was_mode, is_increasing;
  double mode_shift = 0.;

  /* check arguments */
  COOKIE_CHECK(par,CK_TDR_PAR,0);
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* use mode as construction point ? */
  use_mode = (par->method & UNUR_MASK_MODE) ? TRUE : FALSE;
  is_mode = was_mode = FALSE;    /* initialize boolean */

  /* check mode */
  if (use_mode &&
      ( PAR.mode < PAR.bleft || PAR.mode > PAR.bright ) ) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"mode out of domain.");
    use_mode = 0;
  }

  /* reset counter of intervals */
  GEN.n_ivs = 0;

  /* prepare for computing construction points */
  if (!PAR.starting_cpoints) {
    /* move mode into  x = 0 ? */
    mode_shift = use_mode ? PAR.mode : 0.;
    /* angles of boundary of domain */
    left_angle =  ( PAR.bleft <= -INFINITY ) ? -M_PI/2. : atan(PAR.bleft - mode_shift);  
    right_angle = ( PAR.bright >= INFINITY )  ? M_PI/2.  : atan(PAR.bright - mode_shift);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR.n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = PAR.bleft;
  fx = fx_last = (x <= -INFINITY) ? 0. : PDF(x);
  iv = GEN.iv = _unur_tdr_interval_new( gen, x, fx, FALSE );
  CHECK_NULL(iv,0);        /* case of error */
  is_increasing = 1;       /* assume pdf(x) is increasing for the first construction points */

  /* now all the other points */
  for( i=0; i<=PAR.n_starting_cpoints; i++ ) {
    was_mode = is_mode;

    /* starting point */
    if (i < PAR.n_starting_cpoints) {
      if (PAR.starting_cpoints) {   
	/* construction points provided by user */
	x = PAR.starting_cpoints[i];
	/* check starting point */
	if (x <= PAR.bleft || x >= PAR.bright) {
	  _unur_warning(gen->genid,UNUR_ERR_INIT,"starting point out of domain!");
	  continue;
	}
	if (x<=x_last) {
	  _unur_warning(gen->genid,UNUR_ERR_INIT,"starting points are not strictly monotonically increasing! skip!");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;                /** TODO: angle >= M_PI/2. !! **/
	x = tan( angle ) + mode_shift;      /** TODO: possible over/underflow **/
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store 
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = PAR.bright;
    }

    /* insert mode ? */
    if (use_mode && x >= PAR.mode) {
      use_mode = FALSE;   /* we use the mode only once (of course) */
      is_mode = TRUE;     /* the next construction point is the mode */
      if (x>PAR.mode) {
	x = PAR.mode;     /* use the mode now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!PAR.starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x==PAR.mode --> nothing to do */
    }
    else
      is_mode = FALSE;

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of p.d.f. at starting point */
    fx = (x >= INFINITY) ? 0. : PDF(x);

    /* check value of p.d.f. at starting point */
    if (!is_increasing && fx > fx_last) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"p.d.f. not unimodal!");
      return 0;
    }
    if (is_mode && (fx < fx_last)) {
      _unur_warning(gen->genid,UNUR_ERR_INIT,"wrong mode. ignore mode.");
      continue;
    }
    if (was_mode && (fx > fx_last)) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"wrong mode.");
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
    CHECK_NULL(iv,0);     /* case of error */

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
  struct unur_tdr_interval *iv, *iv_new; 
  double x,fx;              /* construction point, value of p.d.f. at x */

  /* check arguments */
  COOKIE_CHECK(par,CK_TDR_PAR,0);
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

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
      iv->next = iv->next->next;
      --GEN.n_ivs;   /* update counter. we do not free the corresponding 
			memory block; too much book keeping. */
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
      _unur_error(gen->genid,UNUR_ERR_INIT,"cannot create bounded hat!");
      return 0;
    }
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    CHECK_NULL(iv_new,0);     /* case of error */
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
  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"pdf(x) < 0.!");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_tdr_iv_stack_pop(gen);
  COOKIE_CHECK(iv,CK_TDR_IV,NULL); 

  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->fx = fx;            /* value of p.d.f. at x */

  if (fx<=0.) {           /* --> -INFINITY */
    /** TODO: fx < very small number ( ?? ) **/
    iv->Tfx = -INFINITY;  /* transformed density */
    iv->dTfx = INFINITY;  /* derivative of transformed density */
    return iv;
  }

  switch( gen->method & TDR_MASK_T ) {
  case TDR_METH_LOG:
    iv->Tfx = log(fx);                                     /** TODO: possible over/underflow **/
    iv->dTfx = (is_mode) ? 0. : (1./fx * dPDF(x));         /** TODO: possible over/underflow **/
    /* we can set dPDF(x) = 0. for the mode */
    break;
  case TDR_METH_SQRT:
    iv->Tfx = -1./sqrt(fx);                                /** TODO: possible over/underflow **/
    iv->dTfx = (is_mode) ? 0. : (0.5/pow(fx,1.5) * dPDF(x));  /** TODO: possible over/underflow **/
    if (isnan(iv->dTfx)) iv->dTfx = INFINITY;              /** TODO: isnan is not ANSI C **/
    break;
  case TDR_METH_POW:
    /** TODO **/
    iv->Tfx = -pow(fx,GEN.c_T);                           /** TODO: possible over/underflow **/
    iv->dTfx = (is_mode) ? 0. : (GEN.c_T*pow(fx,GEN.c_T-1.) * dPDF(x)); /** TODO: possible over/underflow **/
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
     /*  -2 ... construction points too close                                */
     /*   0 ... error (not p.d.f. T-concave)                                 */
     /*----------------------------------------------------------------------*/
{
  double ipt;   /* point at which the interval iv is divided into two parts */

  /* check arguments */
  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  COOKIE_CHECK(iv,CK_TDR_IV,0); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,0);
  COOKIE_CHECK(iv->next,CK_TDR_IV,0); 

  /* get division point of interval 
     (= intersection point of tangents in almost all cases) */
  if ( !_unur_tdr_interval_division_point(gen,iv,&ipt) )
    return 0;

  /* squeeze and area below squeeze */
  if (iv->Tfx > -INFINITY && iv->next->Tfx > -INFINITY) {
    /* slope of transformed squeeze */
    if (iv->x == iv->next->x)    /** TODO **/
      return -2;   /* construction points too close */
    iv->sq = (iv->next->Tfx - iv->Tfx) / (iv->next->x - iv->x);   /** TODO: possible over/underflow **/
    /* check squeeze */
    if ( (iv->sq > iv->dTfx      && iv->dTfx       > -INFINITY) ||
	 (iv->sq < iv->next->dTfx && iv->next->dTfx  < INFINITY) ) {
      _unur_warning(gen->genid,UNUR_ERR_INIT,"Squeeze too steep/flat. p.d.f. not T-concave!");
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
    _unur_warning(gen->genid,UNUR_ERR_INIT,"A(squeeze) > A(hat). p.d.f. not T-concave!");
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
     /* (3) use left or right boundary if (left-right) is nearly the same as */
     /*     left or right, respectively.                                     */
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
  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  COOKIE_CHECK(iv,CK_TDR_IV,0); 

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
  if (iv->dTfx < iv->next->dTfx) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"dTfx0 < dTfx1 (x0>x1). p.d.f. not T-concave!");
    return 0;
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
    _unur_warning(gen->genid,UNUR_ERR_INIT,"intersection point of tangents not in interval. p.d.f. not T-concave!");
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
  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  COOKIE_CHECK(iv,CK_TDR_IV,0); 

  /* length of interval > 0 ? */
  if (x == iv->x)
    return 0.;

  /* unbounded? */
  if ( (slope >= INFINITY)         ||
       (x<=-INFINITY && slope<=0.) ||
       (x>= INFINITY && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  switch( gen->method & TDR_MASK_T ) {

  case TDR_METH_LOG:
    /* T(x) = log(x) */
    if (slope != 0.) {                         
      /** TODO: cannot use this case if slope is `very small' **/
      if (x<=-INFINITY || x>= INFINITY)
	area = iv->fx / slope;                                     /** TODO: possible over/underflow **/
      else
	area = iv->fx / slope * ( exp(slope*(x - iv->x)) - 1. );   /** TODO: possible over/underflow **/
    }
    else /* hat/squeeze almost constant */  /** TODO !! **/
      area = iv->fx * (x - iv->x);
    break;

  case TDR_METH_SQRT:
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
    else /* hat/squeeze almost constant */  /** TODO !! **/
      area = iv->fx * (x - iv->x);
    break;

  case TDR_METH_POW:
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

  /* check arguments */
  COOKIE_CHECK(gen,CK_TDR_GEN,0);
  COOKIE_CHECK(iv_oldl,CK_TDR_IV,0);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug & TDR_DB_SPLIT) 
    _unur_tdr_debug_split_start( gen,iv_oldl,x,fx );
#endif

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN.n_ivs * (iv_oldl->Ahatl + iv_oldl->Ahatr - iv_oldl->Asqueeze) / (GEN.Atotal - GEN.Asqueeze))
       < GEN.bound_for_adding)
    return 0;

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
      _unur_error(gen->genid,UNUR_ERR_INIT,"p.d.f. not T-concave");
      return 0;
    }

    /* compute parameters for chopped interval */
    _unur_tdr_interval_parameter(gen,iv_oldl);
    /** TODO: check for errors! **/

    /* for _unur_tdr_debug_split_stop only */
    iv_newr = iv_oldl;

  }
  else {

    /* we need a new interval */
    iv_newr = _unur_tdr_interval_new( gen, x, fx, FALSE );
    CHECK_NULL(iv_newr,0);     /* case of error */
    
    /* link into list */
    iv_newr->next = iv_oldl->next;
    iv_oldl->next = iv_newr;
    
    /* compute parameters for interval */
    if ( _unur_tdr_interval_parameter(gen, iv_oldl) <= 0 ||
	 _unur_tdr_interval_parameter(gen, iv_newr) <= 0 ) {
      
      /* p.d.f. not T-concave, or new interval to narrow, 
	 or area below hat not bounded */

      /* new construction point not suitable --> do not add */
      _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"Cannot split interval at given point.");
#if UNUR_DEBUG & UNUR_DB_INFO
      /* write info into log file */
      if (gen->debug & TDR_DB_SPLIT) 
	_unur_tdr_debug_split_stop( gen,iv_oldl,iv_newr );
#endif
      
      /* remove from linked list */
      iv_oldl->next = iv_newr->next;
      _unur_tdr_iv_stack_push(gen);
      /* we have to restore the old interval.
	 (this case should not happen, so it is faster not to make a 
	 backup of the old interval) */
      if ( !_unur_tdr_interval_parameter(gen, iv_oldl) ) {
	_unur_error(gen->genid,UNUR_ERR_SAMPLE,"Cannot restore interval. PANIK.");
	exit (-1);
      }
      return 0;
    }
  }
  
  /* update guide table */ 
  _unur_tdr_make_guide_table(gen);
  
#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug & TDR_DB_SPLIT) 
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
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

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
	_unur_warning(gen->genid,UNUR_ERR_INIT,"roundoff error while making guide table!");
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

/*****************************************************************************/

static struct unur_tdr_interval *
_unur_tdr_iv_stack_pop( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* pop free interval from stack; allocate memory block if necessary.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to interval                                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* look for an unused segment */
  if( ! GEN.iv_free ) {
    /* no more unused segments. make some. */
    GEN.iv_stack = _unur_malloc( UNUR_MALLOC_SIZE * sizeof(struct unur_tdr_interval) );

    /* reset counter */
    GEN.iv_free = UNUR_MALLOC_SIZE;

    /* set cookies */
    COOKIE_SET_ARRAY( GEN.iv_stack, CK_TDR_IV, UNUR_MALLOC_SIZE);

    /* add to list of allocated blocks */
    _unur_add_mblocks( &(GEN.mblocks), GEN.iv_stack);

  }

  /* update ....                                   */
  --(GEN.iv_free);   /* pointer to free segments  */
  ++(GEN.n_ivs);     /* counter for used segments */

  /* return pointer to segment */
  return (GEN.iv_stack + GEN.iv_free);

} /* end of _unur_tdr_iv_stack_pop() */

/*---------------------------------------------------------------------------*/

static void
_unur_tdr_iv_stack_push( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* push useless segment back stack                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  /* update counters and pointers */
  --(GEN.n_ivs);
  ++(GEN.iv_free);
} /* end of _unur_tdr_iv_stack_push() */

/*-----------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

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
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(par,/*void*/);
  COOKIE_CHECK(par,CK_TDR_PAR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = transformed density rejection\n",gen->genid);
  fprintf(log,"%s: transformation T_c(x) = ",gen->genid);
  switch( gen->method & TDR_MASK_T ) {
  case TDR_METH_LOG:
    fprintf(log,"log(x)  ... c = 0");                   break;
  case TDR_METH_SQRT:
    fprintf(log,"-1/sqrt(x)  ... c = -1/2");            break;
  case TDR_METH_POW:
    fprintf(log,"-x^(%g)  ... c = %g",PAR.c_T,PAR.c_T); break;
  }
  _unur_print_if_default(par,UNUR_SET_TDR_C);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_tdr_sample_",gen->genid);
  if (par->method & UNUR_MASK_SCHECK)
    fprintf(log,"check()\n");
  else
    switch( gen->method & TDR_MASK_T ) {
    case TDR_METH_LOG:   
      fprintf(log,"log()\n");   break;
    case TDR_METH_SQRT:
      fprintf(log,"sqrt()\n");  break;
    case TDR_METH_POW:
      fprintf(log,"pow()\n");   break;
    }
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: p.d.f with %d arguments\n",gen->genid,PAR.n_pdf_param);
  if (PAR.n_pdf_param)
    for( i=0; i<PAR.n_pdf_param; i++ )
      fprintf(log,"%s:\tparam[%d] = %g\n",gen->genid,i,PAR.pdf_param[i]);
  fprintf(log,"%s:\n",gen->genid);

  if (par->method & UNUR_MASK_MODE )
    fprintf(log,"%s: mode = %g\n",gen->genid,PAR.mode);
  else
    fprintf(log,"%s: mode unknown\n",gen->genid);
  fprintf(log,"%s: domain = (%g, %g)",gen->genid,GEN.bleft,GEN.bright);
  _unur_print_if_default(par,UNUR_SET_DOMAIN);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR.max_ivs);
  _unur_print_if_default(par,UNUR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,UNUR_SET_MAX_RATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
  _unur_print_if_default(par,UNUR_SET_FACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,UNUR_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & UNUR_SET_STP)
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
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

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
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:Intervals: %d\n",gen->genid,GEN.n_ivs);
  if (gen->debug & TDR_DB_IV) {
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

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }
    
  /* print and sum areas below squeeze and hat */
  if (gen->debug & TDR_DB_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\tbelow squeeze\t\t  below hat (left and right)\t\t  cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
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
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  COOKIE_CHECK(iv,CK_TDR_IV,/*void*/);
  COOKIE_CHECK(pt,CK_TDR_IV,/*void*/);

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

/*****************************************************************************/
#endif
