/*****************************************************************************
 *                                                                           *
 *          unuran -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      arou.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    ratio-of-uniforms with enveloping polygon                    *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f of a T-concave distribution;                             *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Leydold J. (2000): Automatic Sampling with the ratio-of-uniforms    *
 *       method, ACM TOMS, forthcoming                                       *
 *                                                                           *
 *   [2] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [2] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (V,U) is uniformly distributed in                                   *
 *       A = {(v,u): 0 < u <= sqrt(g(v/u)), x_0 < v/u < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (V,U) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 * The implemented algorithm uses the fact, that for many distribtions,      *
 * A is convex. Then we easily can construct an enveloping polygon by means  *
 * of tangent lines and a squeeze region by means of secants.                *
 * The resulting algorithm is very fast, since we can sample from the        *
 * squeeze region with immedate acceptance (thus only one uniform random     *
 * number is necessary). The region between envelope and squeeze consists    *
 * of triangles and can be made arbitrarily small (see [1] for details).     *
 *                                                                           *
 * Distributions with a convex set A are characterized by the following      *
 * theorem that shows a connection to transformed density rejection TDR.     *
 *                                                                           *
 * THEOREM:                                                                  *
 *    A is convex if and only if g is T-concave with transformation          *
 *    T(x) = -1/sqrt(x), i.e., -1/sqrt(g(x)) is a concave function.          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * The region A is divided into segments, that contain the origin. Each      *
 * segments consist of an inner triangle (the squeeze region) and the        *
 * outer triangle (the region between envelope and squeeze). We have to      *
 * compute the areas of these triangles.                                     *
 *                                                                           *
 * To generate from the distribution we have to sample from a discrete       *
 * random variate with probability vector proportional to these areas to get *
 * one of these triangles (inner or outer). We use indexed search (or guide  *
 * tables) to perform this task ([1], see also description of DIS).          *
 * When we have an inner triangle (squeeze), we reuse the uniform random     *
 * variate to get a point (v,u) uniformly distributed on the edge opposite   *
 * to the origin and return the ratio x = v/u (Thus generation from the      *
 * squeeze region is equivalent to the inversion method.)                    *
 * When whe have an outer triangle, we have to sample a point (v,u)          *
 * uniformly distributed in this triangle. If u <= g(v/u) then return the    *
 * ratio x = v/u, otherwise reject.                                          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Algorithm AROU                                                            *
 *                                                                           *
 * [Required]                                                                *
 * p.d.f. f(x), construction points c_1,...,c_n                              *
 *                                                                           *
 * [Setup]                                                                   *
 *  1: Construct inner triangles S_i^s and outer triangles S_i^o             *
 *  2: Foreach triangle Do                                                   *
 *  3:    Compute respective areas A_i^s and A_i^o.                          *
 *                                                                           *
 * [Generate]                                                                *
 *  4: Generate I proportional to (A_1^s, A_1^o; A_2^s, A_2^o; ...).         *
 *  5: If inner triangle S_i^s Then                                          *
 *  6:    Compute point (v,u) on edge (Reuse u.r.n. from 4).                 *
 *  7:    Return v/u.                                                        *
 *  8: Else                                                                  *
 *  9:    Generate point (V,U) uniformly distributed in outer triangle S_i^o.*
 * 10:    If u <= f(v/u) Return V/U.                                         *
 * 11:    Goto 4.                                                            *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define AROU_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */
#define AROU_VARFLAG_USECENTER  0x02u   /* flag whether center is used as cpoint or not */

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

#define AROU_DEBUG_SPLIT        0x100u   /* trace splitting of segments      */
#define AROU_DEBUG_SEGMENTS     0x200u   /* print list of segments           */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define AROU_SET_CENTER         0x001u
#define AROU_SET_STP            0x002u
#define AROU_SET_N_STP          0x004u
#define AROU_SET_GUIDEFACTOR    0x010u
#define AROU_SET_MAX_SQHRATIO   0x020u
#define AROU_SET_MAX_SEGS       0x040u

/*---------------------------------------------------------------------------*/

#define GENTYPE "AROU"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_arou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_arou_get_starting_segments( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute segments from given starting construction points.                 */
/*---------------------------------------------------------------------------*/

static struct unur_arou_segment *_unur_arou_segment_new( struct unur_gen *gen, double x, double fx );
/*---------------------------------------------------------------------------*/
/* make a new segment with left construction point x = v/u.                  */
/*---------------------------------------------------------------------------*/

static int _unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for segment.                                   */
/* return 0 if p.d.f. is not T-concave.                                      */
/*---------------------------------------------------------------------------*/

static int _unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_old, double x, double fx );
/*---------------------------------------------------------------------------*/
/* split a segment at point (direction) x. return 0 if not successful.       */                                           
/*---------------------------------------------------------------------------*/

static int _unur_arou_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static double _unur_arou_segment_arcmean( struct unur_arou_segment *seg );
/*---------------------------------------------------------------------------*/
/* compute the "arcmean" of the two construction points of a segement.       */
/*---------------------------------------------------------------------------*/

static struct unur_arou_segment *_unur_arou_segment_stack_pop( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* pop a segment from the stack of free segments.                            */
/*---------------------------------------------------------------------------*/

static void _unur_arou_segment_stack_push( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* push the last popped segment back onto the stack.                         */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_segments( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for segments.                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_split_start( struct unur_gen *gen, struct unur_arou_segment *seg, double x, double fx );
static void _unur_arou_debug_split_stop( struct unur_gen *gen, 
					 struct unur_arou_segment *seg_left, struct unur_arou_segment *seg_right );
/*---------------------------------------------------------------------------*/
/* print before and after a segment has been split (not / successfully).     */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_printratio( double v, double u, char *string );
/*---------------------------------------------------------------------------*/
/* print the ratio of two double if possible, and Inf or NaN otherwise.      */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR   distr->data.cont
#define PAR     par->data.arou
#define GEN     gen->data.arou
#define SAMPLE  gen->sample.cont

#define PDF(x) ((*(GEN.pdf))((x),GEN.pdf_param,GEN.n_pdf_param))
#define dPDF(x) ((*(GEN.dpdf))((x),GEN.pdf_param,GEN.n_pdf_param))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_arou_new( struct unur_distr *distr )
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
  CHECK_NULL(distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.d.f.");
    return NULL;
  }
  if (DISTR.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of p.d.f.");
    return NULL;
  }

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_AROU_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR.guide_factor        = 3.;     /* size of guide table / number of intervals */

  PAR.starting_cpoints    = NULL;   /* pointer to array of starting points   */
  PAR.n_starting_cpoints  = 30;     /* number of starting points             */
  PAR.max_segs            = 50;     /* maximum number of segments            */
  PAR.max_ratio           = 0.95;   /* do not add construction points if
				       ratio r_n = |P^s| / |P^e| > max_ratio */
  PAR.bound_for_adding    = 0.5;    /* do not add a new construction point in a segment, 
				       where abiguous region is too small, i.e. if 
				       the area / (|S^e\S^s|/number of segments) < bound_for_adding */

  par->method   = UNUR_METH_AROU;          /* method                         */
  par->variant  = AROU_VARFLAG_USECENTER;  /* default variant                */
  par->set      = 0u;                      /* inidicate default parameters   */    
  par->urng     = unur_get_default_urng(); /* use default urng               */

  _unur_set_debugflag_default(par); /* set default debugging flags           */

  /* we use the mode (if known) as center of the distribution */
  if (distr->set & UNUR_DISTR_SET_MODE) {
    PAR.center = DISTR.mode;
    par->set |= AROU_SET_CENTER;
  }
  else
    PAR.center = 0.;        /* the default */

  /* routine for starting generator */
  par->init = unur_arou_init;

  return par;

} /* end of unur_arou_new() */

/*****************************************************************************/

int
unur_arou_set_cpoints( struct unur_par *par, int n_stp, double *stp )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"number of starting points < 0");
    return 0;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_SET,"starting points not strictly monotonically increasing");
	return 0;
      }

  /* store date */
  PAR.starting_cpoints = stp;
  PAR.n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= AROU_SET_N_STP | ((stp) ? AROU_SET_STP : 0);

  return 1;

} /* end of unur_arou_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_guidefactor( struct unur_par *par, double factor )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"relative table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= AROU_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_arou_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_max_sqhratio( struct unur_par *par, double max_ratio )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"ratio Atotal / Asqueeze not in [0,1]");
    return 0;
  }

  /* store date */
  PAR.max_ratio = max_ratio;

  /* changelog */
  par->set |= AROU_SET_MAX_SQHRATIO;

  return 1;

} /* end of unur_arou_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_max_segments( struct unur_par *par, int max_segs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of segments                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_segs  ... maximum number of segments                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* check new parameter for generator */
  if (max_segs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"maximum number of segments < 1");
    return 0;
  }

  /* store date */
  PAR.max_segs = max_segs;

  /* changelog */
  par->set |= AROU_SET_MAX_SEGS;

  return 1;

} /* end of unur_arou_set_max_segments() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_center( struct unur_par *par, double center )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* store data */
  PAR.center = center;

  /* changelog */
  par->set |= AROU_SET_CENTER;

  /* o.k. */
  return 1;

} /* end of unur_arou_set_center() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_usecenter( struct unur_par *par, int usecenter )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* we use a bit in variant */
  par->variant = (usecenter) ? (par->variant | AROU_VARFLAG_USECENTER) : (par->variant & (~AROU_VARFLAG_USECENTER));

  /* o.k. */
  return 1;

} /* end of unur_arou_set_usecenter() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_verify( struct unur_par *par, int verify )
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
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( AROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | AROU_VARFLAG_VERIFY) : (par->variant & (~AROU_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_arou_set_verify() */

/*****************************************************************************/

struct unur_gen *
unur_arou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_AROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_AROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_arou_create(par);
  if (!gen) { free(par); return NULL; }

  /* get starting points */
  if (!_unur_arou_get_starting_cpoints(par,gen) ) {
    free(par); unur_arou_free(gen);
    return NULL;
  }

  /* compute segments for given starting points */
  if ( !_unur_arou_get_starting_segments(par,gen) ) {
    free(par); unur_arou_free(gen);
    return NULL;
  }

  /* we have to update the maximal number of segments,
     if the user wants more starting points. */
  if (GEN.n_segs > GEN.max_segs) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"maximal number of segments too small. increase.");
    GEN.max_segs = GEN.n_segs;
  }

  /* make initial guide table */
  _unur_arou_make_guide_table(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_arou_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  /* is there any envelope at all ? */
  if (GEN.Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_INIT_FAILED,"cannot construct envelope. bad construction points.");
    unur_arou_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of unur_arou_init() */

/*****************************************************************************/

double
unur_arou_sample( struct unur_gen *gen )
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
  /** TODO: check uniform random number: u != 0 and u != 1 ??  **/

  struct unur_arou_segment *seg;
  double R,R1,R2,R3,tmp,x,fx,u;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_AROU_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    R = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    seg =  GEN.guide[(int) (R * GEN.guide_size)];
    R *= GEN.Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,0.);

    /* reuse of uniform random number */
    R = seg->Acum - R;

    /* inside or outside squeeze */
    if (R < seg->Ain) {
      /* inside */
      /* reuse of random number.                       */
      /* We can avoid R = (seg->Ain - R) / seg->Ain    */
      return( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	      ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );
    }

    else {
      /* outside */

      /* three uniform random numbers with R1 + R2 + R3 = 1 */
      R1 = (R - seg->Ain) / seg->Aout;  /* reuse of random number (good ?? ) */
      R2 = _unur_call_urng(gen);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  /* swap */
      R3 = 1.-R2;
      R2 -= R1;

      /* point (v,u) and ratio x = v/u */
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;

      /* density at x */
      fx = PDF(x);

      /* being outside the squeeze is bad. improve the situation! */
      if (GEN.n_segs < GEN.max_segs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze)
	_unur_arou_segment_split(gen,seg,x,fx);

      /* if inside region of acceptance, return ratio x */
      if (u*u <= fx) 
	return x;
    }
  }
} /* end of unur_arou_sample() */

/*****************************************************************************/

double
unur_arou_sample_check( struct unur_gen *gen )
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
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  /** TODO: check uniform random number: u != 0 and u != 1 ??  **/

  struct unur_arou_segment *seg;
  double R,R1,R2,R3,tmp,x,fx,u,sqx,a;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_AROU_GEN,0.);

  while (1) {

    /* sample from U(0,1) */
    R = _unur_call_urng(gen);

    /* look up in guide table and search for segment */
    seg =  GEN.guide[(int) (R * GEN.guide_size)];
    R *= GEN.Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,0.);

    /* reuse of uniform random number */
    R = seg->Acum - R;

    /* inside or outside squeeze */
    if (R < seg->Ain) {
      /* inside */
      /* reuse of random number.                       */
      /* We can avoid R = (seg->Ain - R) / seg->Ain    */
      x = ( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	    ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );

      /* density at x */
      fx = PDF(x);

      /* compute value of squeeze at x, i.e., we have to solve 
	 a*ltp[0] + (1-a)*rtp[0] == a*ltp[1] + (1-a)*rtp[1] */
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];

      /* test for T-concavity */
      if (sqx*sqx > fx)
	_unur_error(gen->genid,UNUR_ERR_SAMPLE,"p.d.f. not T-concave.");

      return x;
    }

    else {
      /* outside */

      /* three uniform random numbers with R1 + R2 + R3 = 1 */
      R1 = (R - seg->Ain) / seg->Aout;  /* reuse of random number (good ?? ) */
      R2 = _unur_call_urng(gen);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  /* swap */
      R3 = 1.-R2;
      R2 -= R1;

      /* point (v,u) and ratio x = v/u */
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;

      /* density at x */
      fx = PDF(x);

      /* compute value of squeeze at x, i.e., we have to solve 
	 a*ltp[0] + (1-a)*rtp[0] == a*ltp[1] + (1-a)*rtp[1] */
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];

      /* test for T-concavity */
      if (sqx*sqx > fx)
	_unur_error(gen->genid,UNUR_ERR_SAMPLE,"p.d.f. not T-concave.");

      /* being outside the squeeze is bad. improve the situation! */
      if (GEN.n_segs < GEN.max_segs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze)
	_unur_arou_segment_split(gen,seg,x,fx);

      /* if inside region of acceptance, return ratio x */
      if (u*u <= fx) 
	return x;
    }
  }

} /* end of unur_arou_sample_check() */

/*****************************************************************************/

void
unur_arou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_AROU ) {
    _unur_warning(GENTYPE,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_AROU_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  if (gen->debug) _unur_arou_debug_free(gen);
#endif

  /* free linked list of segments and others */
  _unur_free_mblocks(GEN.mblocks);

  /* free other memory not stored in list */
  _unur_free_genid(gen);
  free(GEN.guide);
  free(gen);

} /* end of unur_arou_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_arou_create( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_AROU_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_AROU_GEN);

  /* set generator identifier */
  _unur_set_genid(gen,GENTYPE);

  /* copy pointer to distribution object */
  /* (we do not copy the entire object)  */
  gen->distr = par->distr;

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & AROU_VARFLAG_VERIFY) ? unur_arou_sample_check : unur_arou_sample;
  gen->destroy = unur_arou_free;

  /* set all pointers to NULL */
  GEN.seg         = NULL;
  GEN.n_segs      = 0;
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.seg_stack   = NULL;
  GEN.seg_free    = 0;
  GEN.mblocks     = NULL;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN.pdf = gen->DISTR.pdf;           /* p.d.f. of distribution              */
  GEN.dpdf = gen->DISTR.dpdf;         /* derivative of p.d.f.                */

  GEN.pdf_param   = gen->DISTR.params;
  GEN.n_pdf_param = gen->DISTR.n_params;

  GEN.bleft = gen->DISTR.domain[0];   /* left boundary of domain             */
  GEN.bright = gen->DISTR.domain[1];  /* right boundary of domain            */

  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN.max_segs = PAR.max_segs;      /* maximum number of segments            */
  GEN.max_ratio = PAR.max_ratio;    
  GEN.bound_for_adding = PAR.bound_for_adding;

  gen->method = par->method;        /* indicates method and variant          */
  gen->variant = par->variant;      /* indicates variant                     */
  _unur_copy_urng_pointer(par,gen); /* copy pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);    /* copy debugging flags into generator object */

  /* center known ?? */
  if (!(par->set & AROU_SET_CENTER))
    /* we cannot use the center as construction point */
    par->variant = par->variant & (~AROU_VARFLAG_USECENTER);
  else {
    /* center must be in domain */
    PAR.center = max(PAR.center,GEN.bleft);
    PAR.center = min(PAR.center,GEN.bright);
  }

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_arou_create() */

/*****************************************************************************/

static int
_unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting segments.                   */
     /* if not provided as arguments compute these                           */
     /* by means of the "equi-angle rule".                                   */
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
  struct unur_arou_segment *seg, *seg_new;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int i, use_center, is_center, was_center, is_increasing;

  /* check arguments */
  COOKIE_CHECK(par,CK_AROU_PAR,0);
  COOKIE_CHECK(gen,CK_AROU_GEN,0);

  /* initialize boolean */
  is_center = was_center = FALSE;

  /* use center as construction point ? */
  use_center = (par->variant & AROU_VARFLAG_USECENTER) ? TRUE : FALSE;

  /* reset counter of segments */
  GEN.n_segs = 0;

  /* prepare for computing construction points */
  if (!PAR.starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  ( GEN.bleft <= -INFINITY ) ? -M_PI/2. : atan(GEN.bleft - PAR.center);  
    right_angle = ( GEN.bright >= INFINITY )  ? M_PI/2.  : atan(GEN.bright - PAR.center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR.n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = GEN.bleft;
  fx = fx_last = (x <= -INFINITY) ? 0. : PDF(x);
  seg = GEN.seg = _unur_arou_segment_new( gen, x, fx );
  CHECK_NULL(seg,0);       /* case of error */
  is_increasing = 1;       /* assume pdf(x) is increasing for the first construction points */

  /* now all the other points */
  for( i=0; i<=PAR.n_starting_cpoints; i++ ) {
    was_center = is_center;

    /* starting point */
    if (i < PAR.n_starting_cpoints) {
      if (PAR.starting_cpoints) {   
	/* construction points provided by user */
	x = PAR.starting_cpoints[i];
	/* check starting point */
	if (x <= GEN.bleft || x >= GEN.bright) {
	  _unur_warning(gen->genid,UNUR_ERR_INIT,"starting point out of domain!");
	  continue;
	}
	if (x<=x_last) {
	  _unur_warning(gen->genid,UNUR_ERR_INIT,"starting points are not strictly monotonically increasing! skip!");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equidistance" rule */
	angle += diff_angle;                /** TODO: angle >= M_PI/2. !! **/
	x = tan( angle ) + PAR.center;      /** TODO: possible over/underflow **/
      }
    }
    else {
      /* the very last segment. it is rather a "virtual" segment to store 
	 the right vertex of the last segment, i.e., the right boundary point. */
      x = GEN.bright;
    }

    /* insert center ? */
    if (use_center && x >= PAR.center) {
      use_center = FALSE;   /* we use the center only once (of course) */
      is_center = TRUE;     /* the next construction point is the center */
      if (x>PAR.center) {
	x = PAR.center;   /* use the center now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!PAR.starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x == PAR.center --> nothing to do */
    }
    else
      is_center = FALSE;

    /** TODO: check if two construction points are too close ??
	check if a point is too close to center ??  */

    /* value of p.d.f. at starting point */
    fx = (x >= INFINITY) ? 0. : PDF(x);

    /* check value of p.d.f. at starting point */
    if (!is_increasing && fx > fx_last) {
      _unur_error(gen->genid,UNUR_ERR_INIT_FAILED,"p.d.f. not unimodal!");
      return 0;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such points */
      if (is_increasing) {
	/* p.d.f. is still increasing, i.e., constant 0 til now */
	if (i<PAR.n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the p.d.f. is constant 0 on all construction points.
	     then we need both boundary points. */
	  /* we only have to change tangent line v/u = x,
	     everything else remains unchanged */
	  seg->dltp[1] = x;
	  /* seg->dltp[0] = -1; seg->dltp[2] = 0.;  not changed */
	  x_last = x;
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with pdf(x) > 0 */
	break;
    }
    
    /* need a new segment */
    seg_new = _unur_arou_segment_new( gen, x, fx );
    CHECK_NULL(seg,0);     /* case of error */

    /* append to linked list */
    seg->next =seg_new;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;
    
    /* next step */
    seg = seg_new;

    /* p.d.f. still increasing ? */
    if (is_increasing && fx < fx_last)
      is_increasing = 0;

    /* store last computed values */
    x_last = x;
    fx_last = fx;
  }

  /* we have left the loop with the right boundary of the support of p.d.f.
     make shure that we will never use seg for sampling. */
  seg->Ain = seg->Aout = 0.;
  seg->Acum = INFINITY;
  seg->next = NULL;         /* terminate list */
  --(GEN.n_segs);           /* we do not count this segment */

  /* o.k. */
  return 1;

} /* end of _unur_arou_get_starting_cpoints() */

/*****************************************************************************/

static int
_unur_arou_get_starting_segments( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute segments for starting points                                 */
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
  struct unur_arou_segment *seg, *seg_new; 
  double x,fx;              /* construction point, value of p.d.f. at x */

  /* check arguments */
  COOKIE_CHECK(par,CK_AROU_PAR,0);
  COOKIE_CHECK(gen,CK_AROU_GEN,0);

  /* compute paramters for all segments */
  for( seg=GEN.seg; seg->next != NULL; ) {

    /* compute parameters for semgent */
    switch (_unur_arou_segment_parameter(gen, seg)) {
    case 0:     /* p.d.f. not T-concave */
      return 0;
    case 1:     /* computation of parameters for segment successful */
      /* skip to next segment. */
      seg = seg->next;
      continue;
    case -1:    /* segment unbounded */
      /* split segment */
      break;
    }

    /* area in segment infinite. insert new construction point. */
    x = _unur_arou_segment_arcmean(seg);  /* use mean point in segment */

    /* value of p.d.f. at x */
    fx = PDF(x);

    /* add a new segment, but check if we had to used too many segments */
    if (GEN.n_segs >= GEN.max_segs) {
      /* we do not want to create too many segments */
      _unur_error(gen->genid,UNUR_ERR_INIT_FAILED,"cannot create bounded envelope!");
      return 0;
    }
    seg_new = _unur_arou_segment_new( gen, x, fx );
    CHECK_NULL(seg_new,0);     /* case of error */
    /* insert into linked list */
    seg_new->next = seg->next;
    seg->next = seg_new;
    /* right vertices */
    seg_new->rtp = seg->rtp;
    seg_new->drtp = seg->drtp;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;

  }

  /* o.k. */
  return 1;

} /* end of _unur_arou_get_starting_segments() */

/*****************************************************************************/

static struct unur_arou_segment *
_unur_arou_segment_new( struct unur_gen *gen, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* get new segment and compute left construction point at x.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left point of new segment                                  */
     /*   fx  ... value of p.d.f. at x                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new segment                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg;
  double u,v,dfx;

  /* check arguments */
  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_warning(gen->genid,UNUR_ERR_INIT,"pdf(x) < 0.");
    return NULL;
  }

  /* we need new segment */
  seg = _unur_arou_segment_stack_pop(gen);
  COOKIE_CHECK(seg,CK_AROU_SEG,NULL); 

  /* make left construction point in segment */

  /* case: x out of support */
  if (fx == 0. ) {
    seg->ltp[0] = 0.;   /* vertex == origin */
    seg->ltp[1] = 0.;
    if (x <= -INFINITY || x >= INFINITY ) {
      /* tangent line == line u = 0 (i.e., v-axis) */
      seg->dltp[0] = 0.;   /* dv */
      seg->dltp[1] = 1.;   /* du */
      seg->dltp[2] = 0.;   /* v * dv + u * du */
    }
    else {
      /* tangent line == line v/u = x */
      seg->dltp[0] = -1.;  /* dv */
      seg->dltp[1] = x;    /* du */
      seg->dltp[2] = 0.;   /* v * dv + u * du */
    }
    return seg;
  }

  /* case: x in support */
  /* boundary point */
  u = sqrt( fx );
  v = x * u;
  seg->ltp[0] = v;
  seg->ltp[1] = u; 

  /* tangent line at tp */

  /* compute derivative of p.d.f. at tp x */
  dfx = dPDF(x);

  /* subcase: derivative bounded 
     use derivative for tangent line */
  if ( dfx > -INFINITY && dfx < INFINITY ) {
    seg->dltp[0] = -dfx / u;             /* dv */    /** TODO: possible overflow **/
    seg->dltp[1] = 2 * u + dfx * x / u;  /* du */    /** TODO: possible overflow **/
    seg->dltp[2] = seg->dltp[0] * v + seg->dltp[1] * u;
    return seg;
  }

  /* subcase: derivative unbounded.
     use straight line through origin and vertex */
  seg->dltp[0] = -u;   /* dv */
  seg->dltp[1] = v;    /* du */
  seg->dltp[2] = 0.;
  return seg;

} /* end of _unur_arou_segment_new() */

/*****************************************************************************/

#define MAX_NORM_INTERSECTION  1.e6    /* maximal distance of intersection point
					  from origin compared to distance of
					  construction points to origin      */

/*---------------------------------------------------------------------------*/

static int
_unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg )
     /*----------------------------------------------------------------------*/
     /* compute all parameters for a segment.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   seg ... pointer to segment                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*  -1 ... area = INFINITY                                              */
     /*   0 ... error (not p.d.f. T-concave)                                 */
     /*----------------------------------------------------------------------*/
{
  /** TODO: roundoff errors ?? **/ 

  double coeff_det, cramer_det[2];
  double norm_vertex;
  double det_bound;

  /* check arguments */
  COOKIE_CHECK(gen,CK_AROU_GEN,0);
  COOKIE_CHECK(seg,CK_AROU_SEG,0);

  /* area inside the squeeze */
  seg->Ain = (seg->ltp[1] * seg->rtp[0] - seg->ltp[0] * seg->rtp[1]) / 2.;
  /* due to our ordering of construction points, seg->Ain must be >= 0 ! */
  if( seg->Ain < 0. ) {
    /* This should not happen ! */
    _unur_error(gen->genid,UNUR_ERR_INIT_FAILED,"non-ascending ordering of construction points");
    return 0;
  }

  /* we use Cramer's rule to compute intersection of tangent lines.
     (well, we could save one multiplication otherwise)             */
  coeff_det     = seg->dltp[0] * seg->drtp[1] - seg->dltp[1] * seg->drtp[0];
  cramer_det[0] = seg->dltp[2] * seg->drtp[1] - seg->dltp[1] * seg->drtp[2];
  cramer_det[1] = seg->dltp[0] * seg->drtp[2] - seg->dltp[2] * seg->drtp[0];

  /* we there are two possibilities for singular coefficent matrix:

     either the outer triangle is unlimited.
            then the two tangents are distinct but parallel and 
	    the corresponding linear equation has no solution, i.e.
	    coeff_det == 0 but cramer_det[0] != 0 or cramer_det[1] != 0.

     or     the outer triangle degenerates to a line segment and has area 0.
            then the two tangents are equal and
	    the corresponding linear equation has a nonunique solution, i.e.
	    coeff_det == cramer_det[0] == cramer_det[1] == 0.
  */

  /* sum of 1-norm of vertices (notice that u > 0 by assumption) */
  norm_vertex = fabs(seg->ltp[0]) + seg->ltp[1] + fabs(seg->rtp[0]) + seg->rtp[1];

  /* we to not allow that the outer triangles becomes too large.
     so if the sup-norm of intersection point is too large compared
     to norm_vertex we assume that this triangle is unbounded.
     we thus avoid numerical errors.
     (we use the sup-norm here since it much easier to handle.)

     (However this might also happen due to roundoff errors,
     when the real position is extremely close to the secant.
     But at least we are on the save side.)
  */
  det_bound = fabs(coeff_det) * norm_vertex * MAX_NORM_INTERSECTION;

  if ( fabs(cramer_det[0]) > det_bound || fabs(cramer_det[1]) > det_bound ) {
    /* case: triangle is assumed to be unbounded */	     
/*      _unur_warning(gen->genid,UNUR_ERR_INIT,"outer triangle assumed unbounded"); */
    seg->Aout = INFINITY;
    return -1;
  }

  /* case: intersection point exists and is unique */
  if (coeff_det != 0.) {

    /* compute intersection point */
    seg->mid[0] = cramer_det[0] / coeff_det;
    seg->mid[1] = cramer_det[1] / coeff_det;

    /* area outside the squeeze */
    seg->Aout = ( (seg->ltp[0] - seg->mid[0]) * (seg->rtp[1] - seg->mid[1])
		  - (seg->ltp[1] - seg->mid[1]) * (seg->rtp[0] - seg->mid[0])) / 2.;
    /* due to our ordering of construction points, seg->Aout must be >= 0
       for a regular triangle.
       Thus if seg->aout < 0, then the intersection point of tangents is on
       the WRONG side of the secant through vertices of segment,
       i.e. the "triangle" outside the squeeze region is unbounded.

       However this might also happen due to roundoff errors when 
       the real position is extremely close to the secant.
       We can distinguish between these two case by means of the u-coordinate
       of the intersection point. If it really is on the wrong side if secant 
       then seg->mid[1] < 0.
    */
    if( seg->mid[1] < 0. ) {
/*        _unur_warning(gen->genid,UNUR_ERR_INIT,"outer triangle unbounded"); */
      seg->Aout = INFINITY;
      return -1;
    }

    /* at last check result.
       we must have:
         (*) seg->aout > 0
         (*) intersection point right of left construction point
         (*) intersection point left of right construction point
    */
    if ( seg->Aout > 0. &&
	 seg->mid[0] * seg->ltp[1] >= seg->ltp[0] * seg->mid[1] &&
	 seg->mid[0] * seg->rtp[1] <= seg->rtp[0] * seg->mid[1]	) {
      /* everything o.k. */
      return 1;
    }

    /* there are two case is a second reason why the above check failed:
       (1) the p.d.f. is not T-concave
       (2) small roundoff errors.
    */
    /** TODO: check for roundoff-errors !!! **/
    _unur_error(gen->genid,UNUR_ERR_INIT_FAILED,"p.d.f. not T-concave");
    return 0;
  }

  /* remaining case: triangle degenerates to a line segment, i.e.
                     intersection point exists but is not unique */

  /* boundary of region is (almost) a straight line
     and area outside squeeze is (almost) zero.
     use middle point as intersection point and 
     set area outside the squeeze to 0.
  */
/*    _unur_warning(gen->genid,UNUR_ERR_INIT,"outer triangle is line"); */
  seg->mid[0] =  0.5 * (seg->ltp[0] + seg->rtp[0]);
  seg->mid[1] =  0.5 * (seg->ltp[1] + seg->rtp[1]);
  seg->Aout = 0.;
  
  /* now it should be o.k. */
  return 1;
  
} /* end if _unur_arou_segment_parameter() */

/*---------------------------------------------------------------------------*/

#undef MAX_NORM_INTERSECTION

/*****************************************************************************/

static int
_unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_oldl, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* insert new segment                                                   */
     /*   old segment -> left hand side                                      */
     /*   new segment -> right hand side                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   seg_oldl ... pointer to segment                                    */
     /*   x        ... left point of new segment                             */
     /*   fx       ... value of p.d.f. at x                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1  ... if successful                                               */
     /*   0  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg_newr;
  double backup;

  /* check arguments */
  COOKIE_CHECK(gen,CK_AROU_GEN,0);
  COOKIE_CHECK(seg_oldl,CK_AROU_SEG,0);

#if UNUR_DEBUG & UNUR_DB_INFO
    /* write info into log file */
    if (gen->debug & AROU_DEBUG_SPLIT) 
      _unur_arou_debug_split_start( gen,seg_oldl,x,fx );
#endif

  /* we only add a new construction point, if the relative area is large enough */
  if (GEN.n_segs * seg_oldl->Aout / (GEN.Atotal - GEN.Asqueeze) < GEN.bound_for_adding )
    return 0;

  /* p.d.f. at x is 0. */
  if (fx <= 0.) {
    if (seg_oldl->rtp[1] <= 0. && seg_oldl->rtp[0] <= 0. ) {
      /* just chop off the right part of segment */
      /* we only have to change tangent line v/u = x
	 at the right hand vertex */
      backup = seg_oldl->drtp[1];
      seg_oldl->drtp[1] = x;    /* du */
    }
    else if (seg_oldl->ltp[1] <= 0. && seg_oldl->ltp[0] <= 0. ) {
      /* just chop off the left part of segment */
      /* we only have to change tangent line v/u = x
	 at the left hand vertex */
      backup = seg_oldl->dltp[1];
      seg_oldl->dltp[1] = x;    /* du */
    }
    else {
      _unur_warning(gen->genid,UNUR_ERR_ADAPT,"This should not happen!!");
      return 0;
    }
    
    /* parameters of new segment */
    if( _unur_arou_segment_parameter(gen,seg_oldl) <= 0 ) {
      /* p.d.f. not T-concave or area in segment not bounded */

      /* error, restore */
      _unur_warning(gen->genid,UNUR_ERR_ADAPT,"Cannot chop segment at given point");
      seg_oldl->drtp[1] = backup;
      if ( !_unur_arou_segment_parameter(gen,seg_oldl) ) {
	_unur_error(gen->genid,UNUR_ERR_ADAPT,"Cannot restore segment. PANIK.");
	exit (-1);
      }
      return 0;
    }

    /* for _unur_arou_debug_split_stop only */
    seg_newr = seg_oldl;
  }

  else {  /* fx > 0 */

    /* need new segment */
    seg_newr = _unur_arou_segment_new(gen,x,fx);
    CHECK_NULL(seg_newr,0);     /* case of error */
    
  /* link into list */
    seg_newr->next = seg_oldl->next;
    seg_oldl->next = seg_newr;
    
    /* right vertices */
    seg_newr->rtp = seg_oldl->rtp;
    seg_newr->drtp = seg_oldl->drtp;
    seg_oldl->rtp = seg_newr->ltp;
    seg_oldl->drtp = seg_newr->dltp;
    
    /* parameters of new segments */
    if( _unur_arou_segment_parameter(gen,seg_oldl) <= 0 ||       /* p.d.f. not T-concave or */
	_unur_arou_segment_parameter(gen,seg_newr) <= 0  ) {     /* area in segment not bounded */
    
      /* new construction point not suitable --> do not add */
      _unur_warning(gen->genid,UNUR_ERR_ADAPT,"Cannot split segment at given point.");
#if UNUR_DEBUG & UNUR_DB_INFO
      /* write info into log file */
      if (gen->debug & AROU_DEBUG_SPLIT) 
	_unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
      
      /* remove from linked list */
      seg_oldl->rtp  = seg_newr->rtp;
      seg_oldl->drtp = seg_newr->drtp;
      seg_oldl->next = seg_newr->next;
      /* push seg_newr back to stack */
      _unur_arou_segment_stack_push(gen);
      /* we have to restore the old segment.
	 (this case should not happen, so it is easier not to make a 
	 backup of the old segment) */
      if ( !_unur_arou_segment_parameter(gen,seg_oldl) ) {
	_unur_error(gen->genid,UNUR_ERR_ADAPT,"Cannot restore segment. PANIK.");
	exit (-1);
      }
      return 0;
    }
  }
    
  /* update guide table */ 
  _unur_arou_make_guide_table(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug & AROU_DEBUG_SPLIT) 
    _unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
  
  /* o.k. */
  return 1;

} /* end of _unur_arou_segment_split() */

/*****************************************************************************/

static int
_unur_arou_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1  ... if successful                                               */
     /*   0  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg;
  double Acum, Aincum, Astep;
  int j;

  /* check arguments */
  COOKIE_CHECK(gen,CK_AROU_GEN,0);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN.guide) {
    int max_guide_size = (GEN.guide_factor > 0.) ? (GEN.max_segs * GEN.guide_factor) : 1;
    GEN.guide = _unur_malloc( max_guide_size * sizeof(struct unur_arou_segment*) );
  }

  /* first we need cumulated areas in segments */
  Acum = 0.;       /* area in enveloping polygon */
  Aincum = 0.;     /* area in squeeze */
  for (seg = GEN.seg; seg != NULL; seg = seg->next ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,0);
    Acum += seg->Ain + seg->Aout;
    Aincum += seg->Ain;
    seg->Acum = Acum;
  }

  /* total area below hat */
  GEN.Atotal = Acum;
  GEN.Asqueeze = Aincum;

  /* actual size of guide table */
  GEN.guide_size = (int)(GEN.n_segs * GEN.guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN.Atotal / GEN.guide_size;
  Acum=0.;
  for( j=0, seg=GEN.seg; j < GEN.guide_size; j++ ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,0);
    while( seg->Acum < Acum )
      if( seg->next != NULL )    /* skip to next segment if it exists */
        seg = seg->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_INIT,"roundoff error while making guide table!");
	break;
      }
    GEN.guide[j] = seg;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide[j] = seg;

  return 1;
} /* end of _unur_arou_make_guide_table() */

/*****************************************************************************/

static double
_unur_arou_segment_arcmean( struct unur_arou_segment *seg )
     /*----------------------------------------------------------------------*/
     /* compute "arctan mean" of two numbers expressed as v/u, u>=0          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   seg ... pointer to segment                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   mean                                                               */
     /*                                                                      */
     /* comment:                                                             */
     /*   "arctan mean" = tan(0.5*(arctan(x0)+arctan(x1)))                   */
     /*                                                                      */
     /*   a combination of arithmetical mean (for x0 and x1 close to 0)      */
     /*   and the harmonic mean (for |x0| and |x1| large).                   */
     /*----------------------------------------------------------------------*/
{
  double xl, xr;

  /* if u != 0 ... x is stored in tp (= v/u)   */
  /* else      ... x is stored in tangent dltp */
  xl = (seg->ltp[0] > 0.) ? (seg->ltp[1] / seg->ltp[0]) :
    ( (seg->dltp[0] == 0.) ? -INFINITY : (seg->dltp[1]) );

  xr = (seg->rtp[0] > 0.) ? (seg->rtp[1] / seg->rtp[0]) :
    ( (seg->drtp[0] == 0.) ? INFINITY : (seg->drtp[1]) );

  return _unur_arcmean(xl,xr);

} /* end of _unur_arou_segment_arcmean() */

/*****************************************************************************/

static struct unur_arou_segment *
_unur_arou_segment_stack_pop( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* pop free segment from stack; allocate memory block if necessary.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to segment                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);
  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);

  /* look for an unused segment */
  if( ! GEN.seg_free ) {
    /* no more unused segments. make some. */
    GEN.seg_stack = _unur_malloc( UNUR_MALLOC_SIZE * sizeof(struct unur_arou_segment) );

    /* reset counter */
    GEN.seg_free = UNUR_MALLOC_SIZE;

    /* set cookies */
    COOKIE_SET_ARRAY( GEN.seg_stack, CK_AROU_SEG, UNUR_MALLOC_SIZE);

    /* add to list of allocated blocks */
    _unur_add_mblocks( &(GEN.mblocks), GEN.seg_stack);
  }

  /* update ....                                   */
  --(GEN.seg_free);   /* pointer to free segments  */
  ++(GEN.n_segs);     /* counter for used segments */

  /* return pointer to segment */
  return (GEN.seg_stack + GEN.seg_free);

} /* end of _unur_arou_segment_stack_pop() */

/*---------------------------------------------------------------------------*/

static void
_unur_arou_segment_stack_push( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* push the last popped segment back onto the stack.                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_AROU_GEN,/*void*/);

  /* update counters and pointers */
  --(GEN.n_segs);
  ++(GEN.seg_free);
} /* end of _unur_arou_segment_stack_push() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
_unur_arou_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_AROU_GEN,/*void*/);
  CHECK_NULL(par,/*void*/);
  COOKIE_CHECK(par,CK_AROU_PAR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = ratio-of-uniforms method with enveloping polygon\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = unur_arou_sample",gen->genid);
  if (par->variant & AROU_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: center = %g",gen->genid,PAR.center);
  _unur_print_if_default(par,AROU_SET_CENTER);
  if (par->variant & AROU_VARFLAG_USECENTER)
    fprintf(log,"\n%s: use center as construction point",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of segments         = %d",gen->genid,PAR.max_segs);
  _unur_print_if_default(par,AROU_SET_MAX_SEGS);
  fprintf(log,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,AROU_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of segments: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
  _unur_print_if_default(par,AROU_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,AROU_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & AROU_SET_STP)
    for (i=0; i<PAR.n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",PAR.starting_cpoints[i]);
    }
  else
    fprintf(log," use \"equdistribution\" rule [default]");
  fprintf(log,"\n%s:\n",gen->genid);
  
  _unur_arou_debug_segments(gen);

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_arou_debug_init() */

/*****************************************************************************/

static void
_unur_arou_debug_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_AROU_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_arou_debug_free() */

/*****************************************************************************/

static void
_unur_arou_debug_segments( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of segments into logfile                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_arou_segment *seg;
  double sAin, sAout;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_AROU_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:Segments: %d\n",gen->genid,GEN.n_segs);
  if (gen->debug & AROU_DEBUG_SEGMENTS) {
    fprintf(log,"%s: Nr.\t    left touching point\t\t   intersection point\t\t tangent at left touching point\n",gen->genid);
    for (seg = GEN.seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,/*void*/); 
      fprintf(log,"%s:[%3d]: (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g,%-12.6g)\n", gen->genid, i,
	      seg->ltp[0],seg->ltp[1],
	      seg->mid[0],seg->mid[1],
	      seg->dltp[0],seg->dltp[1],seg->dltp[2]);
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,/*void*/); 
    fprintf(log,"%s:[...]: (%-12.6g,%-12.6g)\n", gen->genid,seg->ltp[0],seg->ltp[1]);
  }
  fprintf(log,"%s:\n",gen->genid);

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of enveloping polygon not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }
    
  /* print and sum areas inside and outside of squeeze */
  if (gen->debug & AROU_DEBUG_SEGMENTS) {
    fprintf(log,"%s:Areas in segments:\n",gen->genid);
    fprintf(log,"%s: Nr.\t inside squeeze\t\t   outside squeeze\t     total segment\t\tcumulated\n",gen->genid);
    sAin = sAout = 0.;
    for (seg = GEN.seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,/*void*/); 
      sAin += seg->Ain;
      sAout += seg->Aout;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      seg->Ain, seg->Ain * 100. / GEN.Atotal,
	      seg->Aout, seg->Aout * 100. / GEN.Atotal,
	      seg->Ain + seg->Aout, (seg->Ain + seg->Aout) * 100. / GEN.Atotal,
	      seg->Acum, seg->Acum * 100. / GEN.Atotal);
    }
    fprintf(log,"%s:\t----------  ---------  |  ----------  ---------  |  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-11.6g(%6.3f%%)\n",
	    gen->genid,
	    sAin, sAin * 100./GEN.Atotal,
	    sAout, sAout * 100./GEN.Atotal,
	    sAin + sAout, (sAin + sAout) * 100./GEN.Atotal);
    fprintf(log,"%s:\n",gen->genid);
  }

  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_arou_debug_segments() */

/*****************************************************************************/

static void
_unur_arou_debug_split_start( struct unur_gen *gen,
		       struct unur_arou_segment *seg, 
		       double x, double fx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting segment                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   seg ... pointer to segment                                         */
     /*   x   ... split at this point                                        */
     /*   fx  ... value of p.d.f. at x                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  char ratio[14];

  log = unur_get_stream();

  fprintf(log,"%s: split segment at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(log,"%s: old segment:\n",gen->genid);

  _unur_arou_debug_printratio(seg->ltp[0],seg->ltp[1],ratio);
  fprintf(log,"%s:   left  construction point = (%-12.6g,%-12.6g)\t x = v/u = %s\tf(x) = %-12.6g\n",
	  gen->genid, seg->ltp[0], seg->ltp[1], ratio, sqrt(seg->ltp[1]) );

  _unur_arou_debug_printratio(seg->mid[0],seg->mid[1],ratio);
  fprintf(log,"%s:   intersection point       = (%-12.6g,%-12.6g)\t x = v/u = %s\n",
	  gen->genid, seg->mid[0], seg->mid[1], ratio);

  _unur_arou_debug_printratio(seg->rtp[0],seg->rtp[1],ratio);
  fprintf(log,"%s:   right construction point = (%-12.6g,%-12.6g)\t x = v/u = %s\tf(x) = %-12.6g\n",
	  gen->genid, seg->rtp[0], seg->rtp[1], ratio, sqrt(seg->rtp[1]) );

  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Ain, seg->Ain * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Aout, seg->Aout * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  (seg->Ain + seg->Aout), (seg->Ain +seg->Aout) * 100./GEN.Atotal);

  fflush(log);

} /* end of _unur_arou_debug_split_start() */

/*****************************************************************************/

static void
_unur_arou_debug_split_stop( struct unur_gen *gen, 
		      struct unur_arou_segment *seg_left,
		      struct unur_arou_segment *seg_right )
     /*----------------------------------------------------------------------*/
     /* write info about new splitted segments                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv_left  ... pointer to new left hand segment                      */
     /*   iv_right ... pointer to new right hand segment                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  char ratio[14];

  log = unur_get_stream();

  fprintf(log,"%s: new segments:\n",gen->genid);

  _unur_arou_debug_printratio(seg_left->ltp[0],seg_left->ltp[1],ratio);
  fprintf(log,"%s:   left  construction point  = (%-12.6g,%-12.6g)\t x = v/u = %s\tf(x) = %-12.6g\n",
	  gen->genid, seg_left->ltp[0], seg_left->ltp[1], ratio, sqrt(seg_left->ltp[1]) );

  _unur_arou_debug_printratio(seg_left->mid[0],seg_left->mid[1],ratio);
  fprintf(log,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %s\n",
	  gen->genid, seg_left->mid[0], seg_left->mid[1], ratio );

  _unur_arou_debug_printratio(seg_left->rtp[0],seg_left->rtp[1],ratio);
  fprintf(log,"%s:   middle construction point = (%-12.6g,%-12.6g)\t x = v/u = %s\tf(x) = %-12.6g\n",
	  gen->genid, seg_left->rtp[0], seg_left->rtp[1], ratio, sqrt(seg_left->rtp[1]) );

  _unur_arou_debug_printratio(seg_right->mid[0],seg_right->mid[1],ratio);
  fprintf(log,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %s\n",
	  gen->genid, seg_right->mid[0], seg_right->mid[1], ratio );

  _unur_arou_debug_printratio(seg_right->rtp[0],seg_right->rtp[1],ratio);
  fprintf(log,"%s:   right construction point  = (%-12.6g,%-12.6g)\t x = v/u = %s\tf(x) = %-12.6g\n",
	  gen->genid, seg_right->rtp[0], seg_right->rtp[1], ratio, sqrt(seg_right->rtp[1]) );

  fprintf(log,"%s: left segment:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg_left->Ain, seg_left->Ain * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg_left->Aout, seg_left->Aout * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  (seg_left->Ain + seg_left->Aout), (seg_left->Ain +seg_left->Aout) * 100./GEN.Atotal);

  fprintf(log,"%s: right segment:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg_right->Ain, seg_right->Ain * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg_right->Aout, seg_right->Aout * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  (seg_right->Ain + seg_right->Aout), (seg_right->Ain +seg_right->Aout) * 100./GEN.Atotal);

  fprintf(log,"%s: total areas:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_arou_debug_split_stop() */

/*****************************************************************************/

static void
_unur_arou_debug_printratio( double v, double u, char *string )
     /*----------------------------------------------------------------------*/
     /* evaluate ratio v/u and write result on string.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   v      ... numerator                                               */
     /*   u      ... denominator                                             */
     /*   string ... character string of length 14                           */
     /*                                                                      */
     /* comment:                                                             */
     /*   necessary to avoid floating point exception, when u = 0.           */
     /*----------------------------------------------------------------------*/
{
  if (u!=0.)
    sprintf(string,"%-12.6g",v/u);   /** TODO: possible overflow ?? **/
  else    /* u == 0 */
    if (v==0.)
      sprintf(string,"NaN         ");
    else {
      if (v>0.)
	sprintf(string,"Inf         ");
      else
	sprintf(string,"-Inf        ");
    }

} /* end of _unur_arou_debug_printratio() */

/*****************************************************************************/
#endif

