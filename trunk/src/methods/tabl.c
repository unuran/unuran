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
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
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

#define TABL_VARFLAG_VERIFY       0x001UL   /* flag for verifying mode           */

/* indicate how to split interval */
#define TABL_VARMASK_SPLIT        0x0f0UL /* split at        computation     convergence of hat */
#define TABL_VARFLAG_SPLIT_POINT  0x000UL /* sampled point    none            slowest          */
#define TABL_VARFLAG_SPLIT_MEAN   0x010UL /* mean point       slower          better           */
#define TABL_VARFLAG_SPLIT_ARC    0x020UL /* "arcmean"        very slow       very good for almost unbounded domain */

/* indicate if starting intervals have to be split */
#define TABL_VARMASK_STP          0xf00UL
#define TABL_VARFLAG_STP_A        0x100UL /* use equal area rule (SPLIT A in [1])   */
#define TABL_VARFLAG_STP_B        0x200UL /* use main subdivisions (SPLIT B in [1]) */

/*---------------------------------------------------------------------------*/
/* Debugging flags (do not use first 8 bits)                                 */

#define TABL_DEBUG_IV             0x010UL /* show intervals                        */
#define TABL_DEBUG_A_IV           0x020UL /* show intervals after split A, before split B */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TABL_SET_GUIDEFACTOR      0x01UL
#define TABL_SET_VARIANT          0x02UL
#define TABL_SET_SLOPES           0x04UL
#define TABL_SET_AREAFRACTION     0x08UL
#define TABL_SET_MAX_IVS          0x10UL
#define TABL_SET_MAX_SQHRATIO     0x20UL
#define TABL_SET_N_STP            0x40UL
#define TABL_SET_BOUNDARY         0x80UL

/*---------------------------------------------------------------------------*/

#define GENTYPE "TABL"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
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

static int 
_unur_tabl_split_b_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* split starting intervals according to [1]                                 */
/* SPLIT B, but instead of the iteration in [1] use "arcmean".               */
/*---------------------------------------------------------------------------*/

static int
_unur_tabl_split_interval( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			   double x, double fx, unsigned int split_mode );
/*---------------------------------------------------------------------------*/
/* split interval (replace old one by two new ones in same place)            */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *_unur_tabl_iv_stack_pop( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* pop an interval from the stack of free intervals.                         */
/*---------------------------------------------------------------------------*/

/*  static void _unur_tabl_iv_stack_push( struct unur_gen *gen ); */
/*---------------------------------------------------------------------------*/
/* push the last popped interval back onto the stack.                         */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_intervals( struct unur_gen *gen, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals.                                                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR   distr->data.cont
#define PAR     par->data.tabl
#define GEN     gen->data.tabl
#define SAMPLE  gen->sample.cont

#define PDF(x) ((*(GEN.pdf))((x),GEN.pdf_param,GEN.n_pdf_param))

/*---------------------------------------------------------------------------*/

#define min(x,y)   (((x)<(y)) ? (x) : (y))
#define max(x,y)   (((x)>(y)) ? (x) : (y))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_tabl_new( struct unur_distr *distr )
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

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_TABL_PAR);

  /* copy input */
  par->distr        = distr;     /* pointer to distribution object           */

  /* set default values */
  PAR.slopes        = NULL;      /* pointer to slopes of p.d.f.              */
  PAR.n_slopes      = 0;         /* number of slopes                         */

  PAR.n_starting_cpoints = 30;   /* number of starting points                */
  PAR.area_fract    = 0.25;      /* parameter for equal area rule (default from [1] ) */

  PAR.max_ivs       = 100;       /* maximum number of intervals              */
  PAR.max_ratio     = 0.95;      /* bound for ratio  Atotal / Asqueeze       */

  PAR.guide_factor  = 1.; /* guide table has same size as array of intervals */

  par->method       = UNUR_METH_TABL;         /* indicate method             */
  par->variant      = (TABL_VARFLAG_SPLIT_ARC   |     /* variant: split at arc_mean  */
		       TABL_VARFLAG_STP_A |     /* run SPLIT A on slopes       */
		       TABL_VARFLAG_STP_B  );   /* run SPLIT B on slopes       */


  par->set          = 0UL;       /* inidicate default parameters             */    
  par->urng         = unur_get_default_urng(); /* use default urng           */

  _unur_set_debugflag_default(par); /* set default debugging flags           */

  /* we use the domain to get the boundaries for our computation limits      */
  /* (but only if the domain is bounded!)                                    */
  if ( (distr->set & UNUR_DISTR_SET_DOMAIN)
       && DISTR.domain[0] > -INFINITY
       && DISTR.domain[1] <  INFINITY ) {
    PAR.bleft = DISTR.domain[0];
    PAR.bright = DISTR.domain[1];
    par->set |= TABL_SET_BOUNDARY;
  }
  else {               /* the default */
    PAR.bleft  = 0.;  /* left boundary of domain (no useful default) */
    PAR.bright = 0.;  /* right boundary of domain (left = right --> cannot make hat) */
  }

  /* routine for starting generator */
  par->init = unur_tabl_init;

  return par;

} /* end of unur_tabl_new() */

/*****************************************************************************/

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"number of starting points < 0");
    return 0;
  }

  /* store date */
  PAR.n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TABL_SET_N_STP;

  return 1;

} /* end of unur_tabl_set_nstp() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"ratio Atotal / Asqueeze not in [0,1]");
    return 0;
  }

  /* store date */
  PAR.max_ratio = max_ratio;

  /* changelog */
  par->set |= TABL_SET_MAX_SQHRATIO;

  return 1;

} /* end of unur_tabl_set_max_sqhratio() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"maximum number of intervals < 1");
    return 0;
  }

  /* store date */
  PAR.max_ivs = max_ivs;

  /* changelog */
  par->set |= TABL_SET_MAX_IVS;

  return 1;

} /* end of unur_tabl_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_areafraction( struct unur_par *par, double fraction )
     /*----------------------------------------------------------------------*/
     /* set parameter for equal area rule                                    */
     /* (each bar has size fraction * area(pdf))                             */           
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   fraction  ... fraction of area for bar                             */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if (fraction < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"area factor < 0");
    return 0;
  }

  /* store date */
  PAR.area_fract = fraction;

  /* changelog */
  par->set |= TABL_SET_AREAFRACTION;

  return 1;

} /* end of unur_tabl_set_areafraction() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_slopes( struct unur_par *par, double *slopes, int n_slopes )
     /*----------------------------------------------------------------------*/
     /* set slopes of p.d.f.                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   slopes   ... pointer to list of slopes                             */
     /*   n_slopes ... number of slopes                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a]                        */
     /*   such that pdf(a) >= pdf(b).                                        */
     /*   slopes must be decreasing, non-overlapping and sorted              */
     /*----------------------------------------------------------------------*/
{
  int i;
  double al, bl;

  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if( n_slopes <= 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"number of slopes <= 0");
    return 0;
  }

  /* check slopes */
  al = slopes[0];
  bl = slopes[1];
  for( i=1; i<n_slopes; i++ ) {
    /* we do not check here if f(a) >= f(b), since we make no calculations heres */
    if( al > slopes[2*i] || bl > slopes[2*i+1] ) {
      _unur_warning(GENTYPE,UNUR_ERR_SET,"slopes overlapping or not in ascending order");
      return 0;
    }
    al = slopes[2*i];
    bl = slopes[2*i+1];
  }

  /* store date */
  PAR.slopes = slopes;
  PAR.n_slopes = n_slopes;

  /* changelog */
  par->set |= TABL_SET_SLOPES;

  return 1;

} /* end of unur_tabl_set_slopes() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   new boundary points must not be +/- INFINITY                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"domain, left >= right");
     /*                                                                      */
    return 0;
  }
  if (left <= -INFINITY || right >= INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"domain, +/- INFINITY not allowed");
    return 0;
  }

  /* store date */
  PAR.bleft = left;
  PAR.bright = right;

  /* changelog */
  par->set |= TABL_SET_BOUNDARY;

  return 1;

} /* end of unur_tabl_set_boundary() */

/*---------------------------------------------------------------------------*/

int 
unur_tabl_set_variant( struct unur_par *par, unsigned int variant )
     /*----------------------------------------------------------------------*/
     /* set variant of method                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   variant ... indicator for variant                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* store date */
  par->variant = variant;

  /* changelog */
  par->set |= TABL_SET_VARIANT;

  return 1;

} /* end if unur_tabl_set_variant() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check input */
  _unur_check_par_object( TABL );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_SET,"relative table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= TABL_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_tabl_set_guidefactor() */

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
  _unur_check_par_object( TABL );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TABL_VARFLAG_VERIFY) : (par->variant & (~TABL_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_tabl_set_verify() */

/*****************************************************************************/

struct unur_gen *
unur_tabl_init( struct unur_par *par )
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
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_TABL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tabl_create(par);
  if (!gen) { free(par); return NULL; }

  /* get slopes for starting generator */
  if (!_unur_tabl_get_starting_intervals(par,gen)) {
    _unur_error(gen->genid,UNUR_ERR_INIT,"Cannot make hat function.");
    free(par); unur_tabl_free(gen);
    return NULL;
  }

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_tabl_debug_init(par,gen);
  if (gen->debug & TABL_DEBUG_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);
#endif

  /* split according to [1], run SPLIT B */
  if (par->variant & TABL_VARFLAG_STP_B)
    while ( (GEN.n_ivs < PAR.n_starting_cpoints)
	    && (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) )
      if (!_unur_tabl_split_b_starting_intervals(par,gen)) {
	free(par); unur_tabl_free(gen);
	return NULL;
      }
  
  /* we have to update the maximal number of intervals,
     if the user wants more starting points. */
  if (GEN.n_ivs > GEN.max_ivs) {
    /** TODO: do not allow too many intervals ?? **/
    _unur_warning(gen->genid,UNUR_ERR_INIT,"maximal number of intervals too small. increase.");
    GEN.max_ivs = GEN.n_ivs;
  }

  /* make initial guide table */
  _unur_tabl_make_guide_table(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug)
    _unur_tabl_debug_intervals(gen,TRUE);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_tabl_init() */

/*****************************************************************************/

double
unur_tabl_sample_adaptive( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    do {
      /* since we do not upgrade the guide table every time we 
	 split an interval, there exists "black intervals".
	 The total area below the hat decreases, but the values 
	 of iv->Acum cannot be upgraded every time.
	 So there arise holes, i.e., these "black intervals",
	 that represent area above the hat function (cut off by
	 the splitting process).
	 if we hit such a hole, we have to reject the generated 
	 interval, and try again. 
      */

      /* sample from U(0,1) */
      u = _unur_call_urng(gen);
    
      /* look up in guide table and search for interval */
      iv =  GEN.guide[(int) (u * GEN.guide_size)];
      COOKIE_CHECK(iv,CK_TABL_IV,0.);
      u *= GEN.Atotal;
      while (iv->Acum < u) {
	if (iv->next == NULL)
	  /* we have hit an imaginary "black interval" at the end of the list. try again. */
	  break;
	iv = iv->next;
	COOKIE_CHECK(iv,CK_TABL_IV,0.);
      }
      
      /* reuse of uniform random number */
      u = iv->Acum - u;
    
    } while (u > iv->Ahat); /* check whether we have hit a "black interval" */

    /* generation w.r.t. squeeze should be inversion */
    if (iv->slope>0)
      u = iv->Ahat - u;

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
      /** TODO: possible overflow/underflow ?? **/
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (gen->variant & TABL_VARMASK_SPLIT) );
  	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample_adaptive() */

/*****************************************************************************/

double
unur_tabl_sample( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for interval */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,0.);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u < iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
      /** TODO: possible overflow/underflow ?? **/
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (gen->variant & TABL_VARMASK_SPLIT) );
	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample() */

/*****************************************************************************/

double
unur_tabl_sample_check( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for interval */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,0.);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      x = iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze;
      /** TODO: possible overflow/underflow ?? **/
      /* test whether p.d.f. is monotone */
      fx = PDF(x);
      if (fx > iv->fmax)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf > hat. pdf not monotone in interval");
      if (fx < iv->fmin)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf < squeeze. pdf not monotone in interval");
      /* at last return number */
      return x;
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* test whether p.d.f. is monotone */
      if (fx > iv->fmax)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf > hat. pdf not monotone in interval");
      if (fx < iv->fmin)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf < squeeze. pdf not monotone in interval");
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (gen->variant & TABL_VARMASK_SPLIT) );
	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
  
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample_check() */

/*****************************************************************************/

void
unur_tabl_free( struct unur_gen *gen )
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
    _unur_warning(GENTYPE,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  if (gen->debug) _unur_tabl_debug_free(gen);
#endif

  /* free linked list of intervals and others */
  _unur_free_mblocks(GEN.mblocks);

  /* free other memory */
  _unur_free_genid(gen);
  free(GEN.guide);
  free(gen);

} /* end of unur_tabl_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

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
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TABL_GEN);

  /* set generator identifier */
  _unur_set_genid(gen,GENTYPE);

  /* copy pointer to distribution object */
  /* (we do not copy the entire object)  */
  gen->distr = par->distr;

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & TABL_VARFLAG_VERIFY) ? unur_tabl_sample_check : unur_tabl_sample;
  gen->destroy = unur_tabl_free;

  /* set all pointers to NULL */
  GEN.pdf_param   = NULL;
  GEN.n_pdf_param = 0;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.iv          = NULL;
  GEN.n_ivs       = 0;
  GEN.iv_stack    = NULL;
  GEN.iv_free     = 0;
  GEN.mblocks     = NULL;

  /* copy some parameters into generator object */
  GEN.pdf = gen->DISTR.pdf;           /* p.d.f. of distribution              */

  GEN.pdf_param   = gen->DISTR.params;
  GEN.n_pdf_param = gen->DISTR.n_params;

  GEN.bleft       = PAR.bleft;         /* left boundary of domain            */
  GEN.bright      = PAR.bright;        /* right boundary of domain           */

  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN.max_ivs   = PAR.max_ivs;         /* maximum number of intervals        */
  GEN.max_ratio = PAR.max_ratio;       /* bound for ratio  Atotal / Asqueeze */

  gen->method = par->method;           /* indicates method                   */
  gen->variant = par->variant;         /* indicates variant                  */

  _unur_copy_urng_pointer(par,gen);    /* pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);       /* copy debugging flags into generator object */

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_tabl_create() */

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
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   pdf(a) >= pdf(b)                                                   */
     /*----------------------------------------------------------------------*/
{

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* we have two cases: 
     (1) we are given slopes --> check these, compute domain if necessary
     (2) we are given domain and mode --> compute slopes */

  if (PAR.n_slopes > 0 )
    /* slopes are given */
    return _unur_tabl_get_starting_intervals_from_slopes(par,gen);

  if ( (par->set & TABL_SET_BOUNDARY) && (par->distr->set & UNUR_DISTR_SET_MODE) )
    /* no slopes given. need domain and mode */
    /* compute slopes */
    return _unur_tabl_get_starting_intervals_from_mode(par,gen);

  /* else */
  _unur_error(gen->genid,UNUR_ERR_INIT,"number of slopes <= 0, domain or mode not given.");
  return 0;

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
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   pdf(a) >= pdf(b)                                                   */
     /*----------------------------------------------------------------------*/
{
  /** TODO: check for slopes out of support !! **/

  struct unur_tabl_interval *iv;
  int i;

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* init counter of intervals */
  GEN.n_ivs = 0;
  iv = NULL;

  /* compute initial intervals */
  for ( i=0; i < 2*PAR.n_slopes; i+=2 ) {
    /* get a new interval and link into list */
    if (i==0)
      iv = GEN.iv = _unur_tabl_iv_stack_pop(gen);    /* the first interval */
    else
      iv = iv->next = _unur_tabl_iv_stack_pop(gen);  /* all the other intervals */
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    /* max and min of p.d.f. in interval */
    iv->xmax = PAR.slopes[i];      
    iv->fmax = PDF(iv->xmax);
    iv->xmin = PAR.slopes[i+1];    
    iv->fmin = PDF(iv->xmin);

    /* check slopes */
    if (iv->fmax < iv->fmin) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"slopes non-decreasing");
      return 0;
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
    if (!(par->set & TABL_SET_BOUNDARY)) {
      if (iv->slope > 0) {
	GEN.bleft = min(GEN.bleft,iv->xmin);
	GEN.bright = max(GEN.bright,iv->xmax);
      }
      else {
	GEN.bleft = min(GEN.bleft,iv->xmax);
	GEN.bright = max(GEN.bright,iv->xmin);
      }
    }

    /* split interval following [1], split A */
    if (par->variant & TABL_VARFLAG_STP_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return 0;
    }
  }

  /* terminate list */
  iv->next = NULL;

  /* o.k. */
  return 1;

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
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*----------------------------------------------------------------------*/
{
  /** TODO: check for slopes out of support !! **/

  struct unur_tabl_interval *iv;
  
  double mode = par->DISTR.mode;   /* (exact!) location of mode of p.d.f. */

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* init linked list of intervals */
  GEN.n_ivs = 0;

  /* compute initial intervals */
  while (1) {
    /* the first interval */
    iv = GEN.iv = _unur_tabl_iv_stack_pop(gen);
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    if (mode <= GEN.bleft) {
      /* only one ascending interval <a,b> = [a,b] */
      iv->xmax = mode;
      iv->xmin = GEN.bright;
      break;
    }

    if (mode >= GEN.bright) {
      /* only one descending interval <a,b> = [b,a] */
      iv->xmax = mode;
      iv->xmin = GEN.bleft;
      break;
    }

    /* one descending and one ascending interval */
    iv->xmax = mode;
    iv->xmin = GEN.bleft;

    /* the second interval */
    iv = iv->next = _unur_tabl_iv_stack_pop(gen);  /* all the other intervals */
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    iv->xmax = mode;
    iv->xmin = GEN.bright;
    break;
  }

  /* terminate list */
  iv->next = NULL;

  /* compute parameters */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    /* max and min of p.d.f. in interval */
    iv->fmax = PDF(iv->xmax);
    iv->fmin = PDF(iv->xmin);

    /* area of slope and sign of slope (increasing/decreasing) */
    iv->slope = (iv->xmax > iv->xmin) ? 1 : -1;
    iv->Ahat = iv->slope * (iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = iv->slope * (iv->xmax - iv->xmin) * iv->fmin;
    /** TODO: possible overflow/underflow ?? **/
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* split interval following [1], split A */
    if (par->variant & TABL_VARFLAG_STP_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return 0;
    }

  }

  /* o.k. */
  return 1;

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
  double bar_area, x;
  
  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  COOKIE_CHECK(iv_slope,CK_TABL_IV,NULL);
  if (iv_slope->slope != 1 && iv_slope->slope != -1 ) {
    _unur_warning( gen->genid, UNUR_ERR_INIT, "invalid slope.");
    return iv_slope;
  }

  iv = iv_slope;        /* pointer to actual interval */
  iv_last = iv_slope;   /* pointer to last interval in list */
  /* (maximal) area of bar (= hat in one interval) */
  bar_area = par->DISTR.area * PAR.area_fract;

  switch (iv->slope) {
  case +1:
    /* move from right to left */
    while (iv->Ahat > bar_area) {
      x = iv->xmax - bar_area / iv->fmax;
      switch (_unur_tabl_split_interval( gen, iv, x, PDF(x), TABL_VARFLAG_SPLIT_POINT )) {
      case 1:  /* splitting succesful */
	if (iv_last == iv_slope)
	  iv_last = iv->next;
	break;
      case -1: /* interval chopped */
	break; /* nothing to do */
      case 0:  /* error (slope not monotonically increasing) */
	return NULL;
      }
    }
    break;

  case -1:
    /* move from left to right */
    while (iv->Ahat > bar_area) {
      x = iv->xmax + bar_area / iv->fmax;
      switch (_unur_tabl_split_interval( gen, iv, x, PDF(x), TABL_VARFLAG_SPLIT_POINT )) {
      case 1:  /* splitting succesful */
	iv = iv->next; break;
      case -1: /* interval chopped */
	break; /* nothing to do */
      case 0:  /* error (slope not monotonically decreasing) */
	return NULL;
      }
    }
    break;
  }

  /* pointer to last interval */
  return ((iv->slope < 0) ? iv : iv_last);

} /* end of _unur_tabl_split_a_starting_intervals() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_split_b_starting_intervals( struct unur_par *par, 
				       struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* split starting intervals according to [1]                            */
     /* SPLIT B, but instead of the iteration in [1] use "arcmean".          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter list                                  */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Amean;  /* mean area between hat and squeeze in slope */
  
  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* compute mean area between squeeze and hat for each interval */
  Amean = 0.;
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    Amean += iv->Ahat - iv->Asqueeze;
  }
  Amean /= GEN.n_ivs;

  /* now split intervals */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    if ((iv->Ahat - iv->Asqueeze) >= Amean) {
      /* new point instead of the interation of [1] we use "arcmean" */
      switch (_unur_tabl_split_interval( gen, iv, 0., 0., TABL_VARFLAG_SPLIT_ARC )) {
      case 1:  /* splitting succesful */
      case -1: /* interval chopped */
	break; /* nothing to do */
      case 0:  /* error (slope not monotonically decreasing) */
	return 0;
      }
      if (GEN.n_ivs >= PAR.n_starting_cpoints)
	/* no more intervals, yet */
	break;
    }
  }

  /* o.k. */
  return 1;

} /* end of _unur_tabl_split_b_starting_intervals() */

/*****************************************************************************/

/*  static struct unur_tabl_interval * */
static int
_unur_tabl_split_interval( struct unur_gen *gen,
				struct unur_tabl_interval *iv_old, 
				double x, double fx, 
				unsigned int split_mode )
     /*----------------------------------------------------------------------*/
     /* split interval (replace old one by two new ones in same place)       */
     /* new interval is inserted immedately after old one.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval that has to be split                   */
     /*   x   ... splitting point                                            */
     /*   fx  ... value of p.d.f. at splitting point                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... splitting successful                                         */
     /*  -1 ... interval chopped off domain (not part of support of p.d.f.)  */
     /*   0 ... error: p.d.f. not monoton in interval                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv_new;

  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,0);
  COOKIE_CHECK(iv_old,CK_TABL_IV,0);

  /* There are three possibilities for the splitting point:
     (1) use x and avoid computation of pdf(x). 
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
    /* this should not happen */
    _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"Invalid variant, use default");
    break;
  }

  /* check if the new interval is completely outside the support of p.d.f. */
  if (fx <= 0.) {
    /* check montonicity */
    if (iv_old->fmin > 0.) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"p.d.f. not monotone in slope");
      return 0;
    }

    /* chop off part out of support */
    iv_old->xmin = x;

    /* compute new area in interval */
    /** TODO: possible overflow/underflow ?? **/
    iv_old->Ahat = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmax;
    /* iv_old->Asqueeze remains 0 */

    /* interval chopped but not split */
    return -1;
  }

  /* we need a new interval */
  iv_new = _unur_tabl_iv_stack_pop(gen);
  COOKIE_CHECK(iv_new,CK_TABL_IV,0);

  /* iv_new has the same slope as iv_old */
  iv_new->slope = iv_old->slope;

  /* we have to distinguish between two cases:
     pdf is increasing (slope = +1) or
     pdf is decreasing (slope = -1). */

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
    /* this should not happen */
    _unur_error(gen->genid,UNUR_ERR_SAMPLE,"Invalid slope. Cannot split interval.");
    return 0;
  }

  /* compute the areas in both intervals */
  /** TODO: possible overflow/underflow ?? **/
  iv_new->Ahat     = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmax;
  iv_new->Asqueeze = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmin;
  iv_old->Ahat     = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmax;
  iv_old->Asqueeze = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmin;

  /* insert iv_new into linked list of intervals.
     iv_old is stored on the left hand side of iv_new. */
  iv_new->next = iv_old->next;
  iv_old->next = iv_new;

  /* splitting successful */
  return 1;

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
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN.guide) {
    int max_guide_size = (GEN.guide_factor > 0.) ? (GEN.max_ivs * GEN.guide_factor) : 1;
    GEN.guide = _unur_malloc( max_guide_size * sizeof(struct unur_tabl_interval*) );
  }

  /* first we need the cumulated areas of rectangles */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
    
  /* total area below hat */
  GEN.Atotal = Acum;
  GEN.Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN.guide_size = GEN.n_ivs;

  /* make table (use variant 2; see dis.c) */
  Astep = GEN.Atotal / GEN.guide_size;
  Acum=0.;
  for( j=0, iv=GEN.iv; j < GEN.guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    while( iv->Acum < Acum )
      if( iv->next != NULL )    /* skip to next segment if it exists */
        iv = iv->next;
      else {
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
} /* end of _unur_tabl_make_guide_table() */

/*****************************************************************************/

static struct unur_tabl_interval *
_unur_tabl_iv_stack_pop( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);

  /* look for an unused segment */
  if( ! GEN.iv_free ) {
    /* no more unused segments. make some. */
    GEN.iv_stack = _unur_malloc( UNUR_MALLOC_SIZE * sizeof(struct unur_tabl_interval) );

    /* reset counter */
    GEN.iv_free = UNUR_MALLOC_SIZE;

    /* set cookies */
    COOKIE_SET_ARRAY( GEN.iv_stack, CK_TABL_IV, UNUR_MALLOC_SIZE);

    /* add to list of allocated blocks */
    _unur_add_mblocks( &(GEN.mblocks), GEN.iv_stack ); 
  }

  /* update ....                                   */
  --(GEN.iv_free);   /* pointer to free segments  */
  ++(GEN.n_ivs);     /* counter for used segments */

  /* return pointer to segment */
  return (GEN.iv_stack + GEN.iv_free);

} /* end of _unur_tabl_iv_stack_pop() */

/*---------------------------------------------------------------------------*/

#if THIS_IS_NOT_USED_YET & 0
static void
_unur_tabl_iv_stack_push( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* push the last popped interval back onto the stack.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  /* update counters and pointers */
  --(GEN.n_ivs);
  ++(GEN.iv_free);
} /* end of _unur_tabl_iv_stack_push() */
#endif

/*-----------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(log,"%s:\n",gen->genid);

/*---------------------------------------------------------------------------*/

static void
_unur_tabl_debug_init( struct unur_par *par, struct unur_gen *gen )
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

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  empty_line();

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = unur_tabl_sample",gen->genid);
  if (par->variant & TABL_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  empty_line();

  fprintf(log,"%s: computation interval = (%g, %g)",gen->genid,GEN.bleft,GEN.bright);
  if (!(par->set & TABL_SET_BOUNDARY))
      fprintf(log,"   (computed from given slopes)");
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: area fraction for equal area rule = %g ",gen->genid,PAR.area_fract);
  _unur_print_if_default(par,TABL_SET_AREAFRACTION);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR.max_ivs);
  _unur_print_if_default(par,TABL_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,TABL_SET_MAX_SQHRATIO);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
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
    fprintf(log,"%s: slopes = %d\n",gen->genid,PAR.n_slopes);
    for (i=0; i<PAR.n_slopes; i++) {
      if ( PAR.slopes[2*i] > PAR.slopes[2*i+1] )
	fprintf(log,"%s:   (+)  ",gen->genid);
      else
	fprintf(log,"%s:   (-)  ",gen->genid);
      fprintf(log,"< %#g, %#g >\n", PAR.slopes[2*i], PAR.slopes[2*i+1] );
    }
  }
  else
    fprintf(log,"%s: no slopes given. compute from domain and mode.\n",gen->genid);

  if (par->variant & TABL_VARFLAG_STP_A)
    fprintf(log,"%s: split slopes by equal area rule (SPLIT A).\n",gen->genid);
  if (par->variant & TABL_VARFLAG_STP_B)
    fprintf(log,"%s: split slopes by main subdivision rule (SPLIT B).\n",gen->genid);
  empty_line();

  fprintf(log,"%s: number of starting intervals (at least) = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,TABL_SET_N_STP);
  fprintf(log,"\n");
  empty_line();

} /* end of _unur_tabl_debug_init() */

/*****************************************************************************/

static void
_unur_tabl_debug_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  empty_line();
  _unur_tabl_debug_intervals(gen,TRUE);
  empty_line();

  fflush(log);

} /* end of _unur_tabl_debug_free() */

/*****************************************************************************/

static void
_unur_tabl_debug_intervals( struct unur_gen *gen, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tabl_interval *iv;
  double sAsqueeze;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: intervals = %d\n",gen->genid,GEN.n_ivs);
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    for (iv = GEN.iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,/*void*/);
      fprintf(log,"%s:[%3d]: (",gen->genid,i);
      switch (iv->slope) {
      case 1:  fprintf(log,"+"); break;
      case 0:  fprintf(log,"0"); break;
      case -1: fprintf(log,"-"); break;
      default:
	_unur_warning(gen->genid,UNUR_ERR_GENERIC,"invalid value for iv->slope.");
      }
      fprintf(log,")   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }

  if (!print_areas) return;

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN.iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,/*void*/); 
      sAsqueeze += iv->Asqueeze;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / GEN.Atotal,
	      iv->Ahat, iv->Ahat * 100. / GEN.Atotal, 
	      iv->Acum, iv->Acum * 100. / GEN.Atotal);
    }
    fprintf(log,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / GEN.Atotal, GEN.Atotal);
    empty_line();
  }
    
  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  empty_line();

} /* end of _unur_tabl_debug_intervals */

/*****************************************************************************/

#endif
