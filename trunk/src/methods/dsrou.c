/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dsrou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    discrete simple universal method (ratio-of-uniforms method)  *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and mode of a T_{-1/2}-concave distribution                *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      mode of the density                                                  *
 *      area below PDF                                                       *
 *   OPTIONAL:                                                               *
 *      CDF at mode                                                          *
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
 *   [1] Leydold J. (2000): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions, preprint.              *
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
 * A is convex. Then we easily can construct an enveloping rectangle.        *
 * Define                                                                    *
 *    R = {(v,u):  v_l <= v <= v_r, 0 <= u <= u_m},                          *
 *    Q = {(v,u): -v_m <= v <= v_m, 0 <= u <= u_m},                          *
 * where                                                                     *
 *    u_m = sqrt(f(mode)), v_m = (\int f dx) / u_m                           *
 *    v_l = -F(\mode) v_m, v_r = (1-F(mode)) v_m                             *
 * Then                                                                      *
 *    A subset R subset Q                                                    *
 *                                                                           *
 * Thus we can use R to generate whenever the CDF F(mode) at the mode        *
 * is known, and Q otherwise.                                                *
 * Notice, that the rection constant is 2 in the first case and 4 and the    *
 * latter.                                                                   *
 *                                                                           *
 * If F(mode) it known, it is even possible to get an universal squeeze      *
 * (see [1] for details). However its usage is only recommended when         *
 * the PDF is (very) expensive.                                              *
 *                                                                           *
 * When F(mode) is not known the mirror principle can be used. i.e., make    *
 * an enveloping rectangle for f(x)+f(-x). It reduces the rejection constant *
 * to 2 * sqrt(2) at the expense of more evaluations of the PDF.             *
 * Its usage is only recommended when the generation time for the underlying *
 * uniform prng is extreme large.                                            *
 *                                                                           *
 * Distributions with a convex set A are characterized by the following      *
 * theorem that shows a connection to transformed density rejection TDR.     *
 *                                                                           *
 * THEOREM:                                                                  *
 *    A is convex if and only if g is T-concave with transformation          *
 *    T(x) = -1/sqrt(x), i.e., -1/sqrt(g(x)) is a concave function.          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define DSROU_VARFLAG_VERIFY   0x002u  /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DSROU_DEBUG_REINIT    0x00000010u  /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DSROU_SET_CDFBMODE    0x001u   /* CDF at mode-1 is known             */

/*---------------------------------------------------------------------------*/

#define GENTYPE "DSROU"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dsrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute universal bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dsrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_sample( struct unur_gen *gen );
static int _unur_dsrou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_dsrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_dsrou_debug_init( struct unur_gen *gen, int is_reinit );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       par->data.dsrou       /* data for parameter object         */
#define GEN       gen->data.dsrou       /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */     

#define PMF(x)    _unur_discr_PMF((x),&(gen->distr))  /* call to PMF         */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2     (M_SQRT2)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dsrou_new( struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }

  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (!unur_distr_discr_upd_mode(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return NULL; 
    }
  }

  if (!(distr->set & UNUR_DISTR_SET_PMFSUM))
    if (!unur_distr_discr_upd_pmfsum(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return NULL; 
    }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_DSROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.Fbmode    = -1.;                /* CDF at mode-1 (unknown yet)         */

  par->method   = UNUR_METH_DSROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dsrou_init;

  return par;

} /* end of unur_dsrou_new() */

/*****************************************************************************/

int 
unur_dsrou_set_cdfbeforemode( struct unur_par *par, double Fbmode )
     /*----------------------------------------------------------------------*/
     /* set value of cdf at mode-1                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   Fbmode ... cdf at mode-1                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );
 
  /* check input */
  _unur_check_par_object( par,DSROU );

  /* check new parameter for generator */
  if (Fbmode < 0. || Fbmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return 0;
   }
 
  /* store date */
  PAR.Fbmode = Fbmode;

  /* changelog */
  par->set |= DSROU_SET_CDFBMODE;

  return 1;
 
} /* end of unur_dsrou_set_cdfbeforemode() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par,DSROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | DSROU_VARFLAG_VERIFY) : (par->variant & (~DSROU_VARFLAG_VERIFY));

  /* o.k. */
  return 1;
} /* end of unur_dsrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
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
  CHECK_NULL( gen,0 );
  
  /* check input */
  _unur_check_gen_object( gen,DSROU );
 
  if (verify) {
    /* turn verify mode on */
    gen->variant |= DSROU_VARFLAG_VERIFY;
    SAMPLE = _unur_dsrou_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~DSROU_VARFLAG_VERIFY;
    SAMPLE = _unur_dsrou_sample;
  }

  /* o.k. */
  return 1;
 
} /* end of unur_dsrou_chg_verify() */

/*****************************************************************************/

int
unur_dsrou_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );
  
  /* set new parameters in distribution object */
  return unur_distr_discr_set_pmfparams(&(gen->distr),params,n_params);

} /* end of unur_dsrou_chg_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_mode( struct unur_gen *gen, int mode )
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_dsrou_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );

  return unur_distr_discr_upd_mode( &(gen->distr) );
} /* end of unur_dsrou_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_cdfbeforemode( struct unur_gen *gen, double Fbmode )
     /*----------------------------------------------------------------------*/
     /* change value of cdf at mode-1                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   Fbmode ... cdf at mode-1                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );

  /* check new parameter for generator */
  if (Fbmode < 0. || Fbmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode-1)");
    return 0;
  }
  
  /* copy parameters */
  GEN.Fbmode = Fbmode;

  /* changelog */
  gen->set |= DSROU_SET_CDFBMODE;

  /* o.k. */
  return 1;
} /* end of unur_dsrou_chg_cdfbeforemode() */

/*---------------------------------------------------------------------------*/

int 
unur_dsrou_chg_domain( struct unur_gen *gen, int left, int right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }
  
  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;
  
  /* changelog */
  gen->distr.set &= ~(UNUR_DISTR_SET_STDDOMAIN | UNUR_DISTR_SET_MASK_DERIVED );
  gen->distr.set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_dsrou_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_pmfsum( struct unur_gen *gen, double sum )
     /*----------------------------------------------------------------------*/
     /* change sum over PMF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   sum   ... sum                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );
  
  /* check new parameter for generator */
  if (sum <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"sum <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.sum = sum;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_dsrou_chg_pmfsum() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_upd_pmfsum( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute sum over PMF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );

  return unur_distr_discr_upd_pmfsum( &(gen->distr) );
} /* end of unur_dsrou_upd_pmfsum() */

/*****************************************************************************/

struct unur_gen *
_unur_dsrou_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
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
  if ( par->method != UNUR_METH_DSROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dsrou_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute universal bounding rectangle */
  if (! _unur_dsrou_rectangle( gen ) ) {
    free(par); _unur_dsrou_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_dsrou_debug_init(gen, FALSE);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_dsrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dsrou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{ 
  double pm, pbm;               /* PMF at mode and mode-1                    */

  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_DSROU_GEN, 0 );

  /* compute PMF at mode and mode-1 */
  pm = PMF(DISTR.mode);
  pbm = PMF(DISTR.mode-1);
  /* pm and pbm must be positive */
  if (pm <= 0. || pbm <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PMF(mode) <= 0.");
    return 0;
  }
  if (pm >= INT_MAX || pbm >= INT_MAX) { 
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return 0;
  }

  /* heights of rectangles */
  GEN.ul = sqrt(pbm);
  GEN.ur = sqrt(pm);

  /* areas of rectangle */
  if (gen->set & DSROU_SET_CDFBMODE) {
    /* cdf at mode-1 known */
    GEN.al = -GEN.Fbmode * DISTR.sum;
    GEN.ar = DISTR.sum + GEN.al;
  }
  else if (GEN.ul == 0.) {
    /* PMF monotonically decreasing */
    GEN.al = 0.;
    GEN.ar = DISTR.sum;
  }
  else {
    GEN.al = -(DISTR.sum - pm/DISTR.sum);
    GEN.ar = DISTR.sum;
  }    

  /* o.k. */
  return 1;
  
} /* end of _unur_dsrou_rectangle() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_dsrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DSROU_GEN);

  /* copy distribution object into generator object */
  _unur_distr_discr_copy( &(gen->distr), par->distr );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & DSROU_VARFLAG_VERIFY) ? _unur_dsrou_sample_check : _unur_dsrou_sample;
  gen->destroy = _unur_dsrou_free;

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or CDF at mode");
    DISTR.mode = max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = min(DISTR.mode,DISTR.BD_RIGHT);
  }

  /* copy some parameters into generator object */
  GEN.Fbmode = PAR.Fbmode;          /* CDF at mode                           */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */

  /* initialize parameters */
  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_dsrou_create() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DSROU );

  /* compute universal bounding rectangle */
  result = _unur_dsrou_rectangle( gen );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & DSROU_DEBUG_REINIT)
    if (gen->debug) _unur_dsrou_debug_init(gen,TRUE);
#endif

  return result;
} /* end of unur_dsrou_reinit() */

/*****************************************************************************/

int
_unur_dsrou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   int (sample from random variate)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  double U,V;
  int I;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DSROU_GEN,0.);

  while (1) {
    /* generate point uniformly in union of rectangles */
    V = GEN.al + _unur_call_urng(gen->urng) * (GEN.ar - GEN.al);
    V = GEN.al + V * (GEN.ar-GEN.al);
    V /= (V<0.) ? GEN.ul : GEN.ur;    /* if ul==0. then al==0. and thus V>=0. */

    while ( (U = _unur_call_urng(gen->urng)) == 0.);
    U *= (V<0.) ? GEN.ul : GEN.ur;

    /* ratio */
    I = (int)(floor(V/U)) + DISTR.mode;

    /* inside domain ? */
    if ( (I < DISTR.BD_LEFT) || (I > DISTR.BD_RIGHT) )
      continue;

    /* accept or reject */
    if (U*U <= PMF(I))
      return I;
  }
} /* end of _unur_dsrou_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_dsrou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   int (sample from random variate)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,pI,VI;
  int I;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DSROU_GEN,0.);

  while (1) {
    /* generate point uniformly in union of rectangles */
    V = GEN.al + _unur_call_urng(gen->urng) * (GEN.ar - GEN.al);
    V = GEN.al + V * (GEN.ar-GEN.al);
    V /= (V<0.) ? GEN.ul : GEN.ur;

    while ( (U = _unur_call_urng(gen->urng)) == 0.);
    U *= (V<0.) ? GEN.ul : GEN.ur;

    /* ratios */
    I = (int)(floor(V/U)) + DISTR.mode;

    /* inside domain ? */
    if ( (I < DISTR.BD_LEFT) || (I > DISTR.BD_RIGHT) )
      continue;

    /* values of PMF and v-coordinate of point */
    pI = PMF(I);
    VI = V/U * sqrt(pI);

    /* check hat */
    if ( ((2.+4.*DBL_EPSILON) * GEN.ur*GEN.ur < pI)    /* avoid roundoff error with FP registers */
	 || ( VI < (1.+UNUR_EPSILON) * GEN.al/GEN.ul) 
	 || ( VI > (1.+UNUR_EPSILON) * GEN.ar/GEN.ur) 
	 || (-VI < (1.+UNUR_EPSILON) * GEN.al/GEN.ul) 
	 || (-VI > (1.+UNUR_EPSILON) * GEN.ar/GEN.ur) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PMF(x) > hat(x)");
    
    /* accept or reject */
    if (U*U <= pI)
      return I;
  }
} /* end of _unur_dsrou_sample_check() */

/*****************************************************************************/

void
_unur_dsrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DSROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DSROU_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_distr_discr_clear(gen);
  _unur_free_genid(gen);

  free(gen);

} /* end of _unur_dsrou_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dsrou_debug_init( struct unur_gen *gen, int is_reinit )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   is_reinit ... if TRUE the generator has been reinitialized         */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DSROU_GEN,/*void*/);

  log = unur_get_stream();
  
  fprintf(log,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
    fprintf(log,"%s: method  = dsrou (discrete simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(log,"%s: reinit!\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
 
  _unur_distr_discr_debug( &(gen->distr), gen->genid, FALSE );
  
  fprintf(log,"%s: sampling routine = _unur_dsrou_sample",gen->genid);
  if (gen->variant & DSROU_VARFLAG_VERIFY)
    fprintf(log,"_check");
  fprintf(log,"()\n%s:\n",gen->genid);
 
  if (gen->set & DSROU_SET_CDFBMODE)
    fprintf(log,"%s: CDF(mode-1) = %g\n",gen->genid,GEN.Fbmode);
  else
    fprintf(log,"%s: CDFF(mode-1) unknown\n",gen->genid);

  fprintf(log,"%s: no (universal) squeeze\n",gen->genid);
  fprintf(log,"%s: use mirror principle\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: Rectangles:\n",gen->genid);
  fprintf(log,"%s:    left upper point  = (%g,%g)\n",gen->genid,GEN.al/GEN.ul,GEN.ul);
  fprintf(log,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN.ar/GEN.ur,GEN.ur);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_dsrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

