/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      srou.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    simple universal method (ratio-of-uniforms method)           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and mode of a T_{-1/2}-concave distribution              *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      mode of the density                                                  *
 *      area below p.d.f.                                                    *
 *   OPTIONAL:                                                               *
 *      c.d.f. at mode                                                       *
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
 * Thus we can use R to generate whenever the c.d.f. F(mode) at the mode     *
 * is known, and Q otherwise.                                                *
 * Notice, that the rection constant is 2 in the first case and 4 and the    *
 * latter.                                                                   *
 *                                                                           *
 * If F(mode) it known, it is even possible to get an universal squeeze      *
 * (see [1] for details). However its usage is only recommended when         *
 * the p.d.f. is (very) expensive.                                           *
 *                                                                           *
 * When F(mode) is not known the mirror principle can be used. i.e., make    *
 * an enveloping rectangle for f(x)+f(-x). It reduces the rejection constant *
 * to 2 * sqrt(2) at the expense of more evaluations of the p.d.f.           *
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
/* Variants: none                                                            */

#define SROU_VARFLAG_VERIFY   0x002u   /* run verify mode                    */
#define SROU_VARFLAG_SQUEEZE  0x004u   /* use universal squeeze if possible  */
#define SROU_VARFLAG_MIRROR   0x008u   /* use mirror principle               */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define SROU_SET_CDFMODE      0x001u   /* cdf at mode is known               */
#define SROU_SET_PDFMODE      0x002u   /* pdf at mode is set                 */

/*---------------------------------------------------------------------------*/

#define GENTYPE "SROU"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_srou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_srou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute universal bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_srou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_srou_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Re-initialize (existing) generator.                                       */
/*---------------------------------------------------------------------------*/

static double _unur_srou_sample( UNUR_GEN *generator );
static double _unur_srou_sample_mirror( UNUR_GEN *generator );
static double _unur_srou_sample_check( UNUR_GEN *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_srou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_srou_debug_init( struct unur_par *par, struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.srou        /* data for parameter object         */
#define GEN       gen->data.srou        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),&(gen->distr))   /* call to p.d.f.      */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2    1.4142135623731

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_srou_new( struct unur_distr *distr )
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
  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_PDFAREA)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below p.d.f.");
    return NULL; }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_SROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.Fmode     = -1.;                /* c.d.f. at mode (unknown yet)        */
  PAR.um        = -1.;                /* square root of p.d.f. at mode (unknown)*/

  par->method   = UNUR_METH_SROU;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_srou_init;

  return par;

} /* end of unur_srou_new() */

/*****************************************************************************/

int 
unur_srou_set_cdfatmode( struct unur_par *par, double Fmode )
     /*----------------------------------------------------------------------*/
     /* set value of cdf at mode                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   Fmode ... cdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cdf(mode)");
    return 0;
  }

  /* store date */
  PAR.Fmode = Fmode;

  /* changelog */
  par->set |= SROU_SET_CDFMODE;

  return 1;

} /* end of unur_srou_set_cdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_pdfatmode( UNUR_PAR *par, double fmode )
     /*----------------------------------------------------------------------*/
     /* Set pdf at mode. if set the p.d.f. at the mode is never changed.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   fmode ... pdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SROU );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"pdf(mode)");
    return 0;
  }

  /* store date (square root of fmode) */
  PAR.um = sqrt(fmode);

  /* changelog */
  par->set |= SROU_SET_PDFMODE;

  return 1;

} /* end of unur_srou_set_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par,SROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | SROU_VARFLAG_VERIFY) : (par->variant & (~SROU_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_srou_set_verify() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_usesqueeze( struct unur_par *par, int usesqueeze )
     /*----------------------------------------------------------------------*/
     /* set flag for using universal squeeze (default: off)                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   usesqueeze ... 0 = no squeeze,  !0 = use squeeze                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no squeeze is the default                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SROU );

  /* we use a bit in variant */
  par->variant = (usesqueeze) ? (par->variant | SROU_VARFLAG_SQUEEZE) : (par->variant & (~SROU_VARFLAG_SQUEEZE));

  /* o.k. */
  return 1;

} /* end of unur_srou_set_usesqueeze() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_set_usemirror( struct unur_par *par, int usemirror )
     /*----------------------------------------------------------------------*/
     /* set flag for using mirror principle (default: off)                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   usemirror ... 0 = no mirror princ.,  !0 = use mirror principle     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   do not use mirror principle is the default                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SROU );

  /* we use a bit in variant */
  par->variant = (usemirror) ? (par->variant | SROU_VARFLAG_MIRROR) : (par->variant & (~SROU_VARFLAG_MIRROR));

  /* o.k. */
  return 1;

} /* end of unur_srou_set_usemirror() */

/*****************************************************************************/

int
unur_srou_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  register int i;

  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );
  if (n_params>0) CHECK_NULL(params,0);
  
  /* check new parameter for generator */
  if (n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* copy parameters */
  DISTR.n_params = n_params;
  for (i=0; i < n_params; i++)
    DISTR.params[i] = params[i];

  /* changelog */
  /* mode and area might be wrong now! 
     but the user is responsible to change it.
     so we dont say:
     gen->distr.set &= ~(UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_PDFAREA );
     gen->set &= ~SROU_SET_CDFMODE;
  */

  /* o.k. */
  return 1;
} /* end of unur_srou_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_mode( struct unur_gen *gen, double mode )
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
  _unur_check_gen_object( gen,SROU );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_srou_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
     /*----------------------------------------------------------------------*/
     /* change value of cdf at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   Fmode ... cdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"cdf(mode)");
    return 0;
  }
  
  /* copy parameters */
  GEN.Fmode = Fmode;

  /* changelog */
  gen->set |= SROU_SET_CDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_srou_chg_cdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_pdfatmode( struct unur_gen *gen, double fmode )
     /*----------------------------------------------------------------------*/
     /* change value of pdf at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   fmode ... pdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"pdf(mode)");
    return 0;
  }

  /* store date (square root of fmode) */
  GEN.um = sqrt(fmode);

  /* changelog */
  gen->set |= SROU_SET_PDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_srou_chg_pdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_srou_chg_domain( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
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
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_srou_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_srou_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* change area below p.d.f. of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   area  ... area                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );
  
  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"area <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_srou_chg_pdfarea() */

/*****************************************************************************/

struct unur_gen *
_unur_srou_init( struct unur_par *par )
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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_SROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SROU_PAR,NULL);

  if (par->set & SROU_SET_CDFMODE)
    /* cdf at mode known -->
       thus it does not make sense to use the mirror principle */
    par->variant &= ~SROU_VARFLAG_MIRROR;
  else
    /* cdf at mode unknown -->
       thus we cannot use universal squeeze */
    par->variant &= ~SROU_VARFLAG_SQUEEZE;

  /* create a new empty generator object */
  gen = _unur_srou_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute universal bounding rectangle */
  if (! _unur_srou_rectangle( gen ) ) {
    free(par); _unur_srou_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_srou_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_srou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_srou_rectangle( struct unur_gen *gen )
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
  double vm, fm;             /* width of rectangle, pdf at mode              */

  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_SROU_GEN, 0 );

  /* compute pdf at mode (if not given by user) */
  if (!(gen->set & SROU_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"pdf(mode) <= 0.");
      return 0;
    }
    GEN.um = sqrt(fm);    /* height of rectangle */
  }

  /* width of rectangle */
  vm = DISTR.area / GEN.um;

  if (gen->set & SROU_SET_CDFMODE) {
    /* cdf at mode known */
    GEN.vl = -GEN.Fmode * vm;
    GEN.vr = vm + GEN.vl;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
  }
  else {
    /* cdf at mode unknown */
    GEN.vl = -vm;
    GEN.vr = vm;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
    /* we cannot use universal squeeze */
    gen->variant &= ~SROU_VARFLAG_SQUEEZE;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
#endif

  /* o.k. */
  return 1;

} /* end of _unur_srou_rectangle() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_srou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SROU_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_SROU_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  if (par->variant & SROU_VARFLAG_VERIFY)
    SAMPLE = _unur_srou_sample_check;
  else
    SAMPLE = (par->variant & SROU_VARFLAG_MIRROR) ? _unur_srou_sample_mirror : _unur_srou_sample;

  gen->destroy = _unur_srou_free;
  gen->reinit = _unur_srou_reinit;

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or cdf at mode");
    DISTR.mode = max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = min(DISTR.mode,DISTR.BD_RIGHT);
  }

  /* copy some parameters into generator object */
  GEN.Fmode = PAR.Fmode;            /* cdf at mode                           */
  GEN.um    = PAR.um;               /* square root of pdf at mode            */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* initialize parameters */

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_srou_create() */

/*---------------------------------------------------------------------------*/

int
_unur_srou_reinit( struct unur_gen *gen )
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
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SROU );

  /* compute universal bounding rectangle */
  return _unur_srou_rectangle( gen );

} /* end of _unur_srou_reinit() */

/*****************************************************************************/

double
_unur_srou_sample( struct unur_gen *gen )
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
  double U,V,X,x,xx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SROU_GEN,0.);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( (U = _unur_call_urng(gen)) == 0.);
    U *= GEN.um;
    V = GEN.vl + _unur_call_urng(gen) * (GEN.vr - GEN.vl);

    /* ratio */
    X = V/U;

    /* compute x */
    x = X + DISTR.mode;

    /* inside domain ? */
    if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
      continue;

    /* evaluate squeeze */
    if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	 (X >= GEN.xl) && 
	 (X <= GEN.xr ) && 
	 (U < GEN.um) ) {
      xx = V / (GEN.um - U);
      if ( (xx >= GEN.xl) && (xx <= GEN.xr ) )
	return x;
    }

    /* accept or reject */
    if (U*U <= PDF(x))
      return x;
  }

} /* end of _unur_srou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_srou_sample_mirror( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator, use mirror principle                          */
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
  double U,V,X,x,fx,uu;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SROU_GEN,0.);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( (U = _unur_call_urng(gen)) == 0.);
    U *= GEN.um * SQRT2;
    V = GEN.vl + _unur_call_urng(gen) * (GEN.vr - GEN.vl);

    /* ratio */
    X = V/U;

    /* evaluate p.d.f. */
    x = X + DISTR.mode;
    fx = PDF(x);
    uu = U * U;

    /* accept or reject */
    if (uu <= fx) {
      if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	continue; /* not inside domain */
      else
	return x;
    }

    /* try mirrored p.d.f */
    x = -X + DISTR.mode;
    if (uu <= fx + PDF(x)) {
      if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	continue; /* not inside domain */
      else
	return x;
    }
  }

} /* end of _unur_srou_sample_mirror() */

/*---------------------------------------------------------------------------*/

double
_unur_srou_sample_check( struct unur_gen *gen )
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
  double U,uu,V,X,x,fx,sfx,fnx,xfx,xfnx,xx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SROU_GEN,0.);

  if (gen->variant & SROU_VARFLAG_MIRROR) {
    /* use mirror principle */

    while (1) {
      /* generate point uniformly on rectangle */
      while ( (U = _unur_call_urng(gen)) == 0.);
      U *= GEN.um * SQRT2;
      V = GEN.vl + _unur_call_urng(gen) * (GEN.vr - GEN.vl);

      /* ratio */
      X = V/U;

      /* evaluate p.d.f. */
      fx = PDF(X + DISTR.mode);
      fnx = PDF(-X + DISTR.mode);
      uu = U * U;

      /* inside domain ? */
      if (uu <= fx) {
	x = X + DISTR.mode;
	if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	  continue; /* not inside domain */
      }
      else if (uu <= fx + fnx) {
	x = -X + DISTR.mode;
	if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	  continue; /* not inside domain */
      }

      /* check hat */
      xfx = X * sqrt(fx);
      xfnx = -X * sqrt(fnx);
      if ( ((2.+4.*DBL_EPSILON) * GEN.um*GEN.um < fx + fnx)    /* avoid roundoff error with FP registers */
	   || (xfx < (1.-DBL_EPSILON) * GEN.vl) 
	   || (xfx > (1.+DBL_EPSILON) * GEN.vr)
	   || (xfnx < (1.-DBL_EPSILON) * GEN.vl) 
	   || (xfnx > (1.+DBL_EPSILON) * GEN.vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"pdf(x) > hat(x)");
      
      /* accept or reject */
      if (uu <= fx)
	return (X + DISTR.mode);
      
      /* try mirrored p.d.f */
      x = -X + DISTR.mode;
      if (uu <= fx + PDF(x))
	return x;
    }
  }

  else { /* do not use mirror principle */

    while (1) {
      /* generate point uniformly on rectangle */
      while ( (U = _unur_call_urng(gen)) == 0.);
      U *= GEN.um;
      V = GEN.vl + _unur_call_urng(gen) * (GEN.vr - GEN.vl);

      /* ratio */
      X = V/U;

      /* compute x */
      x = X + DISTR.mode;

      /* inside domain ? */
      if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) )
	continue;

      /* evaluate p.d.f. */
      fx = PDF(x);

      /* check hat */
      sfx = sqrt(fx);
      xfx = X * sfx;
      if ( ( sfx > (1.+DBL_EPSILON) * GEN.um )   /* avoid roundoff error with FP registers */
	   || (xfx < (1.-DBL_EPSILON) * GEN.vl) 
	   || (xfx > (1.+DBL_EPSILON) * GEN.vr) )
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"pdf(x) > hat(x)");

      /* evaluate squeeze */
      if ( (gen->variant & SROU_VARFLAG_SQUEEZE) &&
	   (X >= GEN.xl) && 
	   (X <= GEN.xr ) && 
	   (U < GEN.um) ) {
	xx = V / (GEN.um - U);
	if ( (xx >= GEN.xl) && (xx <= GEN.xr ) )
	  return x;
      }

      /* accept or reject */
      if (U*U <= PDF(x))
	return x;
    }
  }

} /* end of _unur_srou_sample_check() */

/*****************************************************************************/

void
_unur_srou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_SROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SROU_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of _unur_srou_free() */

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
_unur_srou_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_SROU_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SROU_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = srou (simple universal ratio-of-uniforms)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = _unur_srou_sample",gen->genid);
  if (par->variant & SROU_VARFLAG_VERIFY)
    fprintf(log,"_check");
  else if (par->variant & SROU_VARFLAG_MIRROR)
    fprintf(log,"_mirror");
  fprintf(log,"()\n%s:\n",gen->genid);

  if (par->set & SROU_SET_CDFMODE)
    fprintf(log,"%s: F(mode) = %g\n",gen->genid,PAR.Fmode);
  else
    fprintf(log,"%s: F(mode) unknown\n",gen->genid);

  if (gen->variant & SROU_VARFLAG_SQUEEZE)
    fprintf(log,"%s: use universal squeeze\n",gen->genid);
  else
    fprintf(log,"%s: no (universal) squeeze\n",gen->genid);

  if (gen->variant & SROU_VARFLAG_MIRROR)
    fprintf(log,"%s: use mirror principle\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: Rectangle:\n",gen->genid);
  fprintf(log,"%s:    left upper point  = (%g,%g)\n",gen->genid,GEN.vl,GEN.um);
  fprintf(log,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN.vr,GEN.um);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_srou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
