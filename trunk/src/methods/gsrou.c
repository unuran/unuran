/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gsrou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generalize simple universal method (ratio-of-uniforms method)*
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and mode of a T_c-concave distribution                     *
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
 *   [1] Leydold J. (2002): Short universal generators via generalized       *
 *       ratio-of-uniforms method, Preprint.                                 *
 *                                                                           *
 *   [2] Wakefield J. C., Gelfand A. E., and Smith A. F. M. (1991):          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method, Statist. Comput. 1(2), pp. 129--133.                        *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define GSROU_VARFLAG_VERIFY  0x002u   /* run verify mode                    */
#define GSROU_VARFLAG_MIRROR  0x008u   /* use mirror principle               */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define GSROU_DEBUG_REINIT    0x00000010u  /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define GSROU_SET_CDFMODE     0x001u   /* CDF at mode is known               */
#define GSROU_SET_R           0x002u   /* parameter r for power transform.   */

/*---------------------------------------------------------------------------*/

#define GENTYPE "GSROU"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gsrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_gsrou_envelope( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute parameters for universal bounding envelope.                       */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gsrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_gsrou_sample_mirror( struct unur_gen *gen );
static double _unur_gsrou_sample_check( struct unur_gen *gen );
static double _unur_gsrou_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_gsrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_gsrou_debug_init( struct unur_gen *gen, int is_reinit );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.gsrou       /* data for parameter object         */
#define GEN       gen->data.gsrou       /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),&(gen->distr))   /* call to PDF         */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_gsrou_new( struct unur_distr *distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (!unur_distr_cont_upd_mode(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return NULL; 
    }
  }

  if (!(distr->set & UNUR_DISTR_SET_PDFAREA))
    if (!unur_distr_cont_upd_pdfarea(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF");
      return NULL; 
    }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_GSROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR.r         = 2.;                 /* parameter for power transformation  */
  PAR.Fmode     = -1.;                /* CDF at mode (unknown yet   )        */

  par->method   = UNUR_METH_GSROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_gsrou_init;

  return par;

} /* end of unur_gsrou_new() */

/*****************************************************************************/

int
unur_gsrou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* set parameter r for power transformation                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   r    ... parameter r                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,GSROU );

  /* check new parameter for generator */
  if (r < 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r < 1");
    return 0;
  }
    
  /* store date */
  PAR.r = r;

  /* changelog */
  par->set |= GSROU_SET_R;

  return 1;

} /* end of unur_gsrou_set_r() */

/*---------------------------------------------------------------------------*/

int 
unur_gsrou_set_cdfatmode( struct unur_par *par, double Fmode )
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
  _unur_check_par_object( par,GSROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return 0;
  }

  /* store date */
  PAR.Fmode = Fmode;

  /* changelog */
  par->set |= GSROU_SET_CDFMODE;

  return 1;

} /* end of unur_gsrou_set_cdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par,GSROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | GSROU_VARFLAG_VERIFY) : (par->variant & (~GSROU_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_gsrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen,GSROU );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= GSROU_VARFLAG_VERIFY;
    SAMPLE = _unur_gsrou_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~GSROU_VARFLAG_VERIFY;
    SAMPLE = (gen->variant & GSROU_VARFLAG_MIRROR) ? _unur_gsrou_sample_mirror : _unur_gsrou_sample;
  }

  /* o.k. */
  return 1;

} /* end of unur_gsrou_chg_verify() */

/*---------------------------------------------------------------------------*/

int 
unur_gsrou_set_usemirror( struct unur_par *par, int usemirror )
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
  _unur_check_par_object( par,GSROU );

  /* we use a bit in variant */
  par->variant = (usemirror) ? (par->variant | GSROU_VARFLAG_MIRROR) : (par->variant & (~GSROU_VARFLAG_MIRROR));

  /* o.k. */
  return 1;

} /* end of unur_gsrou_set_usemirror() */

/*****************************************************************************/

int
unur_gsrou_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  _unur_check_gen_object( gen,GSROU );
  
  /* set new parameters in distribution object */
  return unur_distr_cont_set_pdfparams(&(gen->distr),params,n_params);

} /* end of unur_gsrou_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_chg_mode( struct unur_gen *gen, double mode )
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
  _unur_check_gen_object( gen,GSROU );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_gsrou_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_upd_mode( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,GSROU );

  return unur_distr_cont_upd_mode( &(gen->distr) );
} /* end of unur_gsrou_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
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
  _unur_check_gen_object( gen,GSROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return 0;
  }
  
  /* copy parameters */
  GEN.Fmode = Fmode;

  /* changelog */
  gen->set |= GSROU_SET_CDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_gsrou_chg_cdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_gsrou_chg_domain( struct unur_gen *gen, double left, double right )
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
  _unur_check_gen_object( gen,GSROU );

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
  
} /* end of unur_gsrou_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* change area below PDF of distribution                                */
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
  _unur_check_gen_object( gen,GSROU );
  
  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"area <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_gsrou_chg_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_upd_pdfarea( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute area below PDF of distribution                             */
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
  _unur_check_gen_object( gen,GSROU );

  return unur_distr_cont_upd_pdfarea( &(gen->distr) );
} /* end of unur_gsrou_upd_pdfarea() */

/*****************************************************************************/

struct unur_gen *
_unur_gsrou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_GSROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_GSROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_gsrou_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute parameters for universal bounding envelope */
  if (! _unur_gsrou_envelope( gen ) ) {
    free(par); _unur_gsrou_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_gsrou_debug_init(gen, FALSE);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_gsrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_gsrou_envelope( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute parameters for universal bounding envelope.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{ 
  double fm;             /* PDF at mode               */
  double vm;             /* maximal width of envelope */
  double pr;             /* p^r                       */

  double p;              /* short cuts                */
  double r = GEN.r;

  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_GSROU_GEN, 0 );

  /* compute PDF at mode */
  fm = PDF(DISTR.mode);
  /* fm must be positive */
  if (fm <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
    return 0;
  }
  if (_unur_FP_is_infinity(fm)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return 0;
  }
  GEN.um = pow(fm,1./(r+1.));    /* height of envelope */

  /* maximal width of envelope */
  vm = DISTR.area / (GEN.r*GEN.um);

  if (gen->set & GSROU_SET_CDFMODE) {
    /* cdf at mode known */
    GEN.vl = -GEN.Fmode * vm;
    GEN.vr = vm + GEN.vl;
/*      GEN.xl = GEN.vl/GEN.um; */
/*      GEN.xr = GEN.vr/GEN.um; */
  }
  else {
    /* cdf at mode unknown */
    GEN.vl = -vm;
    GEN.vr = vm;
/*      GEN.xl = GEN.vl/GEN.um; */
/*      GEN.xr = GEN.vr/GEN.um; */
  }

  /* construction point for bounding curve */
  GEN.p = p = 1. - 2.187/pow(r + 5 - 1.28/r, 0.9460 );
  pr = pow(p,r);

  /* parameters for bounding envelope */
  GEN.b = (1. - r * pr/p + (r-1.)*pr) / ((pr-1.)*(pr-1));
  GEN.a = -(p-1.)/(pr-1.) - p * GEN.b;
  GEN.log_ab = log(GEN.a/(GEN.a+GEN.b));   /** TODO **/

  /* o.k. */
  return 1;

} /* end of _unur_gsrou_envelope() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_gsrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_GSROU_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_GSROU_GEN);

  /* copy distribution object into generator object */
  _unur_distr_cont_copy( &(gen->distr), par->distr );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  if (par->variant & GSROU_VARFLAG_VERIFY)
    SAMPLE = _unur_gsrou_sample_check;
  else
    SAMPLE = (par->variant & GSROU_VARFLAG_MIRROR) ? _unur_gsrou_sample_mirror : _unur_gsrou_sample;

  gen->destroy = _unur_gsrou_free;

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
  GEN.r     = PAR.r;                /* parameter for power transformation    */
  GEN.Fmode = PAR.Fmode;            /* CDF at mode                           */

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

} /* end of _unur_gsrou_create() */

/*---------------------------------------------------------------------------*/

int
unur_gsrou_reinit( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,GSROU );

  /* compute parameters for universal bounding envelope */
  result = _unur_gsrou_envelope( gen );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & GSROU_DEBUG_REINIT)
    if (gen->debug) _unur_gsrou_debug_init(gen,TRUE);
#endif

  return result;
} /* end of unur_gsrou_reinit() */

/*****************************************************************************/

double
_unur_gsrou_sample( struct unur_gen *gen )
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
  double U,Ur,V,W,X,Z;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_GSROU_GEN,0.);

  while (1) {
    W = GEN.log_ab *_unur_call_urng(gen->urng);
    Z = GEN.vl + _unur_call_urng(gen->urng) * (GEN.vr - GEN.vl);
    U = (exp(-W)-1.) * GEN.a/GEN.b;
    V = -Z/(GEN.a + GEN.b*U);
    U *= GEN.um;
    Ur = pow(U,GEN.r);
    X = V/Ur + DISTR.mode;
    /* accept or reject */
    if (Ur*U <= PDF(X))
      return X;
  }

} /* end of _unur_gsrou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_gsrou_sample_mirror( struct unur_gen *gen )
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
/*    double U,V,X,x,fx,fnx,uu; */

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_GSROU_GEN,0.);

/*    while (1) { */
    /* generate point uniformly on rectangle */
/*      while ( (U = _unur_call_urng(gen->urng)) == 0.); */
/*      U *= GEN.um * SQRT2; */
/*      V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN.vr; */
    /* vr = vm when the CDF at the mode is not known */

    /* ratio */
/*      X = V/U; */

    /* evaluate PDF */
/*      x = X + DISTR.mode; */
/*      fx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x); */
/*      uu = U * U; */

    /* accept or reject */
/*      if (uu <= fx) */
/*        return x; */

    /* try mirrored PDF */
/*      x = -X + DISTR.mode; */
/*      fnx  = (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) ? 0. : PDF(x); */
/*      if (uu <= fx + fnx) */
/*        return x; */
/*    } */

    return 0.;
} /* end of _unur_gsrou_sample_mirror() */

/*---------------------------------------------------------------------------*/

double
_unur_gsrou_sample_check( struct unur_gen *gen )
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
/*    double U,uu,V,X,x,nx,fx,sfx,fnx,xfx,xfnx,xx; */

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_GSROU_GEN,0.);

/*    if (gen->variant & SROU_VARFLAG_MIRROR) { */
    /* use mirror principle */

/*      while (1) { */
      /* generate point uniformly on rectangle */
/*        while ( (U = _unur_call_urng(gen->urng)) == 0.); */
/*        U *= GEN.um * SQRT2; */
/*        V = 2. * (_unur_call_urng(gen->urng) - 0.5) * GEN.vr; */
      /* vr = vm when the CDF at the mode is not known */	

      /* ratio */
/*        X = V/U; */

      /* x values */
/*        x = X + DISTR.mode; */
/*        nx = -X + DISTR.mode; */

      /* evaluate PDF */
/*        fx  = (x  < DISTR.BD_LEFT || x  > DISTR.BD_RIGHT) ? 0. : PDF(x); */
/*        fnx = (nx < DISTR.BD_LEFT || nx > DISTR.BD_RIGHT) ? 0. : PDF(nx); */
/*        uu = U * U; */

      /* check hat */
/*        xfx  = (x  - DISTR.mode) * sqrt(fx); */
/*        xfnx = (nx - DISTR.mode) * sqrt(fnx); */

/*        if ( ((2.+4.*DBL_EPSILON) * GEN.um*GEN.um < fx + fnx)    * avoid roundoff error with FP registers */
/*  	   || (xfx < (1.+UNUR_EPSILON) * GEN.vl)  */
/*  	   || (xfx > (1.+UNUR_EPSILON) * GEN.vr) */
/*  	   || (xfnx < (1.+UNUR_EPSILON) * GEN.vl)  */
/*  	   || (xfnx > (1.+UNUR_EPSILON) * GEN.vr) ) */
/*  	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)"); */

      /* accept or reject */
/*        if (uu <= fx) */
/*  	return x; */
      
      /* try mirrored PDF */
/*        if (uu <= fx + fnx) */
/*  	return nx; */
/*      } */
/*    } */

/*  else { * do not use mirror principle */

/*      while (1) { */
      /* generate point uniformly on rectangle */
/*        while ( (U = _unur_call_urng(gen->urng)) == 0.); */
/*        U *= GEN.um; */
/*        V = GEN.vl + _unur_call_urng(gen->urng) * (GEN.vr - GEN.vl); */

      /* ratio */
/*        X = V/U; */

      /* compute x */
/*        x = X + DISTR.mode; */

      /* inside domain ? */
/*        if ( (x < DISTR.BD_LEFT) || (x > DISTR.BD_RIGHT) ) */
/*  	continue; */

      /* evaluate PDF */
/*        fx = PDF(x); */

      /* the point on the boundary of the region of acceptance
	 in direction X = V/U has the coordinates
	 ( X * sqrt(fx), sqrt(fx) ). */
/*        sfx = sqrt(fx); */
/*        xfx = X * sfx; */

      /* check hat */
/*      if ( ( sfx > (1.+DBL_EPSILON) * GEN.um )   * avoid roundoff error with FP registers */
/*  	   || (xfx < (1.+UNUR_EPSILON) * GEN.vl)  */
/*  	   || (xfx > (1.+UNUR_EPSILON) * GEN.vr) ) */
/*  	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)"); */

      /* evaluate and check squeeze */
/*        if ( (gen->variant & SROU_VARFLAG_SQUEEZE) && */
/*  	   (X >= GEN.xl) &&  */
/*  	   (X <= GEN.xr ) &&  */
/*  	   (U < GEN.um) ) { */

	/* check squeeze */
/*  	xx = xfx / (GEN.um - sfx); */
/*  	if ( (xx > (1.-UNUR_EPSILON) * GEN.xl) && */
/*  	     (xx < (1.-UNUR_EPSILON) * GEN.xr) ) */
/*  	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)"); */

	/* squeeze acceptance */
/*  	xx = V / (GEN.um - U); */
/*  	if ( (xx >= GEN.xl) && (xx <= GEN.xr ) ) */
/*  	  return x; */
/*        } */
      
      /* accept or reject */
/*        if (U*U <= PDF(x)) */
/*  	return x; */
/*      } */
/*    } */

    return 0.;

} /* end of _unur_gsrou_sample_check() */

/*****************************************************************************/

void
_unur_gsrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_GSROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_GSROU_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_distr_cont_clear(gen);
  _unur_free_genid(gen);

  free(gen);

} /* end of _unur_gsrou_free() */

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
_unur_gsrou_debug_init( struct unur_gen *gen, int is_reinit )
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_GSROU_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
    fprintf(log,"%s: method  = gsrou (generalized simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(log,"%s: reinit!\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = _unur_gsrou_sample",gen->genid);
  if (gen->variant & GSROU_VARFLAG_VERIFY)
    fprintf(log,"_check");
  else if (gen->variant & GSROU_VARFLAG_MIRROR)
    fprintf(log,"_mirror");
  fprintf(log,"()\n%s:\n",gen->genid);

  if (gen->set & GSROU_SET_CDFMODE)
    fprintf(log,"%s: F(mode) = %g\n",gen->genid,GEN.Fmode);
  else
    fprintf(log,"%s: F(mode) unknown\n",gen->genid);

  fprintf(log,"%s: no (universal) squeeze\n",gen->genid);

  if (gen->variant & GSROU_VARFLAG_MIRROR)
    fprintf(log,"%s: use mirror principle\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);

/*    fprintf(log,"%s: Rectangle:\n",gen->genid); */
/*    fprintf(log,"%s:    left upper point  = (%g,%g)\n",gen->genid,GEN.vl,GEN.um); */
/*    fprintf(log,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN.vr,GEN.um); */

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_gsrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
