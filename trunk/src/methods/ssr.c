/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ssr.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection with universal bounds                              *
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
 *   [1] Leydold J. (2001): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions,                        *
 *       ACM Trans. Math. Software 27(1), pp. 66--82.                        *
 *                                                                           *
 *   [2] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * This algorithm is based on transformed density rejection (see [2]), were  *
 * universal upper and lower bounds are used. These are derived via the      *
 * ratio-of-uniforms method. See [1] for details and a description of the    *
 * algorithm. It works for any distribution, where -1/sqrt(PDF(x)) is        *
 * concave. This includes all log-concave distributions.                     *
 *                                                                           *
 * (The mirror principle has not been implemented.)                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define SSR_VARFLAG_VERIFY   0x002u    /* run verify mode                    */
#define SSR_VARFLAG_SQUEEZE  0x004u    /* use universal squeeze if possible  */

#define SSR_VARFLAG_MIRROR   0x008u    /* use mirror principle               */
/* not implemented yet! */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define SSR_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define SSR_SET_CDFMODE      0x001u    /* CDF at mode is known               */
#define SSR_SET_PDFMODE      0x002u    /* PDF at mode is set                 */

/*---------------------------------------------------------------------------*/

#define GENTYPE "SSR"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ssr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_ssr_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute universal hat.                                                    */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ssr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_ssr_sample( struct unur_gen *gen );
static double _unur_ssr_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_ssr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_ssr_debug_init( const struct unur_gen *gen, int is_reinit );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.ssr         /* data for parameter object         */
#define GEN       gen->data.ssr         /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2    1.4142135623731

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_ssr_new( const struct unur_distr *distr )
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

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_SSR_PAR);

  /* copy input */
  par->distr    = distr;              /* pointer to distribution object      */

  /* set default values */
  PAR.Fmode     = -1.;                /* CDF at mode (unknown yet)           */
  PAR.fm        = -1.;                /* PDF at mode (unknown)               */
  PAR.um        = -1.;                /* square of PDF at mode (unknown)     */

  par->method   = UNUR_METH_SSR;      /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_ssr_init;

  return par;

} /* end of unur_ssr_new() */

/*****************************************************************************/

int 
unur_ssr_set_cdfatmode( struct unur_par *par, double Fmode )
     /*----------------------------------------------------------------------*/
     /* set cdf at mode                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   Fmode ... CDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SSR );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return 0;
  }

  /* store date */
  PAR.Fmode = Fmode;

  /* changelog */
  par->set |= SSR_SET_CDFMODE;

  return 1;

} /* end of unur_ssr_set_cdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_ssr_set_pdfatmode( UNUR_PAR *par, double fmode )
     /*----------------------------------------------------------------------*/
     /* Set PDF at mode. if set the PDF at the mode is never changed.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SSR );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode)");
    return 0;
  }
  if (_unur_FP_is_infinity(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return 0;
  }

  /* store date */
  PAR.fm = fmode;
  PAR.um = sqrt(fmode);

  /* changelog */
  par->set |= SSR_SET_PDFMODE;

  return 1;

} /* end of unur_ssr_set_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par,SSR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | SSR_VARFLAG_VERIFY) : (par->variant & (~SSR_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_ssr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen,SSR );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= SSR_VARFLAG_VERIFY;
    SAMPLE = _unur_ssr_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~SSR_VARFLAG_VERIFY;
    SAMPLE = _unur_ssr_sample;
  }

  /* o.k. */
  return 1;

} /* end of unur_ssr_chg_verify() */

/*---------------------------------------------------------------------------*/

int 
unur_ssr_set_usesqueeze( struct unur_par *par, int usesqueeze )
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
  _unur_check_par_object( par,SSR );

  /* we use a bit in variant */
  par->variant = (usesqueeze) ? (par->variant | SSR_VARFLAG_SQUEEZE) : (par->variant & (~SSR_VARFLAG_SQUEEZE));

  /* o.k. */
  return 1;

} /* end of unur_ssr_set_usesqueeze() */

/*****************************************************************************/

int
unur_ssr_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  _unur_check_gen_object( gen,SSR );
  
  /* set new parameters in distribution object */
  return unur_distr_cont_set_pdfparams( gen->distr, params,n_params );

} /* end of unur_ssr_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_mode( struct unur_gen *gen, double mode )
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
  _unur_check_gen_object( gen,SSR );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_ssr_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_upd_mode( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,SSR );

  return unur_distr_cont_upd_mode( gen->distr );
} /* end of unur_ssr_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_cdfatmode( struct unur_gen *gen, double Fmode )
     /*----------------------------------------------------------------------*/
     /* change value of CDF at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   Fmode ... CDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SSR );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return 0;
  }
  
  /* copy parameters */
  GEN.Fmode = Fmode;

  /* changelog */
  gen->set |= SSR_SET_CDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_ssr_chg_cdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_pdfatmode( struct unur_gen *gen, double fmode )
     /*----------------------------------------------------------------------*/
     /* change value of PDF at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   fmode ... PDF at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,SSR );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"PDF(mode)");
    return 0;
  }
  if (_unur_FP_is_infinity(fmode)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
    return 0;
  }

  /* store date */
  GEN.fm = fmode;
  GEN.um = sqrt(fmode);

  /* changelog */
  gen->set |= SSR_SET_PDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_ssr_chg_pdfatmode() */

/*---------------------------------------------------------------------------*/

int 
unur_ssr_chg_domain( struct unur_gen *gen, double left, double right )
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
  _unur_check_gen_object( gen,SSR );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN | UNUR_DISTR_SET_MASK_DERIVED );
  gen->distr->set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_ssr_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_chg_pdfarea( struct unur_gen *gen, double area )
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
  _unur_check_gen_object( gen,SSR );
  
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
} /* end of unur_ssr_chg_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_upd_pdfarea( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,SSR );

  return unur_distr_cont_upd_pdfarea( gen->distr );
} /* end of unur_ssr_upd_pdfarea() */

/*****************************************************************************/

struct unur_gen *
_unur_ssr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_SSR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SSR_PAR,NULL);

  if (! (par->set & SSR_SET_CDFMODE))
    /* cdf at mode unknown -->
       thus we cannot use universal squeeze */
    par->variant &= ~SSR_VARFLAG_SQUEEZE;

  /* create a new empty generator object */
  gen = _unur_ssr_create(par);
  if (!gen) { free(par); return NULL; }

  /* compute universal bounding rectangle */
  if (! _unur_ssr_hat( gen ) ) {
    free(par); _unur_ssr_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_ssr_debug_init(gen, FALSE);
#endif

  /* free parameters */
  free(par);

  return gen;
} /* end of _unur_ssr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_ssr_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal hat                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{ 
  double vm, fm;             /* width of rectangle, PDF at mode              */
  double left,right;

  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_SSR_GEN, 0 );

  /* compute PDF at mode (if not given by user) */
  if (!(gen->set & SSR_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"PDF(mode) <= 0.");
      return 0;
    }
    if (_unur_FP_is_infinity(fm)) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF(mode) overflow");
      return 0;
    }
    GEN.fm = fm;        /* PDF at mode */
    GEN.um = sqrt(fm);  /* square root of PDF at mode */
  }

  /* compute parameters */
  vm = DISTR.area / GEN.um;

  if (gen->set & SSR_SET_CDFMODE) {
    /* cdf at mode known */
    GEN.vl = -GEN.Fmode * vm;
    GEN.vr = vm + GEN.vl;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
    GEN.A  = 2 * DISTR.area;
    GEN.al = (DISTR.BD_LEFT  < DISTR.mode) ? (GEN.Fmode * DISTR.area) : 0.;
    GEN.ar = (DISTR.BD_RIGHT > DISTR.mode) ? (GEN.al + DISTR.area) : GEN.A;
    /* Compute areas below hat in left tails and inside domain of PDF */
    if ( (DISTR.BD_LEFT > -INFINITY) &&
	 (DISTR.BD_LEFT < DISTR.mode) )
      GEN.Aleft = GEN.vl * GEN.vl / (DISTR.mode - DISTR.BD_LEFT);
    else
      GEN.Aleft = 0.;
    
    if ( (DISTR.BD_RIGHT < INFINITY) &&
	 (DISTR.BD_RIGHT > DISTR.mode) )
      GEN.Ain = GEN.A - GEN.vr * GEN.vr / (DISTR.BD_RIGHT - DISTR.mode);
    else
      GEN.Ain = GEN.A;
    GEN.Ain -= GEN.Aleft;
  }

  else {
    /* cdf at mode unknown */
    GEN.vl = -vm;
    GEN.vr = vm;
    GEN.xl = GEN.vl/GEN.um;
    GEN.xr = GEN.vr/GEN.um;
    GEN.A  = 4 * DISTR.area;
    GEN.al = DISTR.area;
    GEN.ar = 3 * DISTR.area;
    /* Compute areas below hat in left tails and inside domain of PDF */
    if (DISTR.BD_LEFT > -INFINITY) {
      left = DISTR.BD_LEFT - DISTR.mode;
      GEN.Aleft = (GEN.xl > left) 
	? (GEN.vl * GEN.vl / (-left)) 
	: (GEN.al + GEN.fm * (left - GEN.xl));
    }
    else 
      GEN.Aleft = 0.;
    
    if (DISTR.BD_RIGHT < INFINITY) {
      right = DISTR.BD_RIGHT - DISTR.mode;
      GEN.Ain = (GEN.xr < right) 
	? (GEN.A - GEN.vr * GEN.vr / right)
	: (GEN.ar - GEN.fm * (GEN.xr - right));
    }
    else 
      GEN.Ain = GEN.A;
    GEN.Ain -= GEN.Aleft;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif

  return 1;
} /* end of _unur_ssr_hat() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_ssr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SSR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_SSR_GEN);

  /* copy distribution object into generator object */
  gen->distr = _unur_distr_cont_clone( par->distr );

  /* check for required data: mode */
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (!unur_distr_cont_upd_mode( gen->distr )) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      _unur_distr_free(gen->distr); free(gen);
      return NULL; 
    }
  }

  /* check for required data: area */
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA))
    if (!unur_distr_cont_upd_pdfarea( gen->distr )) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF");
      _unur_distr_free(gen->distr); free(gen);
      return NULL; 
    }

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & SSR_VARFLAG_VERIFY) ? _unur_ssr_sample_check : _unur_ssr_sample;
  gen->destroy = _unur_ssr_free;
  gen->clone = _unur_ssr_clone;

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
  GEN.Fmode = PAR.Fmode;            /* CDF at mode                           */
  GEN.fm    = PAR.fm;               /* PDF at mode                           */
  GEN.um    = PAR.um;               /* square root of PDF at mode            */

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

} /* end of _unur_ssr_create() */

/*---------------------------------------------------------------------------*/

int
unur_ssr_reinit( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,SSR );

  /* compute universal bounding rectangle */
  result = _unur_ssr_hat( gen );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & SSR_DEBUG_REINIT)
    if (gen->debug) _unur_ssr_debug_init(gen,TRUE);
#endif

  return result;
} /* end of unur_ssr_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ssr_clone( const struct unur_gen *gen )
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
#define CLONE clone->data.ssr

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SSR_GEN,NULL);

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  clone->distr = _unur_distr_cont_clone( gen->distr );

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = unur_gen_clone( gen->gen_aux );

  return clone;

#undef CLONE
} /* end of _unur_ssr_clone() */

/*****************************************************************************/

double
_unur_ssr_sample( struct unur_gen *gen )
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
  double U,V,X,xx,y;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SSR_GEN,0.);

  while (1) {
    /* uniform ~U(0,1) */
    while ( (U = GEN.Aleft + _unur_call_urng(gen->urng) * GEN.Ain) == 0. );

    if (U < GEN.al) {        /* first part */
      X = - GEN.vl * GEN.vl / U;
      y = (U / GEN.vl);
      y = y*y;
    }
    else if (U <= GEN.ar) {  /* second part */
      X = GEN.xl + (U-GEN.al)/GEN.fm;
      y = GEN.fm;
    }
    else {                   /* third part */
      X = GEN.vr * GEN.vr / (GEN.um * GEN.vr - (U-GEN.ar));
      y = (GEN.A - U) / GEN.vr;
      y = y*y;
    }

    /* accept or reject */
    V = _unur_call_urng(gen->urng);
    y *= V;

    /* evaluate squeeze */
    if (gen->variant & SSR_VARFLAG_SQUEEZE) {
      xx = 2 * X;
      if ( xx >= GEN.xl && xx <= GEN.xr && y <= GEN.fm/4. )
	return (X + DISTR.mode);
    }

    /* Compute X */
    X += DISTR.mode;

    /* evaluate PDF */
    if (y <= PDF(X))
      return X;
  }
    
} /* end of _unur_ssr_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_ssr_sample_check( struct unur_gen *gen )
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
  double U,V,X,xx,fx,y;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SSR_GEN,0.);

  while (1) {
    /* uniform ~U(0,1) */
    while ( (U = GEN.Aleft + _unur_call_urng(gen->urng) * GEN.Ain) == 0. );

    if (U < GEN.al) {        /* first part */
      X = - GEN.vl * GEN.vl / U;
      y = (U / GEN.vl);
      y = y*y;
    }
    else if (U <= GEN.ar) {  /* second part */
      X = GEN.xl + (U-GEN.al)/GEN.fm;
      y = GEN.fm;
    }
    else {                   /* third part */
      X = GEN.vr * GEN.vr / (GEN.um * GEN.vr - (U-GEN.ar));
      y = (GEN.A - U) / GEN.vr;
      y = y*y;
    }

    /* compute PDF at x */
    fx = PDF(X + DISTR.mode);

    /* verify hat function */
    if ( (1.+UNUR_EPSILON) * y < fx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");

    /* accept or reject */
    V = _unur_call_urng(gen->urng);
    y *= V;

    /* evaluate and check squeeze */
    if (gen->variant & SSR_VARFLAG_SQUEEZE) {
      xx = 2 * X;
      if ( xx >= GEN.xl && xx <= GEN.xr ) {
	/* check squeeze */
	if ( fx < (1.-UNUR_EPSILON) * GEN.fm/4. )
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");
	/* evaluate squeeze */
	if ( y <= GEN.fm/4. )
	  return (X + DISTR.mode);
      }
    }

    /* Compute X */
    X += DISTR.mode;

    /* evaluate PDF */
    if (y <= fx)
      return X;
  }
    
} /* end of _unur_ssr_sample_check() */

/*****************************************************************************/

void
_unur_ssr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_SSR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SSR_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_distr_free(gen->distr);
  _unur_free_genid(gen);

  COOKIE_CLEAR(gen);
  free(gen);

} /* end of _unur_ssr_free() */

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
_unur_ssr_debug_init( const struct unur_gen *gen, int is_reinit )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   is_reinit ... if TRUE the generator has been reinitialized         */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SSR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
    fprintf(log,"%s: method  = ssr (simple universal transformed density rection)\n",gen->genid);
  }
  else
    fprintf(log,"%s: reinit!\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_ssr_sample",gen->genid);
  if (gen->variant & SSR_VARFLAG_VERIFY)
    fprintf(log,"_check");
  /* else if (gen->variant & SSR_VARFLAG_MIRROR)     not implemented */
  /*   fprintf(log,"_mirror"); */
  fprintf(log,"()\n%s:\n",gen->genid);

  if (gen->set & SSR_SET_CDFMODE)
    fprintf(log,"%s: CDF at mode = %g\n",gen->genid,GEN.Fmode);
  else
    fprintf(log,"%s: CDF at mode unknown\n",gen->genid);

  if (gen->variant & SSR_VARFLAG_SQUEEZE)
    fprintf(log,"%s: use universal squeeze\n",gen->genid);
  else
    fprintf(log,"%s: no (universal) squeeze\n",gen->genid);

  if (gen->variant & SSR_VARFLAG_MIRROR)
    fprintf(log,"%s: use mirror principle\n",gen->genid);

  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: parts:\n",gen->genid);
  fprintf(log,"%s:\txl = %g\n",gen->genid,GEN.xl);
  fprintf(log,"%s:\txr = %g\n",gen->genid,GEN.xr);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: PDF at mode:\n",gen->genid);
  fprintf(log,"%s:\tfm = %g\n",gen->genid,GEN.fm);
  fprintf(log,"%s:\tum = %g\n",gen->genid,GEN.um);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: areas:\n",gen->genid);
  fprintf(log,"%s:\t    al = %g\n",gen->genid,GEN.al);
  fprintf(log,"%s:\t    ar = %g\n",gen->genid,GEN.ar);
  fprintf(log,"%s:\t Aleft = %g\n",gen->genid,GEN.Aleft);
  fprintf(log,"%s:\t   Ain = %g\n",gen->genid,GEN.Ain);
  fprintf(log,"%s:\tAtotal = %g\n",gen->genid,GEN.A);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_ssr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
