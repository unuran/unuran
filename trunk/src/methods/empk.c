/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      empk.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generate from kernel estimation                              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given observed sample.                                               *
 *      Produce a value x consistent with this sample                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to sample                                                    *
 *      alpha factor                                                         *
 *                                                                           *
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
 *   [1] ....                                                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  ... Beschreibung ...                                                     *
 *                                                                           *
 TODO: (bitte auf englisch)

 alpha(K) = Var(K)^(-2/5){ \int K(t)^2 dt}^(1/5)

 bwidth = alpha * beta * (mixture of stdev and interquartilsrange of sample)
 * smoothing_factor

 dann ist der default 1, und man kann, wenn man zB annimmt, dass die
 Verteilung der Daten
 multimodal ist, oder aus anderen Gruenden smoothing factor auch kleiner
 als 1 waehlen.
 
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define EMPK_VARFLAG_VARCOR     0x001u   /* use variance correction          */
#define EMPK_VARFLAG_POSITIVE   0x002u   /* use variance correction          */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define EMPK_SET_KERN           0x010u    /* kernel density                  */
#define EMPK_SET_KERNGEN        0x020u    /* generator for kernel            */
#define EMPK_SET_KERNELVAR      0x001u    /* variance of kernel              */
#define EMPK_SET_ALPHA          0x002u    /* alpha factor                    */
#define EMPK_SET_BETA           0x004u    /* beta factor                     */
#define EMPK_SET_SMOOTHING      0x008u    /* smoothing factor                */

/*---------------------------------------------------------------------------*/

#define GENTYPE "EMPK"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empk_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empk_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_empk_start_kernel( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* start generator for kernel distribution                                   */
/*---------------------------------------------------------------------------*/

/* No reinit() call                                                          */
/*  static int _unur_empk_reinit( struct unur_gen *gen );                    */
/*---------------------------------------------------------------------------*/
/* Re-initialize (existing) generator.                                       */
/*---------------------------------------------------------------------------*/

static double _unur_empk_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_empk_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_empk_debug_init( struct unur_par *par, struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cemp      /* data for distribution object      */

#define PAR       par->data.empk        /* data for parameter object         */
#define GEN       gen->data.empk        /* data for generator object         */
#define DISTR     gen->distr.data.cemp  /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_empk_new( struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CEMP,NULL);

  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_EMPK_PAR);

  /* copy input */
  par->distr    = distr;         /* pointer to distribution object           */

  /* set default values */
  PAR.kernvar   = 1.;            /* variance of used kernel, only used if varcor == 1 */
  PAR.alpha     = 0.7763884;     /* alpha for Gaussian kernel, efficiency = 0.951 */
  PAR.beta      = 1.3637439;     /* optimal BETA if data follow normal distribution */
  PAR.smoothing = 1.;            /* determines how "smooth" the estimated density will be */

  /* kernel */
  PAR.kern      = NULL;          /* kernel distribution                      */
  PAR.kerngen   = NULL;          /* random variate generator for kernel      */

  par->method   = UNUR_METH_EMPK;     /* method and default variant          */
  par->variant  = ( EMPK_VARFLAG_VARCOR );       /* default variant          */

  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_empk_init;

  return par;

} /* end of unur_empk_new() */

/*****************************************************************************/

int 
unur_empk_set_kernel( struct unur_par *par, struct unur_distr *kernel)
     /*----------------------------------------------------------------------*/
     /* set kernel distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   kernel ... distribution object that hold the kernel                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check kernel distribution */
  if (kernel->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return 0; }
  COOKIE_CHECK(kernel,CK_DISTR_CONT,0);

  /* set kernel distribution */
  PAR.kern = kernel;

  /* changelog */
  par->set |= EMPK_SET_KERN;

  /* o.k. */
  return 1;

} /* end of unur_empk_set_kernel() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_kernelgen( struct unur_par *par, struct unur_gen *kernelgen)
     /*----------------------------------------------------------------------*/
     /* set generator for kernel distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   kernelgen ... generator object that hold the kernel                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check kernel distribution */
  if ( (kernelgen->method & UNUR_MASK_TYPE) != UNUR_METH_CONT ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return 0; }

  /* set kernel distribution */
  PAR.kerngen = kernelgen;

  /* changelog */
  par->set |= EMPK_SET_KERNGEN;

  /* o.k. */
  return 1;

} /* end of unur_empk_set_kernelgen() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_alpha( struct unur_par *par, double alpha )
     /*----------------------------------------------------------------------*/
     /* set alpha factor                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   alpha    ... parameter depending on kernel                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check new parameter for generator */
  if (alpha <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"alpha <= 0");
    return 0;
  }

  /* store date */
  PAR.alpha = alpha;

  /* changelog */
  par->set |= EMPK_SET_ALPHA;

  return 1;

} /* end of unur_empk_set_alpha() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_beta( struct unur_par *par, double beta )
     /*----------------------------------------------------------------------*/
     /* set beta factor                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   beta     ... parameter depending on sample                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check new parameter for generator */
  if (beta <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"beta <= 0");
    return 0;
  }

  /* store date */
  PAR.beta = beta;

  /* changelog */
  par->set |= EMPK_SET_BETA;

  return 1;

} /* end of unur_empk_set_beta() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_smoothing( struct unur_par *par, double smoothing )
     /*----------------------------------------------------------------------*/
     /* set smoothing factor                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   smoothing ... smoothing factor                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check new parameter for generator */
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return 0;
  }

  /* store date */
  PAR.smoothing = smoothing;

  /* changelog */
  par->set |= EMPK_SET_SMOOTHING;

  return 1;

} /* end of unur_empk_set_smoothing() */


/*---------------------------------------------------------------------------*/

int
unur_empk_set_varcor( struct unur_par *par, int varcor )
     /*----------------------------------------------------------------------*/
     /* turn variance correction on/off                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   varcor   ... 0 = no variance correction,  !0 = variance correction */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   variance correction is the default                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* we use a bit in variant */
  par->variant = (varcor) 
    ? (par->variant | EMPK_VARFLAG_VARCOR) 
    : (par->variant & (~EMPK_VARFLAG_VARCOR));

  /* o.k. */
  return 1;

} /* end of unur_empk_set_varcor() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_kernelvar( struct unur_par *par, double kernelvar )
     /*----------------------------------------------------------------------*/
     /* set variance of kernel                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   kernelvar... variance of kernel                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* check new parameter for generator */
  if (kernelvar <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"variance <= 0");
    return 0;
  }

  /* store date */
  PAR.kernvar = kernelvar;

  /* changelog */
  par->set |= EMPK_SET_KERNELVAR;

  return 1;

} /* end of unur_empk_set_kernelvar() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_positive( struct unur_par *par, int positive )
     /*----------------------------------------------------------------------*/
     /* turn mirroring (produces positive random numbers only) on/off        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   positive ... 0 = no mirroring,  !0 = mirroring                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no mirroring is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,EMPK );

  /* we use a bit in variant */
  par->variant = (positive) 
    ? (par->variant | EMPK_VARFLAG_POSITIVE) 
    : (par->variant & (~EMPK_VARFLAG_POSITIVE));

  /* o.k. */
  return 1;

} /* end of unur_empk_set_positive() */

/*****************************************************************************/

struct unur_gen *
_unur_empk_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_EMPK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_empk_create(par);
  if (!gen) { free(par); return NULL; }

  if (!_unur_empk_start_kernel( par, gen )) 
    { free(par); _unur_empk_free(gen); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
/*      if (gen->debug) _unur_ssr_debug_init(par,gen); */
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of _unur_empk_init() */

#if 0
  GEN.kernrvg=PAR.kernrvg;/*random variate generator function for the kernel*/

  qksort(GEN.observ,0,(GEN.n-1));   /*the observations must be ordered first*/

#define BETA 1.3637439 /*optimal BETA if data follow normal distribution*/
  GEN.bwidth = PAR.alfa*BETA*
   iqrtilestdev(GEN.observ,GEN.n,&(GEN.xbar),&(GEN.stdev))/exp(0.2*log(GEN.n));
#undef BETA
  GEN.sconst=1./sqrt(1. + PAR.kernvar*SQ(GEN.bwidth/GEN.stdev));
#endif

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_empk_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_EMPK_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_empk_sample;
  gen->destroy = _unur_empk_free;
  gen->reinit = _unur_reinit_error;

  /* copy observed data into generator object */
  GEN.n_observ = par->distr->data.cemp.n_sample;     /* sample size */
  GEN.observ = _unur_malloc( GEN.n_observ * sizeof(double) );
  memcpy( GEN.observ, par->distr->data.cemp.sample, GEN.n_observ * sizeof(double) );
  DISTR.sample = GEN.observ;  /* update pointer in local distribution object */

  /* copy some parameters into generator object */
  GEN.kerngen = PAR.kerngen;        /* generator for kernel distribution     */

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_empk_create() */

/*---------------------------------------------------------------------------*/

int 
_unur_empk_start_kernel( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* start generator for kernel                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* get generator for kernel */

  if (GEN.kerngen == NULL) {
    /* no generator provided by user */

    if (PAR.kern == NULL) {
      /* no kernel distribution given --> use Gaussian kernel (default) */
      PAR.kern = unur_distr_normal(NULL,0);
      if (!PAR.kern) return 0;
    }

    /* check kernel distribution */
    if (PAR.kern->type != UNUR_DISTR_CONT) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_INVALID,""); return 0; }
    COOKIE_CHECK(PAR.kern,CK_DISTR_CONT,0);

    /* get generator object */
    GEN.kerngen = unur_init( unur_cstd_new( PAR.kern ) );
    if (!GEN.kerngen) {
      _unur_error(gen->genid,UNUR_ERR_PAR_INVALID,"cannot init kernel generator");
      return 0;
    }
  }

  /* check generator for kernel distribution */
  if ( (GEN.kerngen->method & UNUR_MASK_TYPE) != UNUR_METH_CONT ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return 0;
  }

  /* clear up */
  if ( !(par->set & EMPK_SET_KERN) )
    /* no kernel distribution given by user.
       delete automatically created distribution object. */
    unur_distr_free( PAR.kern );

  /* o.k. */
  return 1;

} /* _unur_empk_start_kernel() */


/*****************************************************************************/

double
_unur_empk_sample( struct unur_gen *gen )
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
  return 0.;
} /* end of _unur_empk_sample() */

/*****************************************************************************/

void
_unur_empk_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_EMPK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_EMPK_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN.observ) free(GEN.observ);

  if ( !(gen->set & EMPK_SET_KERNGEN) )
    /* no generator for kernel distribution given by user.
       delete automatically created generation object. */
    unur_free( GEN.kerngen );

  _unur_free_genid(gen);
  free(gen);

} /* end of _unur_empk_free() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
