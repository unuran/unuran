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
 *   [1] Leydold J. (2000): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions, preprint.              *
 *                                                                           *
 *   [2] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  ... Beschreibung ...                                                     *
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

int unur_empk_set_kernel( UNUR_PAR *parameters, UNUR_DISTR *kernel);
/* 
   Set kernel distribution.
   Default is a Gaussian kernel.
*/

/*---------------------------------------------------------------------------*/

int unur_empk_set_kernelgen( UNUR_PAR *parameters, UNUR_GEN *kernelgen);
/* 
   Set generator for the kernel used the density estimation.
   Default is a Gaussian kernel.
*/

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
