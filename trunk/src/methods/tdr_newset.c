/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_newset.c                                                 *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a T-concave distribution                                *
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
 *****************************************************************************/

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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL; }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of PDF"); return NULL; }

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
  PAR.max_ivs             = 100;    /* maximum number of intervals           */
  PAR.max_ratio           = 0.95;   /* bound for ratio  Atotal / Asqueeze    */
  PAR.bound_for_adding    = 0.5;    /* do not add a new construction point in an interval,
				       where ambigous region is too small, i.e. if 
				       area / ((A_hat - A_squeeze)/number of segments) < bound_for_adding */
 
  par->method   = UNUR_METH_TDR;                 /* method                   */
  par->variant  = ( TDR_VARFLAG_USECENTER |      /* default variant          */
		    TDR_VARFLAG_USEMODE   |
		    TDR_VARFLAG_PEDANTIC  |
                    TDR_VARIANT_GW );

  par->set      = 0u;               /* inidicate default parameters          */    
  par->urng     = unur_get_default_urng(); /* use default URNG               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

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

double
unur_tdr_get_sqhratio( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio ... on success                                               */
     /*   0     ... on error                                                 */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_gen_object( gen,TDR );

  return (GEN.Asqueeze / GEN.Atotal);

} /* end of unur_tdr_get_sqhratio() */

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
     /* set center (approximate mode) of PDF                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   center ... center of PDF                                           */
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
unur_tdr_set_variant_gw( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use original variant with squeezes as proposed by Gilks & Wild       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
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

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_GW;

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_variant_gw() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_variant_ps( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use squeezes proportional to the hat function                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
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

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_PS;

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_variant_ps() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_variant_ia( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use squeezes proportional to the hat function together with a        */
     /* composition method that required less uniform random numbers.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
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

  /* we use a bit in variant */
  par->variant = (par->variant & ~TDR_VARMASK_VARIANT) | TDR_VARIANT_IA;

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_variant_ia() */

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

/*---------------------------------------------------------------------------*/

int
unur_tdr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen,TDR );

  /* we use a bit in variant */
  gen->variant = (verify) ? (gen->variant | TDR_VARFLAG_VERIFY) : (gen->variant & (~TDR_VARFLAG_VERIFY));

  /* sampling routines */
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    SAMPLE = (verify) ? _unur_tdr_gw_sample_check : _unur_tdr_gw_sample;
    break;
  case TDR_VARIANT_PS:    /* proportional squeeze */
    SAMPLE = (verify) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
    break;
  case TDR_VARIANT_IA:    /* immediate acceptance */
    SAMPLE = (verify) ? _unur_tdr_ia_sample_check : _unur_tdr_ia_sample;
    break;
  default:
    _unur_warning(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

  /* o.k. */
  return 1;

} /* end of unur_tdr_chg_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tdr_set_pedantic( struct unur_par *par, int pedantic )
     /*----------------------------------------------------------------------*/
     /* turn pedantic mode on/off                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   pedantic ... 0 = no pedantic mode, !0 = use pedantic mode          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   pedantic is the default                                            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,TDR );

  /* we use a bit in variant */
  par->variant = (pedantic) ? (par->variant | TDR_VARFLAG_PEDANTIC) : (par->variant & (~TDR_VARFLAG_PEDANTIC));

  /* o.k. */
  return 1;

} /* end of unur_tdr_set_pedantic() */

/*---------------------------------------------------------------------------*/

int 
unur_tdr_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
     /*                                                                      */
     /* This call does not work for variant IA (immediate acceptance).       */
     /* In this case it switches to variant PS!!                             */
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
  CHECK_NULL(gen, 0);
  _unur_check_gen_object(gen, TDR);

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    DISTR.trunc[0] = DISTR.domain[0];
  }
  else
    DISTR.trunc[0] = left;

  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    DISTR.trunc[1] = DISTR.domain[1];
  }
  else
    DISTR.trunc[1] = right;

  /* set bounds of U -- in respect to given bounds */
  GEN.Umin = _unur_tdr_eval_cdfhat(gen,DISTR.trunc[0]);
  GEN.Umax = _unur_tdr_eval_cdfhat(gen,DISTR.trunc[1]);

  /* we have to disable adaptive rejection sampling */
  GEN.max_ivs = GEN.n_ivs;

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
/*    if (gen->debug & NINV_DEBUG_CHG)  */
/*      _unur_ninv_debug_chg_truncated( gen ); */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_tdr_chg_truncated() */

/*---------------------------------------------------------------------------*/

static double
_unur_tdr_eval_cdfhat( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of hat at x (i.e. \int_\infty^x hat(t) dt)              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point at which hat(x) has to be computed                   */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF of hat(x) or                                                   */
     /*   0. in case of error                                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double Aint;
  double cdf;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* the easy case */
  if (DISTR.trunc[0] <= -INFINITY) return 0.;
  if (DISTR.trunc[1] >=  INFINITY) return 1.;

  /* there are differencies between variant GW and variant PS */
  switch (gen->variant & TDR_VARMASK_VARIANT) {

  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */

    /* find interval (sequential search) */
    for (iv = GEN.iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
      /* iv->x is left construction point of interval */
      if (x < iv->next->x) break;
    }

    if (iv->next == NULL)
      /* right boundary of domain */
      return 1.;

    /* now iv->x < x <= iv->next->x */

    /* compute are below hat between construction point of tangent and x. */
    /* we have to cases on either side of the intersection point.         */

    if (x < iv->ip) {
      /* left h.s. of intersection point */
      Aint = _unur_tdr_interval_area( gen, iv, iv->dTfx, x);
      cdf = iv->Acum - iv->Ahat + Aint;
      if (cdf < 0.) cdf = 0.;
    }
    else {
      /* right h.s. of intersection point */
      Aint = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, x);
      cdf = iv->Acum - Aint;
    }

    /* normalize to one (and mind round-off errors) */
    cdf /= GEN.Atotal;
    return ((cdf > 1.) ? 1. : cdf);

    
  case TDR_VARIANT_PS:    /* proportional squeeze */
#if 0
    /* find interval (sequential search) */
    for (iv = GEN.iv; iv->next!=NULL; iv=iv->next) {
      COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 
      if (x > iv->ip) break;
    }
#endif
    return 1.;

  case TDR_VARIANT_IA:    /* immediate acceptance */
    /* this variant is not a pure rejection algorithm, but a
       composition method. Thus it does not make to much sense
       to compute the "CDF" of the hat.
    */
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0.;
  }

  return 1.;


} /* end of _unur_tdr_eval_cdfhat() */

/*****************************************************************************/

