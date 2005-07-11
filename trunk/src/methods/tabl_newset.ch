/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_newset.c                                                *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_tabl_new( const struct unur_distr *distr )
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

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_tabl_par) );
  COOKIE_SET(par,CK_TABL_PAR);

  /* copy input */
  par->distr        = distr;     /* pointer to distribution object           */

  /* set default values */
  PAR->slopes        = NULL;      /* pointer to slopes of PDF                 */
  PAR->n_slopes      = 0;         /* number of slopes                         */

  PAR->starting_cpoints    = NULL;   /* pointer to array of starting points   */
  PAR->n_starting_cpoints  = 30;     /* number of starting points             */
  PAR->area_fract    = 0.1;       /* parameter for equal area rule (default from [1] ) */

  PAR->max_ivs       = 1000;      /* maximum number of intervals              */
  PAR->max_ratio     = 0.90;      /* bound for ratio  Atotal / Asqueeze       */

  PAR->guide_factor  = 1.; /* guide table has same size as array of intervals */

  PAR->darsfactor    = 0.99;   /* factor for (derandomized) ARS.
				 do not add a new construction point in a interval
				 where abiguous region is too small          */

  /* default boundary of compution area */
  PAR->bleft     = -TABL_DEFAULT_COMPUTATION_LIMIT;
  PAR->bright    = TABL_DEFAULT_COMPUTATION_LIMIT;

  par->method   = UNUR_METH_TABL;              /* indicate method            */
  par->variant  = (TABL_VARFLAG_SPLIT_MEAN |   /* variant: split at arc_mean */
		   TABL_VARFLAG_STP_A      |   /* run SPLIT A on slopes      */
		   TABL_VARFLAG_USEDARS    );  /* run DARS (SPLIT B) on slopes */


  par->set      = 0u;                      /* inidicate default parameters   */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_tabl_init;

  return par;

} /* end of unur_tabl_new() */

/*****************************************************************************/

int
unur_tabl_set_usedars( struct unur_par *par, int usedars )
     /*----------------------------------------------------------------------*/
     /* set flag for using DARS (derandomized adaptive rejection sampling).  */
     /* additionally the rule for splitting intervals can be set.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usedars   ... 0 = do not use,  1 = use DARS                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using using DARS is the default                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* we use a bit in variant */
  if (usedars)
    par->variant |= TABL_VARFLAG_USEDARS;
  else
    par->variant &= ~TABL_VARFLAG_USEDARS;

  /* changelog */
  par->set |= TABL_SET_USE_DARS;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_set_usedars() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_darsfactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for derandomized adaptive rejection sampling              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... parameter for DARS                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
    
  /* store date */
  PAR->darsfactor = factor;

  /* changelog */
  par->set |= TABL_SET_DARS_FACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_darsfactor() */

/*---------------------------------------------------------------------------*/

int 
unur_tabl_set_variant_splitmode( struct unur_par *par, unsigned splitmode )
     /*----------------------------------------------------------------------*/
     /* set setup variant for adaptive rejection sampling                    */
     /*                                                                      */
     /* There are three variants for adaptive rejection sampling. These      */
     /* differ in the way how an interval is split:                          */
     /*    splitmode 1: use the generated point to split the interval.       */
     /*    splitmode 2: use the mean point of the interval.                  */
     /*    splitmode 3: use the arcmean point.                               */
     /* Default is splitmode 3.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   splitmode ... indicator for variant                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* store date */
  par->variant &= ~TABL_VARMASK_SPLIT;
  switch (splitmode) {
  case 1:
    par->variant |= TABL_VARFLAG_SPLIT_POINT;
    return UNUR_SUCCESS;
  case 2:
    par->variant |= TABL_VARFLAG_SPLIT_MEAN;
    return UNUR_SUCCESS;
  case 3:
    par->variant |= TABL_VARFLAG_SPLIT_ARC;
    return UNUR_SUCCESS;
  default:
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid variant");
    return UNUR_ERR_PAR_SET;
  }
} /* end if unur_tabl_set_variant_splitmode() */

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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ratio = max_ratio;

  /* changelog */
  par->set |= TABL_SET_MAX_SQHRATIO;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_sqhratio( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio    ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TABL, INFINITY );

  return (GEN->Asqueeze / GEN->Atotal);

} /* end of unur_tabl_get_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_hatarea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below hat                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TABL, INFINITY );

  return GEN->Atotal;

} /* end of unur_tabl_get_hatarea() */

/*---------------------------------------------------------------------------*/

double
unur_tabl_get_squeezearea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below squeeze                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, TABL, INFINITY );

  return GEN->Asqueeze;

} /* end of unur_tabl_get_squeezearea() */

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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= TABL_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_get_n_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get current number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   number of intervals     ... on success                             */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, TABL, 0 );

  return GEN->n_ivs;

} /* end of unur_tabl_get_n_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_areafraction( struct unur_par *par, double fraction )
     /*----------------------------------------------------------------------*/
     /* set parameter for equal area rule                                    */
     /* (each bar has size fraction * area below PDF)                        */           
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   fraction  ... fraction of area for bar                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (fraction <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"area factor <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->area_fract = fraction;

  /* changelog */
  par->set |= TABL_SET_AREAFRACTION;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_areafraction() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_cpoints( struct unur_par *par, int n_stp, const double *stp )
     /*----------------------------------------------------------------------*/
     /* set construction points for hat function                             */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*   stp    ... pointer to array of starting points                     */
     /*              (NULL for changing only the number of default points)   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp <= 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points <= 0");
    return UNUR_ERR_PAR_SET;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }

  /* store date */
  PAR->starting_cpoints = stp;
  PAR->n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TABL_SET_N_STP | ((stp) ? TABL_SET_STP : 0u);

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_nstp( struct unur_par *par, int n_stp )
     /*----------------------------------------------------------------------*/
     /* OBSOLETE !!                                                          */
     /*                                                                      */
     /* set number of construction points for hat at initialization          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= TABL_SET_N_STP;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_nstp() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_set_slopes( struct unur_par *par, const double *slopes, int n_slopes )
     /*----------------------------------------------------------------------*/
     /* set slopes of PDF                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   slopes   ... pointer to list of slopes                             */
     /*   n_slopes ... number of slopes                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a]                        */
     /*   such that PDF(a) >= PDF(b).                                        */
     /*   slopes must be decreasing, non-overlapping and sorted              */
     /*----------------------------------------------------------------------*/
{
  int i;
  double lmin,lmax,rmin,rmax;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if( n_slopes <= 0 ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"number of slopes <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* check slopes */
  lmin = -INFINITY;
  lmax = -INFINITY;
  for( i=0; i<n_slopes; i++ ) {
    rmin = min(slopes[2*i],slopes[2*i+1]);
    rmax = max(slopes[2*i],slopes[2*i+1]);
    if (!(lmax<=rmin || _unur_FP_same(lmax,rmin))) {
      _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes (overlapping or not in ascending order)");
      return UNUR_ERR_PAR_SET;
    }
    lmin = rmin;
    lmax = rmax;
  }

  /* INFINITY is not allowed */
  if (! (_unur_isfinite(slopes[0]) && _unur_isfinite(slopes[2*n_slopes-1])) ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"slopes must be bounded");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->slopes = slopes;
  PAR->n_slopes = n_slopes;

  /* changelog */
  par->set |= TABL_SET_SLOPES;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_slopes() */

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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= TABL_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_guidefactor() */

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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   new boundary points must not be +/- INFINITY                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (left <= -INFINITY || right >= INFINITY) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->bleft = left;
  PAR->bright = right;

  /* changelog */
  par->set |= TABL_SET_BOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_tabl_set_boundary() */

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
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, TABL );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | TABL_VARFLAG_VERIFY) : (par->variant & (~TABL_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_tabl_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, TABL, UNUR_ERR_GEN_INVALID );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= TABL_VARFLAG_VERIFY;
    SAMPLE = _unur_tabl_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~TABL_VARFLAG_VERIFY;
    SAMPLE = _unur_tabl_sample;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_tabl_chg_verify() */

/*---------------------------------------------------------------------------*/
