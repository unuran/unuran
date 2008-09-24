/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_newset.c                                                *
 *                                                                           *
 *   Routines for creating and changing parameter objects.                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_pinv_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pdf == NULL && DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF or CDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_pinv_par) );
  COOKIE_SET(par,CK_PINV_PAR);

  /* copy input */
  par->distr   = distr;           /* pointer to distribution object          */

  /* set default values */
  PAR->order = 5;                /* order of polynomial                      */
  PAR->u_resolution = 1.0e-10;   /* maximal error allowed in u-direction     */
  PAR->bleft = -1.e100;          /* left border of the computational domain  */
  PAR->bright = 1.e100;          /* right border of the computational domain */
  PAR->sleft = TRUE;             /* whether to search for left boundary      */
  PAR->sright = TRUE;            /* whether to search for right boundary     */
  PAR->max_ivs = PINV_DEFAULT_MAX_IVS; /* maximum number of subintervals     */

  par->method   = UNUR_METH_PINV; /* method                                  */
  par->variant  = ( (DISTR_IN.pdf != NULL)  /* default variant:              */
		    ? PINV_VARIANT_PDF      /*   use PDF                     */
		    : PINV_VARIANT_CDF );   /*   use CDF                     */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_pinv_init;

  return par;

} /* end of unur_pinv_new() */

/*****************************************************************************/

int
unur_pinv_set_order( struct unur_par *par, int order)
     /*----------------------------------------------------------------------*/
     /* Set order of Hermite interpolation.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   order  ... order of interpolation polynome                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (order<2 || order>MAX_ORDER) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"order <2 or >19");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->order = order;

  /* changelog */
  par->set |= PINV_SET_ORDER;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_order() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... maximal error in u                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (u_resolution > 1.e-2) {
    /* this is obviously an error */
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too large --> not changed");
    return UNUR_ERR_PAR_SET;
  }
  if (u_resolution < 5.*DBL_EPSILON ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution too small --> corrected");
    u_resolution = 5.*DBL_EPSILON;
  }

  /* store date */
  PAR->u_resolution = u_resolution;

  /* changelog */
  par->set |= PINV_SET_U_RESOLUTION;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_usepdf( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use PDF (if available) to compute approximate inverse CDF.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check function pointer */
  if (par->distr->data.cont.pdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"PDF missing");
    return UNUR_ERR_PAR_SET;
  }

  /* store variant */
  par->variant = PINV_VARIANT_PDF;

  /* changelog */
  par->set |= PINV_SET_VARIANT;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_pinv_set_usepdf() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_usecdf( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Use CDF (if available) to compute approximate inverse CDF.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* check function pointer */
  if (par->distr->data.cont.cdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF missing");
    return UNUR_ERR_PAR_SET;
  }

  /* store variant */
  par->variant = PINV_VARIANT_CDF;

  /* changelog */
  par->set |= PINV_SET_VARIANT;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_pinv_set_usecdf() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_boundary( struct unur_par *par, double left, double right )
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
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (!_unur_FP_less(left,right)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain");
    return UNUR_ERR_PAR_SET;
  }
  if (! (_unur_isfinite(left) && _unur_isfinite(right)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"domain (+/- INFINITY not allowed)");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->bleft = left;
  PAR->bright = right;

  /* changelog */
  par->set |= PINV_SET_BOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_boundary() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_searchboundary( struct unur_par *par, int left, int right )
     /*----------------------------------------------------------------------*/
     /* set flag for boundary searching algorithm                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   left  ... whether to search for left boundary point                */
     /*   right ... whether to search for right boundary point               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, PINV );

  /* store date */
  PAR->sleft  = (left)  ? TRUE : FALSE;
  PAR->sright = (right) ? TRUE : FALSE;

  /* changelog */
  par->set |= PINV_SET_SEARCHBOUNDARY;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_searchboundary() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_set_max_intervals( struct unur_par *par, int max_ivs )
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
  _unur_check_par_object( par, PINV );

  /* check new parameter for generator */
  if (max_ivs < 100 || max_ivs > 1000000) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 100 or > 1000000");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= PINV_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_pinv_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_pinv_get_n_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get number of intervals (or more precisely the number of nodes)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   number of intervals ... on success                                 */
     /*   0     ... on error                                                 */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, 0 );
  _unur_check_gen_object( gen, PINV, 0 );
  return GEN->n_ivs;
} /* end of unur_pinv_get_n_intervals() */

/*****************************************************************************/
