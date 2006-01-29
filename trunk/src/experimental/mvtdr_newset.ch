
/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mvtdr_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  /* first we check parameters */
  
  /* dim (=N) must be an integer >= 2 */
/*   if( (int) N != N || N < 2 || N > 15 ) */
/*     FATAL( "N in config.h not an integer >= 2 and <= 15" ); */


/*   if (!(distr->set & UNUR_DISTR_SET_MEAN)) { */
/*     _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mean"); return NULL; } */
/*   if (!(distr->set & UNUR_DISTR_SET_COVAR)) { */
/*     _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"covariance matrix"); */
/*     return NULL; } */
/*   if (!(distr->set & UNUR_DISTR_SET_STDMARGINAL)) { */
/*     _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"standardized marginals"); */
/*     return NULL; } */

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mvtdr_par) );
  COOKIE_SET(par,CK_MVTDR_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_MVTDR ;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_mvtdr_init;


  /*
    get default values for construction of hat function 
  */
  /** set default values **/
  /* the mode */
#if MODE == 1
#endif

  /* control generation of cones */
  PAR->max_cones = MAX_N_CONES;                /* maximum number of cones (at least 2^(N+T_STEPS_MIN) */
  PAR->steps_min = T_STEPS_MIN;                /* minimum number of triangulation steps */
  PAR->step_tp = OPTIMAL_TP_STEP;              /* triangulation step when optimal touching points is calculated */

#if MODE == 1
  /* parameters for finding mode */
  PAR->mode_to_boundary = MODE_TO_BOUNDARY;    /* move mode to boundary if |mode - boundary| / length < MODE_TO_BOUNDARY */
#endif

#if RECTANGLE == 1
#endif

  return par;

} /* end of unur_mvtdr_new() */

/*---------------------------------------------------------------------------*/
