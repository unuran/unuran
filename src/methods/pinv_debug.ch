/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_debug.c                                                 *
 *                                                                           *
 *   Routines for printing debugging information.                            * 
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile before setup starts.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = PINV (Polynomial interpolation based INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_pinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: order of polynomial = %d",gen->genid,GEN->order);
  _unur_print_if_default(gen,PINV_SET_ORDER);
  fprintf(log,"\n%s: u-resolution = %g",gen->genid,GEN->u_resolution);
  _unur_print_if_default(gen,PINV_SET_U_RESOLUTION);
  fprintf(log,"\n%s: variant = ",gen->genid);
  switch (gen->variant) {
  case PINV_VARIANT_PDF:
    fprintf(log,"use PDF + Lobatto integration"); break;
  case PINV_VARIANT_CDF:
    fprintf(log,"use CDF"); break;
  }
  _unur_print_if_default(gen,PINV_SET_VARIANT);
  fprintf(log,"\n");

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} /* end of _unur_pinv_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_init( const struct unur_gen *gen, int ok )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   ok  ... exitcode of init call                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: domain of computation = [%g,%g]\n",gen->genid, GEN->bleft,GEN->bright);
  fprintf(log,"%s: Umin = 0 [fixed], Umax = %18.16g",gen->genid, GEN->Umax);
  if (_unur_FP_approx(GEN->Umax,1))
    fprintf(log,",  1-Umax = %g",1.-GEN->Umax);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: # Intervals = %d\n",gen->genid,GEN->n_ivs);
  fprintf(log,"%s:\n",gen->genid);

  _unur_pinv_debug_intervals(gen);

  fprintf(log,"%s: initialization %s\n",gen->genid,((ok)?"successful":"*** FAILED ***")); 
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_pinv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_relevant_support (const struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* print relevant domain                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: search for boundaries: left = %s, right = %s\n",gen->genid,
	  GEN->sleft ? "TRUE" : "FALSE", GEN->sright ? "TRUE" : "FALSE");
  fprintf(log,"%s: relevant domain = (%g,%g)   [i.e. where PDF > threshold]\n",gen->genid,
	  GEN->bleft,GEN->bright);
  fprintf(log,"%s: possible support of distribution = (%g,%g)\n",gen->genid,
	  GEN->dleft,GEN->dright);

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);

} /* end of _unur_pinv_debug_relevant_support() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_pdfarea (const struct unur_gen *gen, int approx)
     /*----------------------------------------------------------------------*/
     /* print estimated area below PDF                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   approx ... whether this is an approximate value or accurat up to   */
     /*              requested tolerance                                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: area below PDF %s = %19.16g\n",gen->genid,
	  approx ? "(approx.)" : "(accurate)", GEN->area);

#ifdef PINV_USE_CDFTABLE
  if (GEN->CDFtable && (GEN->CDFtable->n_values > 0)) {
    fprintf(log,"%s: # subintervals in Lobatto integration = %d\n",gen->genid,
	    GEN->CDFtable->n_values-1);
  }
#endif

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} /* end of _unur_pinv_debug_pdfarea() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_computational_domain (const struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* print computational domain                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: computational domain = (%g,%g)\n",gen->genid,
	    GEN->bleft,GEN->bright);

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);

} /* end of _unur_pinv_debug_computational_domain() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print table of intervals into logfile                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int n;
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (gen->debug & PINV_DEBUG_TABLE) {
    for (n=0; n<=GEN->n_ivs; n++) {
      fprintf(log,"%s: [%3d] xi = %g, cdfi = %g\n",gen->genid,
	      n, GEN->iv[n].xi, GEN->iv[n].cdfi);
    }
  }

  fprintf(log,"%s:\n",gen->genid);
  fflush(log);

} /* end of _unur_pinv_debug_intervals() */

/*---------------------------------------------------------------------------*/

void
_unur_pinv_debug_create_table (const struct unur_gen *gen,
			       int iter, int n_incr_h, int n_decr_h)
     /*----------------------------------------------------------------------*/
     /* print data that have been collected while creating polynomials.      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iter     ... total number of iterations                            */
     /*   n_incr_h ... number of steps where h is increased                  */
     /*   n_decr_h ... number of steps where h is decreased                  */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_PINV_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: Create interpolating polynomaials:\n",gen->genid);
  fprintf(log,"%s:\t# iterations   = %d\n",gen->genid,iter);
  fprintf(log,"%s:\t# increasing h = %d  (%g%%)\n",gen->genid,n_incr_h,(100.*n_incr_h)/iter);
  fprintf(log,"%s:\t# decreasing h = %d  (%g%%)\n",gen->genid,n_decr_h,(100.*n_decr_h)/iter);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_pinv_debug_create_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
