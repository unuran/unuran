
/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = MVTDR (Multi-Variate Transformed Density Rejection)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_mvtdr_sample_cvec()\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_mvtdr_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_mvtdr_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVTDR_GEN,RETURN_VOID);

  log = unur_get_stream();

/*   /\* print data *\/ */
/*   fprintf(LOG,"* Setup *\n\n"); */
/*   fprintf(LOG,"time = %g sec\n",time_stop - time_start); */
/*   fprintf(LOG,"Volume: hat =     %f\n",cr->Htot); */
/*   if( cr->volume_density > 0. ) { */
/*     fprintf(LOG,"\tdensity = %f\n",cr->volume_density); */
/*     fprintf(LOG,"\tratio =   %f %%\n",( cr->volume_density / cr->Htot )*100); */
/*   } */

/*   /\* triangulation steps *\/ */
/*   fprintf(LOG,"minimum triangulation level = %d\n",cr->steps_min); */
/*   fprintf(LOG,"maximum triangulation level = %d\n",cr->steps_max); */
/*   fprintf(LOG,"highest level for finding touching points = %d\n",cr->step_tp); */



  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: INIT completed **********************\n",gen->genid);

} /* end of _unur_mvtdr_debug_init_finished() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
