/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_debug.c                                                  *
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

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#include "tdr_gw_debug.ch"
#include "tdr_ps_debug.ch"

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print after (almost empty generator) object has been created.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = transformed density rejection\n",gen->genid);
  fprintf(log,"%s: variant = ",gen->genid);
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    fprintf(log,"original (Gilks & Wild)  ... GW\n"); break;
  case TDR_VARIANT_PS:
    fprintf(log,"proportional squeeze  ... PS\n"); break;
  case TDR_VARIANT_IA:
    fprintf(log,"immediate acceptance  ... IA\n"); break;
  }
  fprintf(log,"%s: transformation T_c(x) = ",gen->genid);
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    fprintf(log,"log(x)  ... c = 0");                   break;
  case TDR_VAR_T_SQRT:
    fprintf(log,"-1/sqrt(x)  ... c = -1/2");            break;
  case TDR_VAR_T_POW:
    fprintf(log,"-x^(%g)  ... c = %g",GEN->c_T,GEN->c_T); break;
  }
  _unur_print_if_default(gen,TDR_SET_C);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tdr_",gen->genid);
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    fprintf(log,"gw"); break;
  case TDR_VARIANT_PS:
    fprintf(log,"ps"); break;
  case TDR_VARIANT_IA:
    fprintf(log,"ia"); break;
  }
  if (gen->variant & TDR_VARFLAG_VERIFY)
    fprintf(log,"_sample_check()\n");
  else
    fprintf(log,"_sample()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: center = %g",gen->genid,GEN->center);
  _unur_print_if_default(gen,TDR_SET_CENTER);
  if (gen->variant & TDR_VARFLAG_USEMODE)
    fprintf(log,"\n%s: use mode as construction point",gen->genid);
  else if (gen->variant & TDR_VARFLAG_USECENTER)
    fprintf(log,"\n%s: use center as construction point",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,TDR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,GEN->max_ratio*100.);
  _unur_print_if_default(gen,TDR_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  if (gen->variant & TDR_VARFLAG_USEDARS) {
    fprintf(log,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(gen,TDR_SET_USE_DARS);
    fprintf(log,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(gen,TDR_SET_DARS_FACTOR);
    fprintf(log,"\n%s:\tDARS rule %d",gen->genid,GEN->darsrule);
    _unur_print_if_default(gen,TDR_SET_USE_DARS);
  }
  else {
    fprintf(log,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(gen,TDR_SET_USE_DARS);
  }
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,TDR_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,GEN->n_starting_cpoints);
  _unur_print_if_default(gen,TDR_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (gen->set & TDR_SET_STP)
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",GEN->starting_cpoints[i]);
    }
  else
    fprintf(log," use \"equidistribution\" rule [default]");
  fprintf(log,"\n%s:\n",gen->genid);
  
  fflush(log);

} /* end of _unur_tdr_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_tdr_debug_intervals(gen,"INIT completed",TRUE);

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void 
_unur_tdr_debug_dars_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print header before runniung DARS into logfile                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (gen->debug & TDR_DEBUG_IV) 
    _unur_tdr_debug_intervals(gen,"Starting Intervals:",TRUE);
  fprintf(log,"%s: DARS started **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(gen,TDR_SET_DARS_FACTOR);
  fprintf(log,"\n%s: DARS rule %d",gen->genid,GEN->darsrule);
  _unur_print_if_default(gen,TDR_SET_USE_DARS);
  fprintf(log,"\n%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tdr_debug_dars_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_dars_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print infor after generator has run DARS into logfile                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS finished **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tdr_debug_dars_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_reinit_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before reinitialization into logfile      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int i;
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: *** Re-Initialize generator object ***\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  if (gen->set & TDR_SET_N_PERCENTILES) {
    fprintf(log,"%s: use percentiles of old hat as starting points for new hat:",gen->genid);
    for (i=0; i<GEN->n_percentiles; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",GEN->percentiles[i]);
    }
    fprintf(log,"\n%s: starting points:",gen->genid);
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",GEN->starting_cpoints[i]);
    }
    fprintf(log,"\n");
  }
  else {
    fprintf(log,"%s: use starting points given at init\n",gen->genid);
  }
  
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_reinit_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_reinit_retry( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before second trial of reinitialization   */
     /* into logfile                                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: *** Re-Initialize failed  -->  second trial ***\n",gen->genid);
  fprintf(log,"%s: use equal-area-rule with %d points\n",gen->genid,GEN->retry_ncpoints);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_reinit_retry() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_reinit_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after reinitialization into logfile       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_tdr_debug_intervals(gen,"*** Generator reinitialized ***",TRUE);

  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_reinit_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    _unur_tdr_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(log,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
    _unur_tdr_debug_intervals(gen,"Intervals after failure:",FALSE);
  }
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_free() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile (orig. variant by Gilks & Wild) */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen         ... pointer to generator object                        */
     /*   header      ... header for table                                   */
     /*   print_areas ... whether table of areas should be printed           */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  log = unur_get_stream();

  if (header) fprintf(log,"%s:%s\n",gen->genid,header);

  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    _unur_tdr_gw_debug_intervals(gen,print_areas);
    return;
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    _unur_tdr_ps_debug_intervals(gen,print_areas);
    return;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return;
  }
} /* end of _unur_tdr_debug_intervals() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/



