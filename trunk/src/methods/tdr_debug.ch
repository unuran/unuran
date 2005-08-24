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
_unur_tdr_debug_init( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TDR_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = transformed density rejection\n",gen->genid);
  fprintf(log,"%s: variant = ",gen->genid);
  switch (par->variant & TDR_VARMASK_VARIANT) {
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
    fprintf(log,"-x^(%g)  ... c = %g",PAR->c_T,PAR->c_T); break;
  }
  _unur_print_if_default(par,TDR_SET_C);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tdr_",gen->genid);
  switch (par->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    fprintf(log,"gw"); break;
  case TDR_VARIANT_PS:
    fprintf(log,"ps"); break;
  case TDR_VARIANT_IA:
    fprintf(log,"ia"); break;
  }
  if (par->variant & TDR_VARFLAG_VERIFY)
    fprintf(log,"_sample_check()\n");
  else
    fprintf(log,"_sample()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: center = %g",gen->genid,PAR->center);
  _unur_print_if_default(par,TDR_SET_CENTER);
  if (par->variant & TDR_VARFLAG_USEMODE)
    fprintf(log,"\n%s: use mode as construction point",gen->genid);
  else if (par->variant & TDR_VARFLAG_USECENTER)
    fprintf(log,"\n%s: use center as construction point",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(par,TDR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,PAR->max_ratio*100.);
  _unur_print_if_default(par,TDR_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  if (par->variant & TDR_VARFLAG_USEDARS) {
    fprintf(log,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(par,TDR_SET_USE_DARS);
    fprintf(log,"\n%s:\tDARS factor = %g",gen->genid,PAR->darsfactor);
    _unur_print_if_default(par,TDR_SET_DARS_FACTOR);
    fprintf(log,"\n%s:\tDARS rule %d",gen->genid,PAR->darsrule);
    _unur_print_if_default(par,TDR_SET_USE_DARS);
  }
  else {
    fprintf(log,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(par,TDR_SET_USE_DARS);
  }
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR->guide_factor);
  _unur_print_if_default(par,TDR_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(par,TDR_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & TDR_SET_STP)
    for (i=0; i<PAR->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",PAR->starting_cpoints[i]);
    }
  else
    fprintf(log," use \"equdistribution\" rule [default]");
  fprintf(log,"\n%s:\n",gen->genid);
  
  _unur_tdr_debug_intervals(gen);

  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_init() */

/*****************************************************************************/

void 
_unur_tdr_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print header before runniung DARS into logfile                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TDR_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: DARS started **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS factor = %g",gen->genid,PAR->darsfactor);
  _unur_print_if_default(par,TDR_SET_DARS_FACTOR);
  fprintf(log,"\n%s: DARS rule %d",gen->genid,PAR->darsrule);
  _unur_print_if_default(par,TDR_SET_USE_DARS);
  fprintf(log,"\n%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tdr_debug_dars_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tdr_debug_dars( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print infor after generator has run DARS into logfile                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TDR_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS finished **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_tdr_debug_intervals(gen);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tdr_debug_dars() */

/*****************************************************************************/

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
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_tdr_debug_intervals(gen);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_free() */

/*****************************************************************************/

void
_unur_tdr_debug_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile (orig. variant by Gilks & Wild) */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    _unur_tdr_gw_debug_intervals(gen);
    return;
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    _unur_tdr_ps_debug_intervals(gen);
    return;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return;
  }
} /* end of _unur_tdr_debug_intervals() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/



