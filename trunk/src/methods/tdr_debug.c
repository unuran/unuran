/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and .... of a T-concave distribution                     *
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
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_tdr_debug_init( struct unur_par *par, struct unur_gen *gen )
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_TDR_PAR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = transformed density rejection\n",gen->genid);
  fprintf(log,"%s: version = ",gen->genid);
  switch (par->variant & TDR_VARMASK_VERSION) {
  case TDR_VAR_VERSION_ORIG:
    fprintf(log,"orignal (Gilks & Wild)  ... ORIG\n"); break;
  case TDR_VAR_VERSION_PS:
    fprintf(log,"proportional squeeze  ... PS\n"); break;
  case TDR_VAR_VERSION_IA:
    fprintf(log,"immeditaly acceptance  ... IA\n"); break;
  }
  fprintf(log,"%s: transformation T_c(x) = ",gen->genid);
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    fprintf(log,"log(x)  ... c = 0");                   break;
  case TDR_VAR_T_SQRT:
    fprintf(log,"-1/sqrt(x)  ... c = -1/2");            break;
  case TDR_VAR_T_POW:
    fprintf(log,"-x^(%g)  ... c = %g",PAR.c_T,PAR.c_T); break;
  }
  _unur_print_if_default(par,TDR_SET_C);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tdr_sample_",gen->genid);
  switch (par->variant & TDR_VARMASK_VERSION) {
  case TDR_VAR_VERSION_ORIG:
    fprintf(log,"orig"); break;
  case TDR_VAR_VERSION_PS:
    fprintf(log,"ps"); break;
  case TDR_VAR_VERSION_IA:
    fprintf(log,"ia"); break;
  }
  if (par->variant & TDR_VARFLAG_VERIFY)
    fprintf(log,"check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: center = %g",gen->genid,PAR.center);
  _unur_print_if_default(par,TDR_SET_CENTER);
  if (par->variant & TDR_VARFLAG_USEMODE)
    fprintf(log,"\n%s: use mode as construction point",gen->genid);
  else if (par->variant & TDR_VARFLAG_USECENTER)
    fprintf(log,"\n%s: use center as construction point",gen->genid);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR.max_ivs);
  _unur_print_if_default(par,TDR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,TDR_SET_MAX_SQHRATIO);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
  _unur_print_if_default(par,TDR_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting points = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,TDR_SET_N_STP);
  fprintf(log,"\n%s: starting points:",gen->genid);
  if (par->set & TDR_SET_STP)
    for (i=0; i<PAR.n_starting_cpoints; i++) {
      if (i%5==0) fprintf(log,"\n%s:\t",gen->genid);
      fprintf(log,"   %#g,",PAR.starting_cpoints[i]);
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

static void
_unur_tdr_debug_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  _unur_tdr_debug_intervals(gen);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_free() */

/*****************************************************************************/

static void
_unur_tdr_debug_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tdr_interval *iv;
  double sAsqueeze, sAhatl, sAhatr;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:Intervals: %d\n",gen->genid,GEN.n_ivs);
  if (GEN.iv) {
    if (gen->debug & TDR_DEBUG_IV) {
      fprintf(log,"%s: Nr.            tp          f(tp)        T(f(tp))    d(T(f(tp)))        squeeze\n",gen->genid);
      for (iv = GEN.iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
	fprintf(log,"%s:[%3d]: %#12.6g   %12.6g   %#12.6g   %#12.6g   %12.6g\n", gen->genid, i,
		iv->x, iv->fx, iv->Tfx, iv->dTfx, iv->sq);
      }
      COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
      fprintf(log,"%s:[...]: %#12.6g   %12.6g   %#12.6g   %#12.6g\n", gen->genid,
	      iv->x, iv->fx, iv->Tfx, iv->dTfx);
    }
    fprintf(log,"%s:\n",gen->genid);
  }
  else
    fprintf(log,"%s: No intervals !\n",gen->genid);

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TDR_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\tbelow squeeze\t\t  below hat (left and right)\t\t  cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
    if (GEN.iv) {
      for (iv = GEN.iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,/*void*/); 
	sAsqueeze += iv->Asqueeze;
	sAhatl += iv->Ahatl;
	sAhatr += iv->Ahatr;
	fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		iv->Asqueeze, iv->Asqueeze * 100. / GEN.Atotal,
		iv->Ahatl, iv->Ahatr, (iv->Ahatl+iv->Ahatr) * 100. / GEN.Atotal, 
		iv->Acum, iv->Acum * 100. / GEN.Atotal);
      }
      fprintf(log,"%s:       ----------  ---------  |  ------------------------  ---------  +\n",gen->genid);
      fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)            %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAsqueeze, sAsqueeze * 100. / GEN.Atotal,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / GEN.Atotal);
      fprintf(log,"%s:\n",gen->genid);
    }
  }

  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_tdr_debug_intervals() */

/*****************************************************************************/

static void
_unur_tdr_debug_sample( struct unur_gen *gen, 
			struct unur_tdr_interval *iv, 
			struct unur_tdr_interval *pt, 
			double x, double fx, double hx, double sqx )
     /*----------------------------------------------------------------------*/
     /* write info about generated point                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   pt  ... pointer to interval that stores construction point         */
     /*   x   ... generated point                                            */
     /*   fx  ... value of p.d.f. at x                                       */
     /*   hx  ... value of hat at x                                          */
     /*   sqx ... value of squeeze at x                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv,/*void*/);   COOKIE_CHECK(iv,CK_TDR_IV,/*void*/);
  CHECK_NULL(pt,/*void*/);   COOKIE_CHECK(pt,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  if (iv == pt)
    fprintf(log,"%s: point generated in left part:\n",gen->genid);
  else
    fprintf(log,"%s: point generated in right part:\n",gen->genid);

  fprintf(log,"%s: construction point: x0 = %g\n",gen->genid,pt->x);
  fprintf(log,"%s: transformed hat     Th(x) = %g + %g * (x - %g)\n",gen->genid,pt->Tfx,pt->dTfx,pt->x);
  fprintf(log,"%s: transformed squeeze Ts(x) = %g + %g * (x - %g)\n",gen->genid,iv->Tfx,iv->sq,iv->x);
  fprintf(log,"%s: generated point: x = %g\n",gen->genid,x);
  fprintf(log,"%s:  h(x) = %.20g\n",gen->genid,hx);
  fprintf(log,"%s:  f(x) = %.20g\n",gen->genid,fx);
  fprintf(log,"%s:  s(x) = %.20g\n",gen->genid,sqx);
  fprintf(log,"%s:    hat: x - x0 = %g",gen->genid,x-pt->x);
  if (x < pt->x && iv == pt) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    h(x) - f(x) = %g",gen->genid,hx-fx);
  if (hx<fx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    squeeze: x - x0 = %g",gen->genid,x-iv->x);
  if (x > pt->x && iv != pt) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    f(x) - s(x) = %g",gen->genid,fx-sqx);
  if (fx<sqx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_sample() */

/*****************************************************************************/

static void
_unur_tdr_debug_split_start( struct unur_gen *gen, struct unur_tdr_interval *iv, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting interval                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   x   ... split at this point                                        */
     /*   fx  ... value of p.d.f. at x                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv,/*void*/);   COOKIE_CHECK(iv,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: split interval at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(log,"%s: old interval:\n",gen->genid);
  fprintf(log,"%s:   left  construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->x,iv->fx);
  fprintf(log,"%s:   right construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->fx);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv->Asqueeze,iv->Asqueeze*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv->Ahatl+iv->Ahatr-iv->Asqueeze),(iv->Ahatl+iv->Ahatr-iv->Asqueeze)*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv->Ahatl,iv->Ahatr,(iv->Ahatl+iv->Ahatr)*100./GEN.Atotal);

  fflush(log);

} /* end of _unur_tdr_debug_split_start() */

/*****************************************************************************/

static void
_unur_tdr_debug_split_stop( struct unur_gen *gen, 
			    struct unur_tdr_interval *iv_left, 
			    struct unur_tdr_interval *iv_right )
     /*----------------------------------------------------------------------*/
     /* write info about new splitted intervals                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv_left  ... pointer to new left hand interval                     */
     /*   iv_right ... pointer to new right hand interval                    */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);       COOKIE_CHECK(gen,CK_TDR_GEN,/*void*/);
  CHECK_NULL(iv_left,/*void*/);   COOKIE_CHECK(iv_left,CK_TDR_IV,/*void*/);
  CHECK_NULL(iv_right,/*void*/);  COOKIE_CHECK(iv_right,CK_TDR_IV,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: inserted point:\n",gen->genid);
  fprintf(log,"%s: x = %g, f(x) = %g, Tf(x)=%g, dTf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->fx, iv_right->Tfx, iv_right->dTfx, iv_right->sq);
  fprintf(log,"%s: new intervals:\n",gen->genid);
  fprintf(log,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(log,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(log,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);

  fprintf(log,"%s: left interval:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv_left->Asqueeze,
	  iv_left->Asqueeze*100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv_left->Ahatl + iv_left->Ahatr - iv_left->Asqueeze),
	  (iv_left->Ahatl + iv_left->Ahatr - iv_left->Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv_left->Ahatl,
	  iv_left->Ahatr,
	  (iv_left->Ahatl + iv_left->Ahatr) * 100./GEN.Atotal);

  if (iv_left == iv_right)
    fprintf(log,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(log,"%s: right interval:\n",gen->genid);
    fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,
	    iv_right->Asqueeze*100./GEN.Atotal);
    fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    (iv_right->Ahatl + iv_right->Ahatr - iv_right->Asqueeze),
	    (iv_right->Ahatl + iv_right->Ahatr - iv_right->Asqueeze) * 100./GEN.Atotal);
    fprintf(log,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	    iv_right->Ahatl,
	    iv_right->Ahatr,
	    (iv_right->Ahatl + iv_right->Ahatr) * 100./GEN.Atotal);
  }

  fprintf(log,"%s: total areas:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  fprintf(log,"%s:\n",gen->genid);

  fflush(log);

} /* end of _unur_tdr_debug_split_stop() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
