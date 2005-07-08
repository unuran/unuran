/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_debug.c                                                 *
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
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(log,"%s:\n",gen->genid);

/*---------------------------------------------------------------------------*/

void
_unur_tabl_debug_init( const struct unur_par *par, const struct unur_gen *gen )
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
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  empty_line();

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_tabl_sample",gen->genid);
  if (par->variant & TABL_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  empty_line();

  fprintf(log,"%s: computation interval = (%g, %g)\n",gen->genid,GEN->bleft,GEN->bright);
  empty_line();

  fprintf(log,"%s: area fraction for equal area rule = %g ",gen->genid,PAR->area_fract);
  _unur_print_if_default(par,TABL_SET_AREAFRACTION);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR->max_ivs);
  _unur_print_if_default(par,TABL_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR->max_ratio*100.);
  _unur_print_if_default(par,TABL_SET_MAX_SQHRATIO);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR->guide_factor);
  _unur_print_if_default(par,TABL_SET_GUIDEFACTOR);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: split intervals at ",gen->genid);
  switch( gen->variant & TABL_VARMASK_SPLIT ) {
  case TABL_VARFLAG_SPLIT_MEAN:
    fprintf(log,"mean point");
    break;
  case TABL_VARFLAG_SPLIT_ARC:
    fprintf(log,"\"arcmean\" point");
    break;
  case TABL_VARFLAG_SPLIT_POINT:
  default: 
    fprintf(log,"sample point");
    break;
  }
  fprintf(log," when using adaptive sampling.\n");
  empty_line();

  if (par->set & TABL_SET_SLOPES) {
    fprintf(log,"%s: slopes = %d\n",gen->genid,PAR->n_slopes);
    for (i=0; i<PAR->n_slopes; i++) {
      if ( PAR->slopes[2*i] > PAR->slopes[2*i+1] )
	fprintf(log,"%s:   (+)  ",gen->genid);
      else
	fprintf(log,"%s:   (-)  ",gen->genid);
      fprintf(log,"< %#g, %#g >\n", PAR->slopes[2*i], PAR->slopes[2*i+1] );
    }
  }
  else
    fprintf(log,"%s: no slopes given. compute from domain and mode.\n",gen->genid);

  if (par->variant & TABL_VARFLAG_STP_A)
    fprintf(log,"%s: split slopes by equal area rule (SPLIT A).\n",gen->genid);

  if (par->variant & TABL_VARFLAG_USEDARS) {
    fprintf(log,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(par,TABL_SET_USE_DARS);
    fprintf(log,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  }
  else {
    fprintf(log,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(par,TABL_SET_USE_DARS);
  }
  fprintf(log,"\n%s:\n",gen->genid);

  fprintf(log,"%s: number of starting intervals (approx.) = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(par,TABL_SET_N_STP);
  fprintf(log,"\n");
  empty_line();

  if (gen->debug & TABL_DEBUG_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);

} /* end of _unur_tabl_debug_init() */

/*****************************************************************************/

void
_unur_tabl_debug_init_finished( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  _unur_tabl_debug_intervals(gen,TRUE);

  empty_line();
  fprintf(log,"%s:  INIT completed **********************\n",gen->genid);
  empty_line();

} /* end of _unur_tabl_debug_init_finished() */

/*****************************************************************************/

void 
_unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen )
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
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: DARS started **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tabl_debug_dars_start() */

/*---------------------------------------------------------------------------*/

void
_unur_tabl_debug_dars( const struct unur_par *par, const struct unur_gen *gen )
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
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS finished **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  if (gen->debug & TABL_DEBUG_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS completed **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fflush(log);
} /* end of _unur_tabl_debug_dars() */

/*****************************************************************************/

void
_unur_tabl_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  empty_line();
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  empty_line();
  _unur_tabl_debug_intervals(gen,TRUE);
  empty_line();

  fflush(log);

} /* end of _unur_tabl_debug_free() */

/*****************************************************************************/

void
_unur_tabl_debug_intervals( const struct unur_gen *gen, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tabl_interval *iv;
  double sAsqueeze, Atotal;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: intervals = %d\n",gen->genid,GEN->n_ivs);
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID);
      fprintf(log,"%s:[%3d]: (",gen->genid,i);
      switch (iv->slope) {
      case 1:  fprintf(log,"+"); break;
      case 0:  fprintf(log,"0"); break;
      case -1: fprintf(log,"-"); break;
      default:
	/* this should not happen:
	   invalid value for iv->slope. */
	fprintf(log,"?"); break;
      }
      fprintf(log,")   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }

  if (!print_areas) {
    /* this just happens, when we call _unur_tabl_debug_intervals(), after
       we have finished split A. */
    fprintf(log,"%s: Split A finished.\n",gen->genid);
    empty_line();
    return;
  }

  if (GEN->Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    Atotal = -1.;   /* to avoid floating point exceptions */
  }
  else {
    Atotal = GEN->Atotal;
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID); 
      sAsqueeze += iv->Asqueeze;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
	      iv->Ahat, iv->Ahat * 100. / Atotal, 
	      iv->Acum, iv->Acum * 100. / Atotal);
    }
    fprintf(log,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / Atotal, Atotal);
    empty_line();
  }
    
  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);

  empty_line();

} /* end of _unur_tabl_debug_intervals */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
