/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl_init.c                                                  *
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
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_tabl_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_TABL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tabl_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tabl_debug_init_start(par,gen);
#endif

  /* we need slopes */
  if (PAR->n_slopes <= 0 ) {
    /* no slopes are given  --> compute from mode */
    if (par->distr->set & UNUR_DISTR_SET_MODE) {
      double slopes[4];
      PAR->n_slopes = _unur_tabl_get_slopes_from_mode(gen,slopes);
      PAR->slopes = slopes; 
    }
  }

  /* now we must have slopes */
  if (PAR->n_slopes <= 0 ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"number of slopes <= 0, cannot compute");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }

  /* get starting intervals from slopes */
  if (_unur_tabl_get_starting_intervals(par,gen)!=UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make hat function");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }

  /* split starting intervals / slopes */
  if (_unur_tabl_compute_intervals(par,gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot split intervals");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }

  /* make initial guide table */
  if (_unur_tabl_make_guide_table(gen) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
    _unur_par_free(par); _unur_tabl_free(gen); return NULL;
  }

  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_tabl_debug_init_finished(par,gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_tabl_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tabl_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_tabl_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TABL_GEN);

  /* check for required data: area */
  if (!(gen->distr->set & UNUR_DISTR_SET_PDFAREA))
    if (unur_distr_cont_upd_pdfarea(gen->distr)!=UNUR_SUCCESS)
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"area below PDF, use default instead");

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & TABL_VARFLAG_VERIFY) ? _unur_tabl_sample_check : _unur_tabl_sample;
  gen->destroy = _unur_tabl_free;
  gen->clone = _unur_tabl_clone;

  /* set all pointers to NULL */
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;

  /* the boundaries for our computation limits are intersection of the       */
  /* domain of the distribution and the given computation boundaries.        */
  if (par->distr->set & UNUR_DISTR_SET_DOMAIN) {
    PAR->bleft  = max(PAR->bleft, DISTR.BD_LEFT);
    PAR->bright = min(PAR->bright,DISTR.BD_RIGHT);
  }
  GEN->bleft       = PAR->bleft;         /* left boundary of domain            */
  GEN->bright      = PAR->bright;        /* right boundary of domain           */

  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN->max_ivs   = PAR->max_ivs;         /* maximum number of intervals        */
  GEN->max_ratio = PAR->max_ratio;       /* bound for ratio  Atotal / Asqueeze */
  GEN->darsfactor = PAR->darsfactor;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_tabl_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tabl_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE  ((struct unur_tabl_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_tabl_interval *iv,*next, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_tabl_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tabl_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
    }
    /* next step */
    next = iv->next;
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;

  /* make new guide table */
  CLONE->guide = NULL;
  if (_unur_tabl_make_guide_table(clone) != UNUR_SUCCESS) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create guide table");
    /* There is no chance to run out of this error! */
    /* however, this should never happen.           */
  }

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tabl_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_tabl_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_TABL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tabl_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_tabl_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free table */
  if (GEN->guide)  free(GEN->guide);

  /* free other memory */
  _unur_generic_free(gen);

} /* end of _unur_tabl_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_tabl_get_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals, slopes are given by user.                */
     /* estimate domain when not given.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   PDF(a) >= PDF(b)                                                   */
     /*----------------------------------------------------------------------*/
{
  /** TODO: check for slopes out of support !! **/

  struct unur_tabl_interval *iv;
  int i;

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* init counter of intervals */
  GEN->n_ivs = 0;
  iv = NULL;

  /* boundary of computational interval are reset by boundaries of slopes */
  GEN->bleft = INFINITY;
  GEN->bright = -INFINITY;

  /* compute initial intervals */
  for ( i=0; i < 2*PAR->n_slopes; i+=2 ) {
    /* get a new interval and link into list */
    if (i==0)  /* the first interval */
      iv = GEN->iv = _unur_xmalloc(sizeof(struct unur_tabl_interval));   
    else       /* all the other intervals */
      iv = iv->next = _unur_xmalloc(sizeof(struct unur_tabl_interval));
    ++(GEN->n_ivs);
    COOKIE_SET(iv,CK_TABL_IV);

    /* max and min of PDF in interval */
    iv->xmax = PAR->slopes[i];      
    iv->fmax = PDF(iv->xmax);
    iv->xmin = PAR->slopes[i+1];    
    iv->fmin = PDF(iv->xmin);

    /* check for overflow */
    if (! (_unur_isfinite(iv->fmax) && _unur_isfinite(iv->fmin)) ) {
      /* overflow */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      iv->next = NULL;  /* terminate list (freeing list) */
      return UNUR_ERR_GEN_DATA;
    }

    /* check slopes */
    if (_unur_FP_less(iv->fmax,iv->fmin)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"slopes non-decreasing");
      iv->next = NULL;  /* terminate list (freeing list) */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area of slope */
    iv->Ahat = fabs(iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = fabs(iv->xmax - iv->xmin) * iv->fmin;
    /** TODO: possible overflow/underflow ?? **/
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* estimate domain */
    if (iv->xmax > iv->xmin) {
      /* increasing slope */
      GEN->bleft = min(GEN->bleft,iv->xmin);
      GEN->bright = max(GEN->bright,iv->xmax);
    }
    else {
      /* decreasing slope */
      GEN->bleft = min(GEN->bleft,iv->xmax);
      GEN->bright = max(GEN->bright,iv->xmin);
    }
  }

  /* terminate list */
  iv->next = NULL;

  /* reset domain of distribution */
  DISTR.trunc[0] = DISTR.BD_LEFT = GEN->bleft;
  DISTR.trunc[1] = DISTR.BD_RIGHT = GEN->bright;

  /* reset area below distribution */
  gen->distr->set &= ~UNUR_DISTR_SET_PDFAREA;
  unur_distr_cont_upd_pdfarea( gen->distr );

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_get_starting_intervals() */

/*---------------------------------------------------------------------------*/

int
_unur_tabl_get_slopes_from_mode( struct unur_gen *gen, double *slopes )
     /*----------------------------------------------------------------------*/
     /* compute slopes from mode and computational domain                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   slopes ... pointer to array of length 4 to store slopes            */
     /*                                                                      */
     /* return:                                                              */
     /*   number of slopes                                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  if (DISTR.mode <= GEN->bleft) {
    /* only one ascending interval <a,b> = [a,b] */
    slopes[0] = GEN->bleft;
    slopes[1] = GEN->bright;
    return 1;
  }

  if (DISTR.mode >= GEN->bright) {
    /* only one descending interval <a,b> = [b,a] */
    slopes[0] = GEN->bright;
    slopes[1] = GEN->bleft;
    return 1;
  }

  else {
    /* two slopes */
    /* one ascending interval */
    slopes[0] = DISTR.mode;
    slopes[1] = GEN->bleft;
    /* one descending */
    slopes[2] = DISTR.mode;
    slopes[3] = GEN->bright;
    return 2;
  }

} /* end of _unur_tabl_get_slopes_from_mode() */

/*---------------------------------------------------------------------------*/

int
_unur_tabl_compute_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute all intervals (by splitting starting intervals/slopes)       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  int i,k;

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* split interval following [1], split A (equal area rule) */
  if (par->variant & TABL_VARFLAG_STP_A) {
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
      iv = _unur_tabl_run_equalarearule( par, gen, iv );
      if (iv == NULL) return UNUR_ERR_GEN_DATA;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* print intervals after starting intervals have been created */ 
  if (gen->debug & TABL_DEBUG_A_IV) _unur_tabl_debug_intervals(gen,FALSE);
#endif

  if (par->variant & TABL_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */
 
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TABL_DEBUG_DARS) _unur_tabl_debug_dars_start(par,gen);
#endif

    for (i=0; i<TABL_N_RETRY_DARS; i++) {
      /* we make several tries */
      
      /* run DARS */
      if (_unur_tabl_run_dars(par,gen)!=UNUR_SUCCESS)
	return UNUR_ERR_GEN_DATA;

      /* check if DARS was completed */
      if (GEN->n_ivs >= GEN->max_ivs)
	break;

      /* else ran ARS instead */
      for (k=0; k<TABL_N_RUN_ARS; k++)
	_unur_sample_cont(gen);
    }

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TABL_DEBUG_DARS) _unur_tabl_debug_dars(par,gen);
#endif
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_tabl_compute_intervals() */

/*---------------------------------------------------------------------------*/

struct unur_tabl_interval *
_unur_tabl_run_equalarearule( struct unur_par *par, struct unur_gen *gen, 
			      struct unur_tabl_interval *iv_slope )
     /*----------------------------------------------------------------------*/
     /* split starting intervals according to [1]                            */
     /* SPLIT A (equal areas rule)                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter list                            */
     /*   gen       ... pointer to generator object                          */
     /*   iv_slope  ... pointer to interval of slope                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to last interval in list of splitted slope                 */
     /*   NULL on error                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv, *iv_last;
  double bar_area, x, fx;
  double slope;

  /* check arguments */
  CHECK_NULL(par,NULL);       COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  CHECK_NULL(gen,NULL);       COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  CHECK_NULL(iv_slope,NULL);  COOKIE_CHECK(iv_slope,CK_TABL_IV,NULL);

  if (GEN->n_ivs >= GEN->max_ivs) 
    /* maximal number of intervals reached */
    return NULL;

  iv = iv_slope;        /* pointer to actual interval */
  iv_last = iv_slope;   /* pointer to last interval in list */
  /* (maximal) area of bar (= hat in one interval) */
  bar_area = DISTR.area * PAR->area_fract;

  while (_unur_FP_greater(iv->Ahat, bar_area)) {

    /* increasing or decreasing slope */
    slope = (iv->xmax > iv->xmin) ? 1. : -1.;

    /* compute splitting point:
       slope == +1 --> move from right to left
       slope == -1 --> move from left to right */
    x = iv->xmax - slope * bar_area / iv->fmax;

    /* compute PDF */
    fx = PDF(x);

    /* check for overflow */
    if (_unur_FP_is_infinity(fx)) {
      /* overflow */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
      return NULL;
    }

    /* now split interval at x */
    switch (_unur_tabl_split_interval( gen, iv, x, fx, TABL_VARFLAG_SPLIT_POINT )) {
    case UNUR_SUCCESS:  /* splitting succesful */
      if (slope > 0.) {
	if (iv_last == iv_slope)
	  iv_last = iv->next;
      }
      else { /* slope < 0 */
	iv = iv->next; break;
      }
      break;
    case UNUR_ERR_SILENT: /* interval chopped */
      break; /* nothing to do */
    default: /* error (slope not monotonically increasing) */
      return NULL;
    }

    /* check number of intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* maximal number of intervals reached */
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"split A stopped, maximal number of intervals reached.");
      break;
    }
  }

  /* pointer to last interval */
  return ((iv->xmax > iv->xmin) ? iv_last : iv);

} /* end of _unur_tabl_run_equalarearule() */

/*---------------------------------------------------------------------------*/

int
_unur_tabl_run_dars( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Atot, Asqueezetot;    /* total area below hat and squeeze, resp. */
  double Alimit;               /* threshhold value for splitting interval */
  int n_splitted = 1;          /* count splitted intervals */
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TABL_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* there is no need to run DARS when the DARS factor is INFINITY */
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_arou_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;

  /* now split intervals */
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_ivs < GEN->max_ivs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN->n_ivs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_ivs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 

    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);

      /* do not exceed the maximum number of intervals */
      if (GEN->n_ivs >= GEN->max_ivs)
	break;

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  /* goto next interval */

      switch (_unur_tabl_split_interval( gen, iv, 0., 0., TABL_VARFLAG_SPLIT_ARC )) {
      case UNUR_SUCCESS:  /* splitting succesful */
      case UNUR_ERR_SILENT: /* interval chopped */
	++n_splitted;
	break; /* nothing to do */
      default: /* error (slope not monotonically decreasing) */
	return UNUR_ERR_GEN_DATA;
      }
    }

    if (n_splitted == 0) {
      /* we are not successful in splitting any inteval.
	 abort to avoid endless loop */
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }

  /* ratio between squeeze and hat o.k. ? */
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_ivs >= GEN->max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN->max_ivs = GEN->n_ivs;
  }
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_run_dars() */

/*****************************************************************************/

int
_unur_tabl_split_interval( struct unur_gen *gen,
			   struct unur_tabl_interval *iv_old, 
			   double x, double fx, 
			   unsigned split_mode )
     /*----------------------------------------------------------------------*/
     /* split interval (replace old one by two new ones in same place)       */
     /* new interval is inserted immedately after old one.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   iv_old     ... pointer to interval that has to be split            */
     /*   x          ... splitting point                                     */
     /*   fx         ... value of PDF at splitting point                     */
     /*   split_mode ... how to split interval                               */
     /*                  TABL_VARFLAG_SPLIT_POINT: split at given point x    */
     /*                  TABL_VARFLAG_SPLIT_MEAN:  at mean point of interval */
     /*                  TABL_VARFLAG_SPLIT_ARC:   at arc mean point         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... splitting successful                           */
     /*   UNUR_ERR_SILENT ... interval chopped off domain                    */
     /*                       (not part of support of PDF)                   */
     /*   others          ... error: PDF not monotone in interval            */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv_new;
  double A_hat_old, A_squ_old;
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);     COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_old,UNUR_ERR_NULL);  COOKIE_CHECK(iv_old,CK_TABL_IV,UNUR_ERR_COOKIE);

  /* There are three possibilities for the splitting point:
     (1) use x and avoid computation of PDF(x). 
     (2) use middle of interval. converges faster in many cases.
     (3) use "arc_mean" of interval. 
         converges faster when domain is almost unbounded. */
  switch( split_mode ) {
  case TABL_VARFLAG_SPLIT_POINT:    /* (1) */
    /* nothing to do (default) */
    break;
  case TABL_VARFLAG_SPLIT_MEAN:     /* (2) */
    x = 0.5 * (iv_old->xmin + iv_old->xmax); 
    fx = PDF(x);
    break;
  case TABL_VARFLAG_SPLIT_ARC:      /* (3) */
    x = _unur_arcmean(iv_old->xmin, iv_old->xmax); 
    fx = PDF(x);
    break;
  default: 
    /* this should not happen:
       Invalid variant, use default n*/
    _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    break;
  }

  /* check for overflow */
  if (!_unur_isfinite(fx)) {
    /* overflow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return UNUR_ERR_GEN_DATA;
  }

  /* store areas of old interval */
  A_hat_old = iv_old->Ahat;
  A_squ_old = iv_old->Asqueeze;

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {
    /* check montonicity */
    if (iv_old->fmin > 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not monotone in slope");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* chop off part out of support */
    iv_old->xmin = x;

    /* compute new area in interval */
    /** TODO: possible overflow/underflow ?? **/
    iv_old->Ahat = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmax;
    /* iv_old->Asqueeze remains 0 */

    /* update total area */
    GEN->Atotal += iv_old->Ahat - A_hat_old;
    /* GEN->Asqueeze remains unchanged */

    /* check result */
    if (!_unur_isfinite(GEN->Atotal)) {
      /* we have decreased the hat, thus ... */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_INF;
    }

    /* interval chopped but not split */
    return UNUR_ERR_SILENT;
  }

  /* we need a new interval */
  iv_new = _unur_xmalloc(sizeof(struct unur_tabl_interval));
  ++(GEN->n_ivs);
  COOKIE_SET(iv_new,CK_TABL_IV);

  /* iv_new has the same slope as iv_old */

  /* we have to distinguish between two cases:
     PDF is increasing (slope = +1) or
     PDF is decreasing (slope = -1). */

  if (iv_old->xmax > iv_old->xmin) {
    /* increasing slope */
    /* (x) The iv_new inherits the maximum of iv_old.
           iv_old keeps the minimum.
       (x) The splitting point is the minimum of iv_new and
           the maximum of iv_old.
    */
    iv_new->xmax  = iv_old->xmax;  
    iv_new->fmax = iv_old->fmax;
    iv_old->xmax  = iv_new->xmin = x; 
    iv_old->fmax = iv_new->fmin = fx; 
  }
  else {
    /* decreasing slope */
    /* (x) The iv_new inherits the minimum of iv_old.
           iv_old keeps the maximum.
       (x) The splitting point is the maximum of iv_new and
           the minimum of iv_old.
    */
    iv_new->xmin  = iv_old->xmin;  
    iv_new->fmin = iv_old->fmin;
    iv_old->xmin  = iv_new->xmax = x; 
    iv_old->fmin = iv_new->fmax = fx; 
  }

  /* compute the areas in both intervals */
  iv_new->Ahat     = fabs(iv_new->xmax - iv_new->xmin) * iv_new->fmax;
  iv_new->Asqueeze = fabs(iv_new->xmax - iv_new->xmin) * iv_new->fmin;
  iv_old->Ahat     = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmax;
  iv_old->Asqueeze = fabs(iv_old->xmax - iv_old->xmin) * iv_old->fmin;

  /* update total areas */
  GEN->Atotal += iv_old->Ahat + iv_new->Ahat - A_hat_old;
  GEN->Asqueeze += iv_old->Asqueeze + iv_new->Asqueeze - A_squ_old;

  /* insert iv_new into linked list of intervals.
     iv_old is stored on the left hand side of iv_new. */
  iv_new->next = iv_old->next;
  iv_old->next = iv_new;

  /* check result */
  if (! (_unur_isfinite(GEN->Atotal) && _unur_isfinite(GEN->Asqueeze)) ) {
    _unur_error(gen->genid,UNUR_ERR_INF,"hat unbounded");
    return UNUR_ERR_INF;
  }

  /* splitting successful */
  return UNUR_SUCCESS;

} /* end of _unur_tabl_split_interval() */

/*****************************************************************************/

int
_unur_tabl_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*   1 (--> successful)                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TABL_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    int max_guide_size = (GEN->guide_factor > 0.) ? (GEN->max_ivs * GEN->guide_factor) : 1;
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tabl_interval*) );
  }

  /* first we need the cumulated areas of rectangles */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
    
  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN->guide_size = GEN->n_ivs;

  /* make table (use variant 2; see dis.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TABL_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      if( iv->next != NULL )    /* skip to next segment if it exists */
        iv = iv->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;

  /* check table */
  if (! (_unur_isfinite(GEN->Atotal) && _unur_isfinite(GEN->Asqueeze)) ) {
    /* in this case the guide table is corrupted and completely worthless. */
    /* this error is unlikely to happen. for this rare case the check is   */
    /* done at the end to have the table at least filled with numbers      */
    /* and thus avoid infinite loops in the sampling routines.             */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"sum of areas not finite");
    return UNUR_ERR_GEN_DATA;
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_tabl_make_guide_table() */

/*---------------------------------------------------------------------------*/

