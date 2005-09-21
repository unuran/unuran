/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_init.c                                                   *
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

/*---------------------------------------------------------------------------*/

#include "tdr_gw_init.ch"
#include "tdr_ps_init.ch"

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  int i,k;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_TDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tdr_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdr_debug_init_start(par,gen);
#endif

  /* get starting points */
  if (_unur_tdr_starting_cpoints(par,gen)!=UNUR_SUCCESS) {
    _unur_par_free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* compute intervals for given starting points */
  if (_unur_tdr_starting_intervals(par,gen)!=UNUR_SUCCESS) {
    _unur_par_free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* update maximal number of intervals */
  if (GEN->n_ivs > GEN->max_ivs) {
    GEN->max_ivs = GEN->n_ivs;
  }
  
  /* set boundaries for U */
  GEN->Umin = 0.;
  GEN->Umax = 1.;


  if (par->variant & TDR_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TDR_DEBUG_DARS) {
      /* make initial guide table (only necessary for writing debug info) */
      _unur_tdr_make_guide_table(gen);
      /* write info into log file */
      _unur_tdr_debug_dars_start(par,gen);
    }
#endif

    for (i=0; i<3; i++) {
      /* we make several tries */

      /* run DARS */
      if (_unur_tdr_run_dars(par,gen)!=UNUR_SUCCESS) {
	_unur_par_free(par); _unur_tdr_free(gen);
	return NULL;
      }
    
      /* make initial guide table */
      _unur_tdr_make_guide_table(gen);

      /* check if DARS was completed */
      if (GEN->n_ivs < GEN->max_ivs) {
	/* ran ARS instead */
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_tdr_debug_dars_finished(gen);
#endif
  }
  
  else { /* do not run DARS */
    /* make initial guide table */
    _unur_tdr_make_guide_table(gen);
  }

  /* free parameters */
  _unur_par_free(par);

  /* is there any hat at all ? */
  if (GEN->Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdr_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_tdr_debug_init_finished(gen);
#endif

  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

  /* o.k. */
  return gen;

} /* end of _unur_tdr_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_TDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_tdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* which transformation */
  if (PAR->c_T == 0.)
    par->variant = (par->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_LOG;
  else if (_unur_FP_same(PAR->c_T, -0.5))
    par->variant = (par->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_SQRT;
  else
    par->variant = (par->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_POW;

  /** TODO: remove this **/
  if ((par->variant & TDR_VARMASK_T) == TDR_VAR_T_POW) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"c != 0. and c != -0.5 not implemented!");
    _unur_distr_free(gen->distr); free(gen);
    return NULL;
  }

  /* routines for sampling and destroying generator */
  gen->destroy = _unur_tdr_free;
  gen->clone = _unur_tdr_clone;

  /* sampling routines */
  switch (par->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    SAMPLE = (par->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_gw_sample_check : _unur_tdr_gw_sample;
    break;
  case TDR_VARIANT_PS:    /* proportional squeeze */
    SAMPLE = (par->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
    break;
  case TDR_VARIANT_IA:    /* immediate acceptance */
    SAMPLE = (par->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ia_sample_check : _unur_tdr_ia_sample;
    break;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_distr_free(gen->distr); free(gen);
    return NULL;
  }
  
  /* set all pointers to NULL */
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */
  GEN->c_T = PAR->c_T;                /* parameter for transformation          */

  /* bounds for adding construction points  */
  GEN->max_ivs = max(2*PAR->n_starting_cpoints,PAR->max_ivs);  /* maximum number of intervals */
  GEN->max_ratio = PAR->max_ratio;    /* bound for ratio  Atotal / Asqueeze    */
  GEN->bound_for_adding = PAR->bound_for_adding;

  /* mode known and in given domain ?? */
  if ( !(par->distr->set & UNUR_DISTR_SET_MODE)
       || (DISTR.mode < DISTR.BD_LEFT)
       || (DISTR.mode > DISTR.BD_RIGHT))
    /* we cannot use the mode as construction point */
    par->variant = par->variant & (~TDR_VARFLAG_USEMODE);

  /* center known ?? */
  if (!(par->set & TDR_SET_CENTER))
    /* we cannot use the center as construction point */
    par->variant = par->variant & (~TDR_VARFLAG_USECENTER);
  else {
    /* center must be in domain */
    PAR->center = max(PAR->center,DISTR.BD_LEFT);
    PAR->center = min(PAR->center,DISTR.BD_RIGHT);
  }

  /* set default for DARS */
  if (!(par->set & TDR_SET_USE_DARS) && !PAR->starting_cpoints)
    /* no starting points given by user
       --> enable derandomized ARS      */
    par->variant |= TDR_VARFLAG_USEDARS;

  /* copy variant */
  gen->variant = par->variant;

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_tdr_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_tdr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_tdr_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_tdr_interval *iv,*next, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tdr_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
      clone_iv->prev = NULL;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
      clone_iv->prev = clone_prev;
    }
    /* next step */
    next = iv->next;
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;

  /* make new guide table */
  CLONE->guide = NULL;
  _unur_tdr_make_guide_table(clone);

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_tdr_clone() */

/*****************************************************************************/

void
_unur_tdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_TDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_tdr_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_tdr_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free table */
  if (GEN->guide)  free(GEN->guide);

  /* free other memory not stored in list */
  _unur_generic_free(gen);

} /* end of _unur_tdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_tdr_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting intervals.                  */
     /* if not provided as arguments compute these                           */
     /* by means of the "equiangular rule" from AROU.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int use_center, use_mode, is_mode, was_mode;
  int i, is_increasing;
  double extra_cpoint;
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  
  /* use mode as construction point ? */
  use_mode = (par->variant & TDR_VARFLAG_USEMODE) ? TRUE : FALSE;

  /* use center as construction point ? */
  use_center = (!use_mode && (par->variant & TDR_VARFLAG_USECENTER)) ? TRUE : FALSE;

  /* add extra construction point        */
  /* (use either mode or center or none) */
  extra_cpoint = use_mode ? DISTR.mode : (use_center ? PAR->center : 0. );

  /* reset counter of intervals */
  GEN->n_ivs = 0;

  /* prepare for computing construction points */
  if (!PAR->starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT  - PAR->center);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT - PAR->center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = DISTR.BD_LEFT;
  if (use_mode && DISTR.mode <= x) {
    /* this is the mode of the distribution */
    is_mode = TRUE;
    use_mode = FALSE;  /* do not use the mode again */
    is_increasing = FALSE;
  }
  else if (use_center && PAR->center <= x) {
    is_mode = FALSE;
    use_center = FALSE;     /* do not use the center again */
    is_increasing = TRUE;   /* the center may be left of (unknown) mode */
  }
  else {
    is_mode = FALSE;
    is_increasing = TRUE;
  }
    
  fx = fx_last = _unur_FP_is_minus_infinity(x) ? 0. : PDF(x);
  iv = GEN->iv = _unur_tdr_interval_new( gen, x, fx, is_mode );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */

  /* terminate beginning of list */
  iv->prev = NULL;

  /* now all the other points */
  for( i=0; i<=PAR->n_starting_cpoints; i++ ) {
    was_mode = is_mode;

    /* construction point */
    if (i < PAR->n_starting_cpoints) {
      if (PAR->starting_cpoints) {   
	/* construction points provided by user */
	x = PAR->starting_cpoints[i];
	/* check starting point */
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle ) + PAR->center;
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store 
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /* insert mode or center ? */
    if ((use_mode || use_center) && x >= extra_cpoint) {
      is_mode = use_mode;              /* the next construction point is the mode */
      use_center = use_mode = FALSE;   /* we use the mode only once (of course) */
      if (x>extra_cpoint) {
	x = extra_cpoint;     /* use the mode now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!PAR->starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x == extra_cpoint --> nothing to do */
    }
    else
      is_mode = FALSE;

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of PDF at starting point */
    fx = _unur_FP_is_infinity(x) ? 0. : PDF(x);

    /* check value of PDF at starting point */
    if (!is_increasing && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }
    if (is_mode && fx < fx_last * (1.-DBL_EPSILON)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"mode -> ignore");
      continue;
    }
    if (was_mode && fx > fx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"mode");
      return UNUR_ERR_GEN_DATA;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such point */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<PAR->n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the PDF is constant 0 on all construction points.
	     then we need both boundary points. */
	  iv->x = x;  /* we only have to change x, everything else remains unchanged */
	  x_last = x;
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with PDF(x) > 0 */
	break;
    }
    
    /* need a new interval */
    iv->next = _unur_tdr_interval_new( gen, x, fx, is_mode );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */
    
    /* link into list and skip pointer to current interval */
    iv->next->prev = iv;
    iv = iv->next;

    /* PDF still increasing ? */
    if (is_increasing && fx < fx_last)
      is_increasing = 0;

    /* store last computed values */
    x_last = x;
    fx_last = fx;

  }

  /* we have left the loop with the right boundary of the support of PDF
     make shure that we will never use iv for sampling. */
  iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
  iv->Acum = INFINITY;
  iv->ip = iv->x;
  iv->fip = iv->fx;
  iv->next = NULL;         /* terminate list */
  --(GEN->n_ivs);           /* we do not count this interval */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_starting_cpoints() */

/*****************************************************************************/

int
_unur_tdr_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
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
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    return _unur_tdr_gw_starting_intervals(par,gen);
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    return _unur_tdr_ps_starting_intervals(par,gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }
} /* end of _unur_tdr_interval_parameter() */

/*****************************************************************************/

int
_unur_tdr_run_dars( struct unur_par *par, struct unur_gen *gen )
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
  struct unur_tdr_interval *iv;
  double Atot, Asqueezetot;    /* total area below hat and squeeze, resp. */

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  
  /* there is no need to run DARS when the DARS factor is INFINITY */
  if (_unur_FP_is_infinity(PAR->darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_tdr_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }

  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;

  /* now run DARS for different variants */
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    return _unur_tdr_gw_dars(par,gen);
  case TDR_VARIANT_PS:    /* proportional squeeze */
  case TDR_VARIANT_IA:    /* immediate acceptance */
    return _unur_tdr_ps_dars(par,gen);
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

} /* end of _unur_tdr_run_dars() */

/*****************************************************************************/

struct unur_tdr_interval *
_unur_tdr_interval_new( struct unur_gen *gen, double x, double fx, int is_mode )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x       ... left point of new interval                             */
     /*   fx      ... value of PDF at x                                      */
     /*   is_mode ... if TRUE, x is a mode of the PDF                        */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double dfx;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return NULL;
  }
  if (_unur_FP_is_infinity(fx)) {
    /* over flow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_xmalloc( sizeof(struct unur_tdr_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN->n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_TDR_IV);

  /* avoid uninitialized variables */
  iv->Acum = iv->Ahat = iv->Ahatr = iv->Asqueeze = 0.;
  iv->ip = iv->fip = iv->sq = 0.;

  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->fx = fx;            /* value of PDF at x */
  
  if (fx<=0.) {           /* --> -INFINITY */
    iv->Tfx = -INFINITY;  /* transformed density */
    iv->dTfx = INFINITY;  /* derivative of transformed density */
    return iv;
  }
  
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    iv->Tfx = log(fx);
    dfx = dPDF(x);
    /* we can set dPDF(x) = 0. for the mode */
    if (is_mode || dfx==0.)
      iv->dTfx = 0.;
    else
      iv->dTfx = (1./fx * dfx);   /* possible overflow ? */
    break;
  case TDR_VAR_T_SQRT:
    iv->Tfx = -1./sqrt(fx);
    dfx = dPDF(x);
    /* we can set dPDF(x) = 0. for the mode */
    if (is_mode || dfx==0.)
      iv->dTfx = 0.;
    else
      iv->dTfx = (dfx<0.) ? -exp( -M_LN2 - 1.5*log(fx) + log(-dfx))
	: exp( -M_LN2 - 1.5*log(fx) + log(dfx));
    break;
  case TDR_VAR_T_POW:
    /** TODO **/
    /*      iv->Tfx = -pow(fx,GEN->c_T); */
    /*      iv->dTfx = 0.; */
    break;
  }

  /* the program requires dTfx > -INFINITY */
  if (iv->dTfx <= -INFINITY)
    iv->dTfx = INFINITY;

  return iv;

} /* end of _unur_tdr_interval_new() */

/*****************************************************************************/

int
_unur_tdr_tangent_intersection_point( struct unur_gen *gen, struct unur_tdr_interval *iv, double *ipt )
     /*----------------------------------------------------------------------*/
     /* compute cutting point of interval into left and right part.          */
     /* (1) use intersection point of tangents of transformed hat.           */
     /* (2) use mean point if (1) is unstable due to roundoff errors.        */
     /* (3) use boundary point which is closer to the mode. this is          */
     /*     important when the transformed tagents are extremely steep.      */
     /*     (This might cause a serious roundoff error while sampling.)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   ipt ... pointer to intersection point                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* 
     case: there is no tangent at one of the boundary points of the interval
           (then the slope is INFINITY)
     or
     case: the tangents are too steep  (--> case (3))
  */
  if ( iv->dTfx > 1.e+140 ) {
    *ipt = iv->x;        /* intersection point = left boundary of interval */
    return UNUR_SUCCESS; 
  }
  if ( iv->next->dTfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dTfx)) {
    *ipt = iv->next->x;   /* intersection point = right boundary of interval */
    return UNUR_SUCCESS; 
  }
  /** TODO: 1.e+140 (= sqrt(DBL_MAX) / 1.e15) is arbitrary  **/

  /* test for T-concavity */
  if ( _unur_FP_less( iv->dTfx, iv->next->dTfx ) ) {

    /* it might happen because of round-off errors 
       that iv->next->dTfx is almost zero although it should be large.
       thus we ignore this case. */
    if ( fabs(iv->dTfx) < DBL_EPSILON * fabs(iv->next->dTfx) ) {
      *ipt = iv->x;        /* intersection point = left boundary of interval */
      iv->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else if ( fabs(iv->next->dTfx) < DBL_EPSILON * fabs(iv->dTfx) ) {
      *ipt = iv->next->x;   /* intersection point = right boundary of interval */
      iv->next->dTfx = INFINITY;
      return UNUR_SUCCESS; 
    }
    else {
/*        fprintf(stdout,"\ndTfx0 = %g < %g = dTfx1 (x0 = %g, x1 = %g)\n", */
/*  	      iv->dTfx,iv->next->dTfx,iv->x,iv->next->x); */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /** TODO: the following test is too sensitve to roundoff errors **/
  /*    if (iv->next->Tfx > iv->x + iv->dTfx*(iv->next->x - iv->x)) { */
  /*      _unur_warning(gen->genid,UNUR_ERR_INIT,"tangent below PDF  not T-concave!"); */
  /*      return UNUR_ERR_INIT; */
  /*    } */
  
  /* case (2): computing intersection of tangents is unstable */
  if (_unur_FP_approx(iv->dTfx, iv->next->dTfx)) {
    /* use mean point */
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }

  /* case (1): compute intersection point of tangents (regular case) */
  *ipt = ( (iv->next->Tfx - iv->Tfx - iv->next->dTfx * iv->next->x + iv->dTfx * iv->x) / 
	   (iv->dTfx - iv->next->dTfx) );

  /* check position of intersection point */
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    /* intersection point of tangents not in interval.
       This is mostly the case for numerical reasons.
       Thus we is the center of the interval instead.
       if the PDF not T-concave, it will catched at a later
       point when we compare slope of tangents and squeeze. */
    *ipt = 0.5 * (iv->x + iv->next->x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_tangent_intersection_point() */

/*****************************************************************************/

double
_unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute area below piece of hat or squeeze in                             */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   area                                                                    */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /* log(x)                                                                    */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(t-x0)) dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /* -1/sqrt(x)                                                                */
     /*   area = | \int_{x0}^x 1/(Tf(x0) + slope*(t-x0))^2 dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = infinity                                      if T(f(x)) >= 0    */
     /*        = | (x-x0) / (Tf(x0)*(Tf(x0)+slope*(x-x0))) |   otherwise          */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double area = 0.;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* if the construction point is at infinity, we cannot compute an area.
     (in this case we should have x == iv->x == INFINITY). */
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return 0.;

  /* unbounded? */
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:
    /* T(x) = log(x) */
    if (slope != 0.) {                         
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = iv->fx / slope;
      else {
	double t = slope * (x - iv->x);
	if (fabs(t) > 1.e-6)
	  area = iv->fx * (x - iv->x) * ( exp(t) - 1. ) / t;
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  area = iv->fx * (x - iv->x) * (1. + t/2. + t*t/6.);
	else
	  area = iv->fx * (x - iv->x) * (1. + t/2.);
      }
    }
    else { /* hat/squeeze almost constant */
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_SQRT:
    /* T(x) = -1./sqrt(x) */
    if (slope != 0.) {
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	area = 1. / ( iv->Tfx * slope );
      else {
	/* compute value of transformed hat at integration boundary */
	double hx = iv->Tfx + slope * (x - iv->x);
	/* the transformed hat must always be below the x-axis.
	   otherwise the area below the hat in unbounded. */
	if (hx>=0.)
	  return INFINITY; 
	else
	  area = (x - iv->x) / ( iv->Tfx * hx );
      }
    }
    else { /* hat/squeeze almost constant */
      if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
	return INFINITY;
      area = iv->fx * (x - iv->x);
    }
    break;

  case TDR_VAR_T_POW:
    /* T(x) = -x^c */
    /** TODO **/
    break;
  }

  return ( (area<0.) ? -area : area );

} /* end of _unur_tdr_interval_area() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_interval_xxarea( struct unur_gen *gen, struct unur_tdr_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute the interal of x times hat or squeeze ("expected value") in       */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   "expected value"                                                        */
     /*   (to get the real expected value, it must be divided by the area below   */
     /*   the function (hat or squeeze).)                                         */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /* log(x)                                                                    */
     /*   ev = \int_{x0}^x t * \exp(Tf(x0) + slope*(t-x0)) dt                     */
     /*      = 0.5 * f(x0) * (x^2 - x0^2)                            if slope = 0 */
     /*      = f(x0)/slope^2 * (\exp(slope*(x-x0))*(slope*x-1) - (slope*x0-1))    */
     /*                                                             if slope != 0 */
     /*                                                                           */
     /* -1/sqrt(x)                                                                */
     /*   ev = \int_{x0}^x t / (Tf(x0) + slope*(t-x0))^2 dt                       */
     /*      = 0.5 * f(x0) * (x^2 - x0^2)                            if slope = 0 */
     /*      = infinity                         if T(f(x)) >= 0 or |x| = infinity */
     /*      = x0 / (slope*Tf(x0)) - x / (slope*u) + log(u/Tf(x0)) / slope^2      */
     /*        where u = Tf(x0) + slope*(x-x0)                          otherwise */
     /*                                                                           */
     /* To get the right sign we have to return sign(x-x0) * ev.                  */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double ev = 0.;
  double hx,u;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* if the construction point is at infinity, we cannot compute the integral
     (in this case we should have x == iv->x == INFINITY). */
  if (_unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x))
    return 0.;

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return 0.;

  /* unbounded? */
  if ( _unur_FP_is_infinity(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;


  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:    /* T(x) = log(x) */
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x)) {
      ev = iv->fx / (slope*slope) * (1-slope*iv->x);
    }
    else {
      u = (x-iv->x) * slope;

      if (fabs(u) > 1.e-6) {
	ev = iv->fx / (slope*slope) * (exp(u)*(slope*x-1.) - slope*iv->x + 1.);
      }
      else {
	/* use Taylor series */
	/* constant term */
	ev = 0.5 * (x+iv->x);
	if (fabs(u) > 0) {
	  /* 1st order expansion */
	  ev += 1./6. * (2.*x+iv->x) * u;
	  /* 2nd order expansion */
	  ev += 1./24. * (3.*x+iv->x) * u * u;
	}
	ev *= iv->fx * (x-iv->x);
      }
    }
    break;

  case TDR_VAR_T_SQRT:    /* T(x) = -1./sqrt(x) */
    if (_unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x))
      /* the integral becomes INFINITY */
      return INFINITY;

    /* compute value of transformed hat at integration boundary */
    hx = iv->Tfx + slope * (x - iv->x);

    if (hx >= 0.)
      /* the transformed hat must always be below the x-axis.
	 otherwise the area below the hat in unbounded. */
      return INFINITY; 

    u = (x-iv->x) * slope / iv->Tfx;

    if (fabs(u) > 1.e-6) {
      ev = ( iv->x / (slope * iv->Tfx) - x / (slope * hx)
	     + log( hx / iv->Tfx ) / (slope*slope) );
    }
    else {
      /* use Taylor series */
      /* constant term */
      ev = 0.5 * (x+iv->x);
      if (fabs(u) > 0) {
	/* 1st order expansion */
	ev -= 1./3. * (2.*x+iv->x) * u;
	/* 2nd order expansion */
	ev += 1./4. * (3.*x+iv->x) * u * u;
      }
      ev *= iv->fx * (x-iv->x);
    }
    break;

  case TDR_VAR_T_POW:    /* T(x) = -x^c */
    /** TODO **/
    break;
  }

  return ((x>iv->x) ? ev : -ev);

} /* end of _unur_tdr_interval_xxarea() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_eval_intervalhat( struct unur_gen *gen, struct unur_tdr_interval *iv, double x )
     /*----------------------------------------------------------------------*/
     /* evaluate hat at x in interval.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval that stores constr. point of tangent   */
     /*   x   ... point at which hat(x) has to be computed                   */
     /*                                                                      */
     /* return:                                                              */
     /*   hat(x) or                                                          */
     /*   0. if x is not finite or                                           */
     /*   INFINITY the hat cannot be computed                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /* comment:                                                             */
     /*   x0    ... construction point of tangent (= iv->x)                  */
     /*                                                                      */
     /* log(x)                                                               */
     /*   hat(x) = f(x0) * exp( Tf'(x0)(x - x_0) )                           */
     /*                                                                      */
     /* 1/sqrt(x)                                                            */
     /*   hat(x) = 1/(Tf(x0) + Tf'(x0)(x - x_0))^2                           */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_TDR_IV,INFINITY); 

  /* we cannot compute the hat at x if any of the parameters are not finite  */
  if ( _unur_FP_is_minus_infinity(iv->Tfx) || _unur_FP_is_infinity(iv->dTfx) )
    return INFINITY;

  /* at +/- infinity the hat should be 0 (or infinity) */
  if ( _unur_FP_is_infinity(x) || _unur_FP_is_minus_infinity(x) ||
       _unur_FP_is_infinity(iv->x) || _unur_FP_is_minus_infinity(iv->x) )
    return 0.;

  /* now evaluate hat at x */
  switch( gen->variant & TDR_VARMASK_T ) {

  case TDR_VAR_T_LOG:
    /* T(x) = log(x) */
    return (iv->fx * exp( iv->dTfx * (x - iv->x) ));

  case TDR_VAR_T_SQRT:
    /* T(x) = -1./sqrt(x) */
    {
      /* compute value of transformed hat at x */
      double hx = iv->Tfx + iv->dTfx * (x - iv->x);
      /* hx must be less than 0 ! */
      return ((hx<0.) ? 1./(hx*hx) : INFINITY);
    }

  case TDR_VAR_T_POW:
    /* T(x) = -1./x^c */
    /** TODO **/
    return INFINITY;

  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }

} /* end of _unur_tdr_eval_intervalhat() */

/*****************************************************************************/

int
_unur_tdr_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    int max_guide_size = (GEN->guide_factor > 0.) ? (GEN->max_ivs * GEN->guide_factor) : 1;
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_tdr_interval*) );
  }

  /* first we need cumulated areas in intervals */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN->guide_size = (int)(GEN->n_ivs * GEN->guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, iv=GEN->iv; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   /* this is the last virtual intervall --> do not use */
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = iv;

  return UNUR_SUCCESS;
} /* end of _unur_tdr_make_guide_table() */

/*****************************************************************************/

