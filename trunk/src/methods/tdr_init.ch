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
  if (!gen) { free(par); return NULL; }

  /* get starting points */
  if (_unur_tdr_starting_cpoints(par,gen)!=UNUR_SUCCESS) {
    free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* compute intervals for given starting points */
  if (_unur_tdr_starting_intervals(par,gen)!=UNUR_SUCCESS) {
    free(par); _unur_tdr_free(gen);
    return NULL;
  }

  /* update maximal number of intervals */
  if (GEN.n_ivs > GEN.max_ivs) {
    GEN.max_ivs = GEN.n_ivs;
  }
  
  /* set boundaries for U */
  GEN.Umin = 0.;
  GEN.Umax = 1.;


  if (par->variant & TDR_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & TDR_DEBUG_DARS) {
      /* make initial guide table (only necessary for writing debug info) */
      _unur_tdr_make_guide_table(gen);
      /* write info into log file */
      _unur_tdr_debug_init(par,gen);
      _unur_tdr_debug_dars_start(par,gen);
    }
#endif

    for (i=0; i<3; i++) {
      /* we make several tries */

      /* run DARS */
      if (_unur_tdr_run_dars(par,gen)!=UNUR_SUCCESS) {
	free(par); _unur_tdr_free(gen);
	return NULL;
      }
    
      /* make initial guide table */
      _unur_tdr_make_guide_table(gen);

      /* check if DARS was completed */
      if (GEN.n_ivs < GEN.max_ivs) {
	/* ran ARS instead */
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) {
      if (gen->debug & TDR_DEBUG_DARS)
	_unur_tdr_debug_dars(par,gen);
      else 
	_unur_tdr_debug_init(par,gen);
    }
#endif
  }
  
  else { /* do not run DARS */
    /* make initial guide table */
    _unur_tdr_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_tdr_debug_init(par,gen);
#endif
  }

  /* free parameters */
  free(par);

  /* is there any hat at all ? */
  if (GEN.Atotal <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_tdr_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of _unur_tdr_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
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
  gen = _unur_generic_create( par );

  /* magic cookies */
  COOKIE_SET(gen,CK_TDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* which transformation */
  if (PAR.c_T == 0.)
    par->variant = (par->variant & (~TDR_VARMASK_T)) | TDR_VAR_T_LOG;
  else if (_unur_FP_same(PAR.c_T, -0.5))
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
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.iv          = NULL;
  GEN.n_ivs       = 0;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */
  GEN.c_T = PAR.c_T;                /* parameter for transformation          */

  /* bounds for adding construction points  */
  GEN.max_ivs = max(2*PAR.n_starting_cpoints,PAR.max_ivs);  /* maximum number of intervals */
  GEN.max_ratio = PAR.max_ratio;    /* bound for ratio  Atotal / Asqueeze    */
  GEN.bound_for_adding = PAR.bound_for_adding;

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
    PAR.center = max(PAR.center,DISTR.BD_LEFT);
    PAR.center = min(PAR.center,DISTR.BD_RIGHT);
  }

  /* set default for DARS */
  if (!(par->set & TDR_SET_USE_DARS) && !PAR.starting_cpoints)
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
#define CLONE clone->data.tdr

  struct unur_gen *clone;
  struct unur_tdr_interval *iv,*next, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN.iv; iv != NULL; iv = next) {
    /* copy segment */
    clone_iv = _unur_malloc( sizeof(struct unur_tdr_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_tdr_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE.iv = clone_iv;
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
  CLONE.guide = NULL;
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
    for (iv = GEN.iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free table */
  if (GEN.guide)  free(GEN.guide);

  /* free other memory not stored in list */
  _unur_generic_free(gen);

} /* end of _unur_tdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static int
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
  extra_cpoint = use_mode ? DISTR.mode : (use_center ? PAR.center : 0. );

  /* reset counter of intervals */
  GEN.n_ivs = 0;

  /* prepare for computing construction points */
  if (!PAR.starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT  - PAR.center);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT - PAR.center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR.n_starting_cpoints + 1);
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
  else if (use_center && PAR.center <= x) {
    is_mode = FALSE;
    use_center = FALSE;     /* do not use the center again */
    is_increasing = TRUE;   /* the center may be left of (unknown) mode */
  }
  else {
    is_mode = FALSE;
    is_increasing = TRUE;
  }
    
  fx = fx_last = _unur_FP_is_minus_infinity(x) ? 0. : PDF(x);
  iv = GEN.iv = _unur_tdr_interval_new( gen, x, fx, is_mode );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */

  /* terminate beginning of list */
  iv->prev = NULL;

  /* now all the other points */
  for( i=0; i<=PAR.n_starting_cpoints; i++ ) {
    was_mode = is_mode;

    /* construction point */
    if (i < PAR.n_starting_cpoints) {
      if (PAR.starting_cpoints) {   
	/* construction points provided by user */
	x = PAR.starting_cpoints[i];
	/* check starting point */
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle ) + PAR.center;
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
	if (!PAR.starting_cpoints)
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
    if (was_mode && fx > fx_last * (1+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"mode");
      return UNUR_ERR_GEN_DATA;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such point */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<PAR.n_starting_cpoints) {
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
  --(GEN.n_ivs);           /* we do not count this interval */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_starting_cpoints() */

/*****************************************************************************/

static int
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

/*---------------------------------------------------------------------------*/

static int
_unur_tdr_gw_starting_intervals( struct unur_par *par, struct unur_gen *gen )
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
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              /* construction point, value of PDF at x */
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN.iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN.iv,CK_TDR_IV,UNUR_ERR_COOKIE); 
  
  /* compute paramters for all intervals */
  for( iv=GEN.iv; iv->next != NULL; ) {

    /* compute parameters for interval */
    switch (_unur_tdr_gw_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:   /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN.n_ivs);
      
      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make shure that we will never use this segment */
	iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      else
	/* we need a pointer to the previous entry in the list */
	iv->next->prev = iv;
      continue;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area below hat infinite.
       insert new construction point. */
    x = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */

    /* value of PDF at x */
    fx = PDF(x);

    /* add a new interval, but check if we had to used too many intervals */
    if (GEN.n_ivs >= GEN.max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */


    /* if fx is 0, then we can cut off the tail of the distribution
       (since it must be T-concave)  */
    if (fx <= 0.) {
      if (iv->fx <= 0.) {
	/* cut off left tail */
	iv_new->next = iv->next;
	free(iv); 
	--(GEN.n_ivs);
	GEN.iv = iv_new;
	iv_new->prev = NULL;
	/* compute the parameters for the new left part */
	iv = iv_new;
      }
      else if (iv->next->fx <= 0.) {
	/* cut off right tail */
	free(iv->next);
	--(GEN.n_ivs);	
	iv->next = iv_new;
	iv_new->prev = iv;
	/* compute the paramters for the new part */
	/* (nothing to do here) */
      }
      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    else {
      /* insert new interval into linked list */
      iv_new->prev = iv;
      iv_new->next = iv->next;
      iv->next->prev = iv_new;
      iv->next = iv_new;
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_starting_intervals() */

/*---------------------------------------------------------------------------*/

static int
_unur_tdr_ps_starting_intervals( struct unur_par *par, struct unur_gen *gen )
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
  struct unur_tdr_interval *iv, *iv_new, *iv_tmp; 
  double x,fx;              /* construction point, value of PDF at x */
  double lb, flb;           /* left boundary point of domain of PDF  */
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN.iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN.iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* first interval in list */
  iv = GEN.iv;

  /* the left boundary of the domain:
     iv->x in the first interval is always the left boundary of 
     the domain of the PDF */
  lb = iv->x;
  flb = iv->fx;

  /* there is no use for a point iv->x that is not used as a 
     construction point in variants PS and IA.
     (In variant GW it is used to store the left boundary of the domain
     of the PDF)
     Thus we remove it from the list. At such points the slope of the
     tangents to the transformed density is set to INFINITY. */
  if (_unur_FP_is_infinity(iv->dTfx)) {
    GEN.iv = iv->next;
    GEN.iv->prev = NULL;
    free (iv);
    --(GEN.n_ivs);
    iv = GEN.iv;
  }

  /* set left boundary:
     it is stored in iv->ip in the first interval */
  iv->ip = lb;
  iv->fip = flb;

  /* compute paramters for all intervals */
  while (iv) {
    if (iv->next == NULL) {
      /* the last interval in the list*/

      /* analogously to variant GW we want to have a last virtual 
	 (stopping) interval, that simply stores the right boundary
	 point and guaranties that the loop in indexed search stops
	 on an existing interval.
	 However in the case where iv->x is used as construction point
	 we have to add such an interval. */
      if (!_unur_FP_is_infinity(iv->dTfx)) {
	/* get interval */
	iv->next = iv_new = _unur_tdr_interval_new( gen, iv->x, 0., FALSE );
	if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */
	/* link into list */
	iv_new->prev = iv;
	/* copy right boundary of domain */
	iv_new->ip = iv->x;
	iv_new->fip = iv->fx;
	/* make shure that we will never use this interval for sampling. */
	iv->next->Asqueeze = iv->next->Ahat = iv->next->Ahatr = 0.;
	iv->Acum = INFINITY;
	iv->next-> sq = 0.;
	/* we even have to to some additional work */
      }
      else
	/* nothing to do any more */
	break;
	
    }

    /* compute parameters for interval */
    switch (_unur_tdr_ps_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:     /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:    /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:    /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN.n_ivs);
      
      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make shure that we will never use this segment */
	iv->Asqueeze = iv->Ahat = iv->Ahatr = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      else
	/* we need a pointer to the previous entry in the list */
	iv->next->prev = iv;
      continue;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* check if we had to used too many intervals */
    if (GEN.n_ivs >= GEN.max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* area below hat infinite --> insert new construction point.
       We have to find out on which side of the construction point
       the area of the hat is unbounded. */
    if (iv->Ahatr >= INFINITY) {
      /* right hand side */

      /* iv should never be the last (virtual) interval */
      CHECK_NULL(iv->next,0);

      /* use mean point between the construction point and the right
	 boundary of the interval.
	 (The right boundary point might be a better choice
	 but cannot be used in every case.) */
      x = _unur_arcmean(iv->x,iv->next->ip);
      fx = PDF(x);

      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */

      /* if fx is 0, then we can cut off the tail of the distribution
	 (since it must be T-concave)  */
      if (fx <= 0.) {
	if (iv->next->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  return UNUR_ERR_GEN_CONDITION;
	}

	/* cut off right tail */
	free(iv->next);
	--(GEN.n_ivs);
	iv->next = iv_new;
	iv_new->prev = iv;

      }

      else {
	/* insert the new interval into the linked list after the old one. */
	iv_new->prev = iv;
	iv_new->next = iv->next;
	iv->next->prev = iv_new;
	iv->next = iv_new;
      }
      /* the old interval has to be recomputed. */

    }

    else {
      /* left hand side */

      x = _unur_arcmean(iv->ip,iv->x);
      fx = PDF(x);

      iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
      if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* PDF(x) < 0 or overflow !! */


      /* if fx is 0, then we can cut off the tail of the distribution
	 (since it must be T-concave)  */
      if (fx <= 0.) {
	if (iv->fx > 0.) {
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	  return UNUR_ERR_GEN_CONDITION;
	}

	/* cut off left tail */
	iv_new->next = iv->next;
	iv_new->prev = iv->prev;
	iv_new->ip = iv->ip;
	iv_new->fip = iv->fip;
	--(GEN.n_ivs);
	GEN.iv = iv_new;
	/* continue with this new interval */
	free(iv);
	iv = iv_new;
      }

      else {

	if (iv->prev) {
	  /* insert new interval in just before the old unbounded one */
	  iv_tmp = iv->prev;
	  iv_new->prev = iv->prev;
	  iv_new->next = iv;
	  iv->prev->next = iv_new;
	  iv->prev = iv_new;
	  
	  /* make sure that _unur_arcmean(iv->ip,iv->x) is never out of range */
	  iv_new->ip = iv->ip;
	  
	  /* continue with the interval before the old one
	     (neccessary since it will change too). */
	  iv = iv_tmp;
	}
	else { /* iv->prev == NULL */
	  /* insert new interval as first entry in list */
	  iv_new->ip = iv->ip;
	  iv_new->fip = iv->fip;
	  iv_new->prev = NULL;
	  iv_new->next = iv;
	  iv->prev = iv_new;
	  GEN.iv = iv_new;
	  
	  /* continue with this new interval */
	  iv = iv_new;
	}
      }
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_starting_intervals() */

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
  if (_unur_FP_is_infinity(PAR.darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_tdr_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Atot += iv->Ahat;
    Asqueezetot += iv->Asqueeze;
  }

  GEN.Atotal = Atot;
  GEN.Asqueeze = Asqueezetot;

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

/*---------------------------------------------------------------------------*/

int
_unur_tdr_gw_dars( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling  (Gilks&Wild)           */
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
  struct unur_tdr_interval *iv, *iv_next;
  double Alimit;               /* threshhold value for splitting interval */
  double x0, x1;               /* boundary of interval */
  double xsp, fxsp;            /* splitting point in interval */
  double xAhatl, xAhatr, xAsqueeze;
  int rule;                    /* id for splitting rule that has been applied */
  int n_splitted;              /* count splitted intervals */
  int splitted;                /* result of splitting routine */

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  
  /* now split intervals */
  while ( (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) &&
	  (GEN.n_ivs < GEN.max_ivs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN.n_ivs > 1)
      Alimit = PAR.darsfactor * ( (GEN.Atotal - GEN.Asqueeze) / GEN.n_ivs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 
    
    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (iv = GEN.iv; iv->next != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
      
      /* do not exceed the maximum number of intervals */
      if (GEN.n_ivs >= GEN.max_ivs)
	break;

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if ((iv->Ahat - iv->Asqueeze) <= Alimit) 
	continue;  /* goto next interval */

      /* store pointer to next interval */
      iv_next = iv->next;

      /* boundary of interval */
      x0 = iv->x;
      x1 = iv->next->x;

      /* get splitting point */
      for (rule = PAR.darsrule; rule <= 3; rule++) {
	switch (rule) {
	case 1:   /* rule 1: expected value */
	  
	  if ( _unur_FP_is_minus_infinity(x0) ||
	       _unur_FP_is_infinity(x1) ||
	       _unur_FP_approx(x0,x1) )
	    /* we do not use the expected value in case of unbounded intervals */
	    continue;  /* try next rule */
	  
	  /* l.h.s. hat */
	  xAhatl = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);
	  /* r.h.s. hat */
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv->next, iv->next->dTfx, iv->ip);
	  /* squeeze */
	  if (iv->Asqueeze > 0.)
	    /* always integrate from point with greater value of transformed density
	       to the other point */
	    xAsqueeze = (iv->Tfx > iv->next->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, iv->sq, x1)
	      : _unur_tdr_interval_xxarea( gen, iv->next, iv->sq, x0);
	  else  /* there is no squeeze */
	    xAsqueeze = 0.;
	  
	  /* check results */
	  if ( _unur_FP_is_infinity(xAhatl) || _unur_FP_is_minus_infinity(xAhatl) ||
	       _unur_FP_is_infinity(xAhatr) || _unur_FP_is_minus_infinity(xAhatr) ||
	       _unur_FP_is_infinity(xAsqueeze) || _unur_FP_is_minus_infinity(xAsqueeze) ||
	       _unur_FP_equal(iv->Ahat,iv->Asqueeze) )
	    continue;  /* try next rule */
	  
	  /* compute expected value */
	  xsp = (xAhatl+xAhatr-xAsqueeze) / (iv->Ahat - iv->Asqueeze);
	  break;
	  
	case 2:   /* rule 2: arcmean */
	  xsp = _unur_arcmean(x0,x1);
	  break;
	  
	case 3:  /* rule 3: mean */
	  if (_unur_FP_is_minus_infinity(x0) || _unur_FP_is_infinity(x1))
	    /* no arithmetic mean for unbounded intervals */
	    continue;  /* try next rule */
	  xsp = 0.5 * (x0 + x1);
	  break;
	  
	default:   /* this should not happen */
	  continue;
	}
	
	/* value of PDF at splitting point */
	fxsp = PDF(xsp);
	
	/* now split interval at given point */
	splitted = _unur_tdr_gw_interval_split(gen, iv, xsp, fxsp);

	if (splitted==UNUR_SUCCESS) {
	  /* splitting successful */
	  ++n_splitted;

	  /* now depending on the location of xps in the interval iv,
	     iv points to the left of the two new intervals,
	     or to the right of the two new intervals.
	     For the first case we have to move the pointer to the
	     new right interval. Then iv will be moved to the next
	     old interval in the list by the for loop.
	     We can distinguish between these two cases by looking 
	     at the iv->next pointer and compare it the pointer in the
	     old unsplitted interval.                                  */
	  if (iv->next != iv_next)
	    iv = iv->next;
	  /* no more splitting points in this interval */
	  break;
	}
	else if (splitted!=UNUR_ERR_SILENT) {
	  /* some serious error occurred */
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  return UNUR_ERR_GEN_CONDITION;
	}
	/* else: could not split construction points: too close (?) */
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
  if ( GEN.max_ratio * GEN.Atotal > GEN.Asqueeze ) {
    if ( GEN.n_ivs >= GEN.max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN.max_ivs = GEN.n_ivs;
  }
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_dars() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_ps_dars( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling  (Gilks&Wild)           */
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
  struct unur_tdr_interval *iv, *iv_next;
  double Adiff;                /* area between hat and squeeze */
  double Alimit;               /* threshhold value for splitting interval */
  double squeeze_gw;           /* slope of squeeze in variant GW */
  double Asqueeze_gw;          /* area below squeeze in variant GW */
  double Ahat_gw;              /* area below hat in variant GW */
  double x0, x1;               /* boundary of interval */
  double xsp, fxsp;            /* splitting point in interval */
  double xAhatl, xAhatr, xAsqueeze_gw;
  int rule;                    /* id for splitting rule that has been applied */
  int n_splitted;              /* count splitted intervals */
  int splitted;                /* result of splitting routine */
  
  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_TDR_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);

  /* split intervals */
  while ( (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) &&
	  (GEN.n_ivs < GEN.max_ivs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN.n_ivs > 1)
      Alimit = PAR.darsfactor * ( (GEN.Atotal - GEN.Asqueeze) / GEN.n_ivs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 
    
    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
      COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);

      /* do not exceed the maximum number of intervals */
      if (GEN.n_ivs >= GEN.max_ivs)
	break;

      /* boundary of interval:
	 in opposition to GW the construction point is in the middle of
	 an interval. The boundary points are the intersection points of
	 the tangents. Notice that it does not make sense to have a
	 splitting point near the construction point in the middle of
	 the interval.
	 Thus we use the same intervals as in GW, i.e., the intervals
	 between the construction points. However, now we have to 
	 compute some data and the first and last interval has to 
	 be treated differently.
      */
      if (iv==GEN.iv) {
	/* the first interval */
	x0 = iv->ip;       /* left boundary of interval */
	x1 = iv->x;        /* right boundary of interval */
	Adiff = iv->Ahat - iv->Ahatr; /* area below hat */
	Adiff *= 1. - iv->sq;         /* area between hat and squeeze */
      }
      else {
	/* all the other intervals:
	   we use the intervals between the preceeding construction point
	   and the construction point in the given interval for
	   finding a splitting point.                                */
	x0 = iv->prev->x;       /* left boundary of interval */
	x1 = iv->x;             /* right boundary of interval */
	/* area between hat and squeeze.
	   (notice: we habe two parts in two intervals.)      */
	Adiff = ( (1. - iv->prev->sq) * iv->prev->Ahatr +
		  (1. - iv->sq) * (iv->Ahat - iv->Ahatr) );
      }

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if (Adiff <= Alimit) 
	continue;  /* goto next interval */

      /* get splitting point */
      for (rule = PAR.darsrule; rule <= 3; rule++) {
	switch (rule) {
	case 1:   /* rule 1: expected value */
	  
	  if ( _unur_FP_is_minus_infinity(x0) ||
	       _unur_FP_is_infinity(x1) ||
	       _unur_FP_approx(x0,x1) )
	    /* we do not use the expected value in case of unbounded intervals */
	    continue;  /* try next rule */

	  /* l.h.s. hat */
	  xAhatl = ( (iv->prev == NULL) ? 0. : 
		     _unur_tdr_interval_xxarea( gen, iv->prev, iv->prev->dTfx, iv->ip) );
	  /* r.h.s. hat */
	  xAhatr = _unur_tdr_interval_xxarea( gen, iv, iv->dTfx, iv->ip);

	  /* area below hat */
	  Ahat_gw = ( ((iv->prev == NULL) ? 0. : iv->prev->Ahatr)
		      + (iv->Ahat - iv->Ahatr) );

	  /* squeeze */
	  if (iv->Asqueeze > 0. && iv->prev != NULL) {
	    /* slope of transformed squeeze for variant GW */
	    squeeze_gw = (iv->Tfx - iv->prev->Tfx) / (iv->x - iv->prev->x);

	    /* always integrate from point with greater value of transformed density
	       to the other point */
	    xAsqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_xxarea( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_xxarea( gen, iv->prev, squeeze_gw, iv->x);
	    Asqueeze_gw = (iv->Tfx > iv->prev->Tfx)
	      ? _unur_tdr_interval_area( gen, iv, squeeze_gw, iv->prev->x)
	      : _unur_tdr_interval_area( gen, iv->prev, squeeze_gw, iv->x);
	  }
	  else { /* there is no squeeze */
	    xAsqueeze_gw = 0.;
	    Asqueeze_gw = 0.;
	  }

	  /* check results */
	  if ( _unur_FP_is_infinity(xAhatl) || _unur_FP_is_minus_infinity(xAhatl) ||
	       _unur_FP_is_infinity(xAhatr) || _unur_FP_is_minus_infinity(xAhatr) ||
	       _unur_FP_is_infinity(xAsqueeze_gw) || _unur_FP_is_minus_infinity(xAsqueeze_gw) ||
	       _unur_FP_is_infinity(Asqueeze_gw)  || _unur_FP_is_minus_infinity(Asqueeze_gw) ||
	       _unur_FP_equal(Ahat_gw,Asqueeze_gw) )
	    continue;  /* try next rule */

	  /* compute expected value */
	  xsp = (xAhatl+xAhatr-xAsqueeze_gw) / (Ahat_gw - Asqueeze_gw);
	  break;
	  
	case 2:   /* rule 2: arcmean */
	  xsp = _unur_arcmean(x0,x1);
	  break;
	  
	case 3:  /* rule 3: mean */
	  if (_unur_FP_is_minus_infinity(x0) || _unur_FP_is_infinity(x1))
	    /* no arithmetic mean for unbounded intervals */
	  continue;  /* try next rule */
	  xsp = 0.5 * (x0 + x1);
	  break;
	  
	default:   /* this should not happen */
	  continue;
	}
	
	/* value of PDF at splitting point */
	fxsp = PDF(xsp);

	/* store pointer to next interval */
	iv_next = iv->next;

	/* now split interval at given point.
	   we have to take care that xsp might be either in
	   interval iv or in its preceeding interval, depending
	   whether xsp is less than iv->ip or not. */

	splitted = _unur_tdr_ps_interval_split(gen, ((xsp<iv->ip)?iv->prev:iv), xsp, fxsp);

	if (splitted==UNUR_SUCCESS) {
	  /* splitting successful */
	  ++n_splitted;

	  /* now depending on the location of xps in the interval iv,
	     iv points to the left of the two new intervals,
	     or to the right of the two new intervals.
	     For the first case we have to move the pointer to the
	     new right interval. Then iv will be moved to the next
	     old interval in the list by the for loop.
	     We can distinguish between these two cases by looking 
	     at the iv->next pointer and compare it the pointer in the
	     old unsplitted interval.                                  */
	  if (iv->next != iv_next)
	    iv = iv->next;
	  /* no more splitting points in this interval */
	  break;
	}
	else if (splitted!=UNUR_ERR_SILENT) {
	  /* some serious error occurred */
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  return UNUR_ERR_GEN_CONDITION;
	}
	/* else: could not split construction points: too close (?) */
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
  if ( GEN.max_ratio * GEN.Atotal > GEN.Asqueeze ) {
    if ( GEN.n_ivs >= GEN.max_ivs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN.max_ivs = GEN.n_ivs;
  }
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_dars() */

/*****************************************************************************/

static struct unur_tdr_interval *
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
  iv = _unur_malloc( sizeof(struct unur_tdr_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN.n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_TDR_IV);

  /* avoid uninitialized variables */
  iv->Acum = iv->Ahat = iv->Ahatr = iv->Asqueeze = 0.;
  
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
    /*      iv->Tfx = -pow(fx,GEN.c_T); */
    /*      iv->dTfx = 0.; */
    break;
  }

  /* the program requires dTfx > -INFINITY */
  if (iv->dTfx <= -INFINITY)
    iv->dTfx = INFINITY;

  return iv;

} /* end of _unur_tdr_interval_new() */

/*****************************************************************************/

static int
_unur_tdr_gw_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (Gilks & Wild variant)                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... do not add this construction point             */
     /*   UNUR_ERR_INF    ... area = INFINITY                                */
     /*   others          ... error (PDF not T-concave)                      */
     /*----------------------------------------------------------------------*/
{
  double Ahatl;    /* area below hat at left side of intersection point */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents.
     used to partition interval into left hand part (construction point of tangent
     on the left hand boundary) and right hand part (construction point of tangent
     on the left hand boundary). */
  if ( _unur_tdr_tangent_intersection_point(gen,iv,&(iv->ip))!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;

  /* squeeze and area below squeeze */
  if (iv->Tfx > -INFINITY && iv->next->Tfx > -INFINITY) {

    /* we do not compute the slope when the construction points
       are too close. at least 8 significant digits should remain. */
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   /* construction points too close */

    /* slope of transformed squeeze */
    iv->sq = (iv->next->Tfx - iv->Tfx) / (iv->next->x - iv->x);

    /* check squeeze */
    /* we have to take care about round off error.
       the following accepts PDFs with might be a little bit not T_concave */
    if ( ( (iv->sq > iv->dTfx       && !_unur_FP_approx(iv->sq,iv->dTfx)) || 
	   (iv->sq < iv->next->dTfx && !_unur_FP_approx(iv->sq,iv->next->dTfx)) )
	 && iv->next->dTfx < INFINITY ) {
      /* There are big troubles when the density is extremely small. 
	 Then round-off errors may cancel out all significant figures and
	 0 remains. Thus we simply ignore all violations when the 
	 slope of the squeeze or tangent is 0.  */
      if ( iv->sq != 0. && iv->dTfx != 0. && iv->next->dTfx != 0. ) {
      	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      	return UNUR_ERR_GEN_CONDITION;
      }
    }

    /* volume below squeeze */
    /* always integrate from point with greater value of transformed density
       to the other point */
    iv->Asqueeze = (iv->Tfx > iv->next->Tfx)
      ? _unur_tdr_interval_area( gen, iv, iv->sq, iv->next->x)
      : _unur_tdr_interval_area( gen, iv->next, iv->sq, iv->x);
  }
  else {  /* no squeeze */
    iv->sq = 0.;
    iv->Asqueeze = 0.;
  }

  /* volume below hat */
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv->next, iv->next->dTfx, iv->ip);

  /* areas below head unbounded ? */
  if (_unur_FP_is_infinity(Ahatl) || _unur_FP_is_infinity(iv->Ahatr))
    return UNUR_ERR_INF;

  /* total area */
  iv->Ahat = iv->Ahatr + Ahatl;

  /* check area */
  /* we cannot be more accurate than in the `check squeeze' section */
  if ( iv->Asqueeze > iv->Ahat && !_unur_FP_approx(iv->Asqueeze, iv->Ahat) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"A(squeeze) > A(hat). PDF not T-concave!");
    return UNUR_ERR_GEN_CONDITION; 
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_interval_parameter() */

/*---------------------------------------------------------------------------*/

static int
_unur_tdr_ps_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (proportional squeezes)                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... do not add this construction point             */
     /*   UNUR_ERR_INF    ... area = INFINITY                                */
     /*   others          ... error (PDF not T-concave)                      */
     /*----------------------------------------------------------------------*/
{
  double Ahatl;    /* area below hat at left side of intersection point */
  double hxl, hxr; /* value of hat at left and right point of interval */
  double sq;       /* ration PDF(x) / hat(x) */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents. it is used as boundaries of intervals.
     it is stored together with the right hand (i.e. next) construction point. */
  if (_unur_tdr_tangent_intersection_point(gen,iv,&(iv->next->ip))!=UNUR_SUCCESS)
    return UNUR_ERR_GEN_CONDITION;
  /* value of PDF at intersection point */
  iv->next->fip = _unur_FP_is_infinity(iv->next->ip) ? 0. : PDF(iv->next->ip);

  /* volume below hat */
  Ahatl = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->ip);
  iv->Ahatr = _unur_tdr_interval_area( gen, iv, iv->dTfx, iv->next->ip);

  /* areas below head unbounded ? */
  if (Ahatl >= INFINITY || iv->Ahatr >= INFINITY)
    return UNUR_ERR_INF;

  /* total area */
  iv->Ahat = iv->Ahatr + Ahatl;

  /* compute squeeze:
     squeeze ration = min_{boundary points} PDF(x) / hat(x) */
  
  /* left boundary point */
  hxl = _unur_tdr_eval_intervalhat(gen,iv,iv->ip);
  if (_unur_FP_greater(iv->fip, hxl) ) {
    /* PDF(x) > hat(x); this should not happen */
    if ( (iv->fip < 1.e-50) || _unur_FP_approx(iv->fip, hxl)) {
      /* hat(x) and PDF(x) are approximatly the same, or extremely small.
	 assume round-off error */
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      /* we really have PDF(x) > hat(x); we do not assume a round-off error */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  iv->sq = (_unur_FP_is_infinity(hxl) || hxl <= 0.) ? 0. : iv->fip / hxl;

  /* right boundary point */
  hxr = _unur_tdr_eval_intervalhat(gen,iv,iv->next->ip);
  if (_unur_FP_greater(iv->next->fip, hxr)) {
    /* PDF(x) > hat(x); this should not happen */
    if ((iv->next->fip < 1.e-50) || _unur_FP_approx(iv->next->fip, hxr)) {
      /* hat(x) and PDF(x) are approximatly the same, or extremely small.
	 assume round-off error */
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) might be < PDF(x)");
    }
    else {
      /* we really have PDF(x) > hat(x); we do not assume a round-off error */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"hat(x) < PDF(x)");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  sq = (_unur_FP_is_infinity(hxr) || hxr <= 0.) ? 0. : iv->next->fip / hxr;

  /* squeeze */
  if (iv->sq > sq) iv->sq = sq;

  /* area below squeeze */
  iv->Asqueeze = iv->Ahat * iv->sq;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_interval_parameter() */

/*****************************************************************************/

static int
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

static double
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

static double
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

static double
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

static int
_unur_tdr_gw_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv_oldl, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* split interval iv_oldl into two intervals at point x                 */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv_oldl ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   fx      ... value of PDF at x                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv_newr;  /* pointer to new interval */
  struct unur_tdr_interval iv_bak;    /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_TDR_IV,UNUR_ERR_COOKIE);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_gw_debug_split_start( gen,iv_oldl,x,fx );
#endif

  /* the splitting point must be inside the interval */
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN.n_ivs * (iv_oldl->Ahat - iv_oldl->Asqueeze) / (GEN.Atotal - GEN.Asqueeze))
       < GEN.bound_for_adding)
    return UNUR_ERR_SILENT;

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }

  /* back up data */
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_tdr_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {
    
    /* one of the two boundary points must be 0, too! */
    if (iv_oldl->fx <= 0.) {
      /* chop off left part (it's out of support) */
      iv_oldl->x = x;
    }
    else if (iv_oldl->next->fx <= 0.) {
      /* chop off right part (it's out of support) */
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* compute parameters for chopped interval */
    success = _unur_tdr_gw_interval_parameter(gen, iv_oldl);
    
    /* we did not add a new interval */
    iv_newr = NULL;
  }

  else {
    
    /* we need a new interval */
    iv_newr = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_newr == NULL) {
      /* PDF(x) < 0 or overflow !! */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* insert into linked list */
    iv_newr->prev = iv_oldl;
    iv_newr->next = iv_oldl->next;
    iv_oldl->next->prev = iv_newr;
    iv_oldl->next = iv_newr;
    
    /* compute parameters for interval */
    success   = _unur_tdr_gw_interval_parameter(gen, iv_oldl);
    success_r = _unur_tdr_gw_interval_parameter(gen, iv_newr);
    
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  
  /* successfull ? */
  if (success!=UNUR_SUCCESS) {
    /* cannot split interval at given point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");

    /* the case of unbounded hat is treated as round-off error for 
       very steep tangents. so we simply do not add this construction point. */

    /* restore old interval */
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_tdr_interval));
    /* remove from linked list; remaines to restore prev pointer in next interval */
    if (iv_oldl->next)
      iv_oldl->next->prev = iv_oldl;

    /* decrement counter for intervals and free unused interval */
    if (iv_newr) {
      --(GEN.n_ivs); 
      free( iv_newr );
    }

  return ( (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF) 
	   ? UNUR_ERR_GEN_CONDITION : UNUR_SUCCESS );
  }

  /* successful */

  /* update total area below hat and squeeze */
  GEN.Atotal   = ( GEN.Atotal - iv_bak.Ahat
		   + iv_oldl->Ahat + ((iv_newr) ? iv_newr->Ahat : 0.) );
  GEN.Asqueeze = ( GEN.Asqueeze - iv_bak.Asqueeze
		   + iv_oldl->Asqueeze + ((iv_newr) ? iv_newr->Asqueeze : 0. ) );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug & TDR_DEBUG_SPLIT)
      _unur_tdr_gw_debug_split_stop( gen,iv_oldl,iv_newr );
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_gw_interval_split() */

/*---------------------------------------------------------------------------*/

static int
_unur_tdr_ps_interval_split( struct unur_gen *gen, struct unur_tdr_interval *iv, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* split interval iv into two intervals at point x                      */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv      ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   fx      ... value of PDF at x                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *oldl, *oldr;  /* pointer to old intervals (left, right) */
  struct unur_tdr_interval *iv_new;       /* pointer to new interval */
  struct unur_tdr_interval oldl_bak, oldr_bak; /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL); COOKIE_CHECK(gen,CK_TDR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);  COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);

  /* we only add a new construction point, if the relative area is large enough */
  if ( (GEN.n_ivs * (iv->Ahat - iv->Asqueeze) / (GEN.Atotal - GEN.Asqueeze))
       < GEN.bound_for_adding)
    return UNUR_ERR_SILENT;

  /* the splitting point must be inside the interval */
  if (x < iv->ip || x > iv->next->ip) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }

  /* which side of construction point */
  if (x < iv->x) {
    /* left hand side */
    oldl = iv->prev;
    oldr = iv;
  }
  else {
    /* right hand side */
    oldl = iv;
    oldr = iv->next;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_start( gen,oldl,oldr,x,fx );
#endif

  /* back up data */
  if (oldl) memcpy(&oldl_bak, oldl, sizeof(struct unur_tdr_interval));
  memcpy(&oldr_bak, oldr, sizeof(struct unur_tdr_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (fx <= 0.) {

    /* one of the two boundary points must be 0, too! */
    if (oldr->fip <= 0. && oldl==NULL) {
      /* chop off left part (it's out of support) */
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else if (oldr->fip <= 0. && oldr->next==NULL) {
      /* chop off right part (it's out of support) */
      oldr->x = x;
      oldr->ip = x;
      oldr->fip = 0.;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* we do not add a new interval */
    iv_new = NULL;
  }

  else {

    /* we need a new interval */
    iv_new = _unur_tdr_interval_new( gen, x, fx, FALSE );
    if (iv_new == NULL) {
      /* PDF(x) < 0 or overflow !! */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return -1;
    }

    /* insert into linked list */
    iv_new->prev = oldl;
    iv_new->next = oldr;
    oldr->prev = iv_new;
    if (oldl) oldl->next = iv_new;
    
  }

  /* compute parameters for intervals */
  success = UNUR_SUCCESS;
  /* left hand interval */
  if (oldl) {
   success = _unur_tdr_ps_interval_parameter(gen, oldl);
  }
  if( iv_new ) {
    /* middle (newly created) interval) */
    if (!oldl) {
      /* we have to copy the left intersection point from the 
         right hand interval */
      iv_new->ip = oldr->ip;
      iv_new->fip = oldr->fip;
    }
    success_r = _unur_tdr_ps_interval_parameter(gen, iv_new);
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }
  if ( oldr->next ) {
    /* right hand interval */
    success_r = _unur_tdr_ps_interval_parameter(gen, oldr);
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }

  /* successfull ? */
  if (success!=UNUR_SUCCESS) {
    /* cannot split interval at given point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");

    /* the case of unbounded hat is treated as round-off error for 
       very steep tangents. so we simply do not add this construction point. */

    /* restore old interval */
    if (oldl) memcpy(oldl, &oldl_bak, sizeof(struct unur_tdr_interval));
    memcpy(oldr, &oldr_bak, sizeof(struct unur_tdr_interval));
    /* remove from linked list; remaines to restore prev pointer in next interval */
    oldr->prev = oldl;
    if (oldl) oldl->next = oldr;

    /* decrement counter for intervals and free unused interval */
    if (iv_new) {
      --(GEN.n_ivs); 
      free( iv_new );
    }

  return ( (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF) 
	   ? UNUR_ERR_GEN_CONDITION : UNUR_SUCCESS );
  }

  /* successful */

  /* we have update the pointer to the list */
  if (oldl == NULL && iv_new)
    /* new first entry */
    GEN.iv = iv_new;

  /* update total area below hat and squeeze */
  GEN.Atotal   = ( GEN.Atotal + (oldr->Ahat - oldr_bak.Ahat)
		   + ((oldl) ? (oldl->Ahat - oldl_bak.Ahat) : 0.)
		   + ((iv_new) ? iv_new->Ahat : 0.) );
  GEN.Asqueeze = ( GEN.Asqueeze + (oldr->Asqueeze - oldr_bak.Asqueeze)
		   + ((oldl) ? (oldl->Asqueeze - oldl_bak.Asqueeze) : 0.)
		   + ((iv_new) ? iv_new->Asqueeze : 0.) );

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & TDR_DEBUG_SPLIT) 
    _unur_tdr_ps_debug_split_stop( gen,oldl,iv_new,oldr );
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_tdr_ps_interval_split() */

/*****************************************************************************/

static int
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
  if (!GEN.guide) {
    int max_guide_size = (GEN.guide_factor > 0.) ? (GEN.max_ivs * GEN.guide_factor) : 1;
    GEN.guide = _unur_malloc( max_guide_size * sizeof(struct unur_tdr_interval*) );
  }

  /* first we need cumulated areas in intervals */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN.Atotal = Acum;
  GEN.Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN.guide_size = (int)(GEN.n_ivs * GEN.guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dis.c) */
  Astep = GEN.Atotal / GEN.guide_size;
  Acum=0.;
  for( j=0, iv=GEN.iv; j < GEN.guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TDR_IV,UNUR_ERR_COOKIE);
    while( iv->Acum < Acum )
      iv = iv->next;
    if( iv->next == NULL ) {   /* this is the last virtual intervall --> do not use */
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN.guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide[j] = iv;

  return UNUR_SUCCESS;
} /* end of _unur_tdr_make_guide_table() */

/*****************************************************************************/

