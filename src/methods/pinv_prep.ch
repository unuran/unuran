/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_prep.c                                                  *
 *                                                                           *
 *   Routines for preprocessing steps.                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/** Preprocessing                                                           **/
/*****************************************************************************/

int
_unur_pinv_relevant_support ( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* 1a. Estimate computationally relevant domain (support) of PDF        */
     /*     (finite interval where PDF is above some threshold value).       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double fb;

  /* check PDF at boudary */
  if(GEN->sleft) {
    /* compute PDF at boundary */
    fb = PDF(GEN->dleft);
    if (fb > 1.e-20 && fb < 1.e300) {
      /* PDF does neither vanish nor is it infinity:      */
      /* use boundary of given domain of distribution for */
      /* computational domain.                            */
      GEN->bleft = GEN->dleft;        
      /* there is no need for further searching           */
      GEN->sleft = FALSE;
    }
    /** FIXME: throw error if fb is not finite! **/
  }
 
  /* same for right hand boundary */
  if(GEN->sright) {
    fb = PDF(GEN->dright);
    if (fb > 1.e-20 && fb < 1.e300) {
      GEN->bright = GEN->dright;        
      GEN->sright = FALSE;
    }
  }

  /* search for interval of computational relevance (if required) */
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_searchborder(gen,DISTR.center, GEN->bleft, 
					 &(GEN->dleft), &(GEN->sleft) );
    if (!_unur_isfinite(GEN->bleft)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get left boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  
  if(GEN->sright) {
    GEN->bright = _unur_pinv_searchborder(gen,DISTR.center, GEN->bright,
					  &(GEN->dright), &(GEN->sright) );
    if (!_unur_isfinite(GEN->bright)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get right boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_relevant_support(gen);
#endif

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_relevant_support() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound,
			 double *dom, int *search)
     /*----------------------------------------------------------------------*/
     /* [1a.] Find left or right hand border of relevant domain.             */
     /*                                                                      */
     /* Calculate domain of computational relevant region.                   */
     /* Start at 'x0' and search towards 'bound'.                            */
     /* The boundary points of this domain are approximately given as        */
     /*      PDF(x0) * PINV_PDFLLIM                                          */
     /*                                                                      */
     /* As a side effect the support of the distribution is shrinked if      */
     /* points with PDF(x)=0 are found.                                      */
     /* If in addition a discontinuity is detected then the exact position   */
     /* (up to machine precision) of the boundary of the support of the PDF  */
     /* is located via interval bisection.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   x0     ... starting point for searching boudary                    */
     /*              PDF(x0) must not be too small                           */
     /*   bound  ... stop searching at this point                            */
     /*   dom    ... pointer to boundary of domain / support of distribution */
     /*   search ... pointer to boolean that indicates whether we have to    */
     /*              for cut-off points. The boolean is set to FALSE if a    */
     /*              discontinuity is found at the boundary.                 */
     /*                                                                      */
     /* return:                                                              */
     /*   boundary point                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double x;         /* current and previous searching point */
  double xs, xl;       /* point where PDF is less than and larger than threshold */
  double fx;              /* PDF at x */
  double fs, fl;          /* PDF at xs and xl */
  double fllim;           /* threshold value */
  double fulim;           /* threshold for detecting discontinuity at boundary */

  /* threshold value where we stop searching */
  fllim = PDF(x0) * PINV_PDFLLIM;
  fulim = 1.e4 * fllim;

  /* starting point */
  xl = x0; 
  fl = INFINITY;
  x = _unur_arcmean(x0,bound);

  /* find a points where PDF values bracket threshold: */
  /*   fs = PDF(xs) <= fllim <= PDF(xl) = fl           */
  while ( (fx=PDF(x)) > fllim ) {
    if (_unur_FP_same(x,bound))
      return bound;
    xl = x; fl = fx;
    x = _unur_arcmean(x,bound);
  }
  xs = x; fs = fx;

  /* decrease length of bracket if necessary */
  while (!_unur_FP_same(xs,xl)) {

    /* truncate domain if possible */
    if (_unur_iszero(fs)) {
      *dom = xs;
    }

    /* new point */
    x = xs/2. + xl/2.;
    fx = PDF(x);

    /* check PDF at new point */
    if (fx < fllim) {
      /* update bracket */
      xs = x; fs = fx;
    }
    else {
      if (fl > fulim) {
	/* assume discontinuity at boundary -> update bracket */
	xl = x; fl = fx;
      }
      else {
	/* assume smooth PDF -> stop and return point */
	return x;      
      }
    }
  }

  /* since we have already found boundary point, we can use this       */
  /* for the tail cut-off point. Thus we switch off further searching. */
  *search = FALSE;

  /* return point */
  return xl;
  
} /* end of _unur_pinv_searchborder() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_approx_pdfarea (struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* 1b. Compute area below PDF over relevant domain approximately.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* Remark:                                                              */
     /*   The user can also provide the area below the PDF.                  */
     /*   However, then we probably need not method PINV                     */
     /*----------------------------------------------------------------------*/
{
  double tol;   /* tolerated integration error */
  int i;        /* number of trials            */
  int res = UNUR_SUCCESS; /* exit code of subroutine */

  /* we might need two trials */
  for (i=1; i<=2; i++) {

    /* we only need a rough approximation of the area below the PDF. */
    /* for the tolerated absolute integration error we assume that   */
    /* it is not "too far" from 1 or the value provides by the user  */
    /* (if given), respectively.                                     */
    tol = PINV_UERROR_AREA_APPROX * GEN->area;

    /* check center of distribution */
    DISTR.center = _unur_max(DISTR.center, GEN->dleft);
    DISTR.center = _unur_min(DISTR.center, GEN->dright);

    /* compute area */
    GEN->area  = _unur_pinv_adaptivelobatto5( gen, GEN->bleft,   DISTR.center - GEN->bleft,
					      tol, NULL );
    GEN->area += _unur_pinv_adaptivelobatto5( gen, DISTR.center, GEN->bright - DISTR.center,
					      tol, NULL );

    /* check estimated area */
    if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
      res = UNUR_FAILURE;
      break;
    }
    
    /* The tolerated integration error is relative to the PDF area.
     * For the first trial we have used area=1 (or the value given by the user) 
     * as a first guess.
     * However, we have to check whether this guess is large enough, 
     * since otherwise the relative integration error is too large.
     */
    if (GEN->area > 1.e-2) {
      /* relative integration error is sufficiently small */
      break;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_pdfarea(gen,TRUE);
#endif

  return res;

} /* end of _unur_pinv_approx_pdfarea() */

/*---------------------------------------------------------------------------*/
#ifdef PINV_USE_CDFTABLE
/*---------------------------------------------------------------------------*/

int
_unur_pinv_pdfarea (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* 1d. Compute area below PDF with requested accuracy and               */
     /*     store intermediate results from adaptive integration.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* Remark:                                                              */
     /*   The user can also provide the area below the PDF.                  */
     /*   However, then we probably need not method PINV                     */
     /*----------------------------------------------------------------------*/
{
  double tol;   /* tolerated integration error */

  /* the tolerated U-error depends on the given u-resolution.        */
  /* however, we have to make some adjustments for integration error */
  /* in order to achieve the desired precison goal.                  */
  tol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION * PINV_UTOL_CORRECTION;

  /* check center of distribution */
  DISTR.center = _unur_max(DISTR.center, GEN->dleft);
  DISTR.center = _unur_min(DISTR.center, GEN->dright);

  /* before we compute the area we have to reset the array for storing CDF values */
  if (GEN->CDFtable)
    _unur_pinv_CDFtable_append(GEN->CDFtable,GEN->bleft,0.);
  
  /* compute area. the ordering of the two integrals is important */
  GEN->area  = _unur_pinv_adaptivelobatto5( gen, GEN->bleft,   DISTR.center - GEN->bleft,
					    tol, GEN->CDFtable );
  GEN->area += _unur_pinv_adaptivelobatto5( gen, DISTR.center, GEN->bright - DISTR.center,
					    tol, GEN->CDFtable );

  /* check estimated area */
  if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
    return UNUR_FAILURE;
  }

  /* reallocate memory for table of CDF values */
  _unur_pinv_CDFtable_resize(&(GEN->CDFtable));

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_pdfarea(gen,FALSE);
#endif

  return UNUR_SUCCESS;

} /* end of _unur_pinv_pdfarea() */

/*---------------------------------------------------------------------------*/
#endif /* defined(PINV_USE_CDFTABLE) */
/*---------------------------------------------------------------------------*/

int
_unur_pinv_computational_domain (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* 1c. Compute computational domain where inverse CDF is approximated   */
     /*     (interval where we safely can compute coefficients of            */
     /*     interpolating polynomial).                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double tailcut_error;    /* threshold values for cut-off points */
  double range;            /* length of current working domain */

  /* parameters for tail cut-off points: maximal area in tails          */
  /* We use the given U-reslution * PINV_TAILCUTOFF_FACTOR * PDFarea.   */
  tailcut_error = GEN->u_resolution * PINV_TAILCUTOFF_FACTOR(GEN->u_resolution);
  tailcut_error = _unur_min( tailcut_error, PINV_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PINV_UERROR_CORRECTION;

  /* length of current working domain */
  range = GEN->bright-GEN->bleft;

  /* compute cut-off points for tails */
  if(GEN->sleft) {
    GEN->bleft = _unur_pinv_cut( gen, GEN->dleft, GEN->bleft, -range, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

  if(GEN->sright) {
    GEN->bright = _unur_pinv_cut( gen, GEN->dright, GEN->bright, range, tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_computational_domain(gen);
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_computational_domain() */

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#if 1
/*---------------------------------------------------------------------------*/

/* new version */

double
_unur_pinv_cut( struct unur_gen *gen, double dom, double w, double dw, double crit )
     /*----------------------------------------------------------------------*/
     /* [1c.] Calculate cut-off points for computational domain of           */
     /* distribution.                                                        */
     /* The area outside the cut-off point is given by 'crit'.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   dom  ... boundary of domain / support of distribution              */
     /*   w    ... starting point for searching cut-off point                */
     /*   dw   ... initial step size for searching,                          */
     /*            sign of dw gives searching direction:                     */
     /*               dw < 0 ... left hand side cut-off point                */
     /*               dw > 0 ... right hand side cut-off point               */
     /*   crit ... u-error criterium for tail cut off                        */
     /*                                                                      */
     /* return:                                                              */
     /*   cut-off point                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double sgn = (dw>0) ? 1. : -1.; /* searching direction */

  double fl,fx,fr;  /* value of PDF at x-dx, x, and x+dx */
  double x = w;     /* current point */
  double dx;        /* step length for calculation of derivative and local concavity */                                   
  double xnew;      /* new point in iteration */

  double df;        /* estimate for derivative of PDF at x */
  double lc;        /* estimate for local concavity of PDF at x */
  double area;      /* estimate for tail probability */

  int i,j;          /* auxiliary variables */

  /* check length of interval */
  if (_unur_iszero(fabs(dw))) return w;

  /* starting point and step size for search */
  x = w;

  /* iteratively search for cut-off point with tail probability approximately 'crit'ical value */
  for (i=1; i<100; i++) {

    /* we need step size 'dx' for computing derivative 'df' and local concavity 'lc' */

    /* first try */
    dx = (fabs(dw) + fabs(x-w)) * 1.e-3;

    /* check boundary of domain */
    if (x-dx < GEN->dleft)  dx = x - GEN->dleft;
    if (x+dx > GEN->dright) dx = GEN->dright - x;

    /* now let us try to find non-zero points for the PDF */
    for (j=1;;j++) {

      /* decrease step size */
      dx = dx/2.;
      
      /* check length and protect against infinite loops */
      if (dx < 128.*DBL_EPSILON*fabs(dw)) {
	/* we are too close to the boundary. So we just return the last value. */
	return x;
      }

      /* compute PDF values */
      fx = PDF(x);
      fl = PDF(x-dx);
      fr = PDF(x+dx);

      /* check values */
      if (! (_unur_iszero(fl) || _unur_iszero(fx) ||_unur_iszero(fr)) )
	break;
    }

    /* derivative */
    df = (fr-fl)/(2.*dx);

    /* local concavity */
    lc = fl/(fl-fx)+fr/(fr-fx) - 1;

    /* tail probability */
    area = fabs(fx*fx / ((lc+1.) * df));

    /* check results */
    if (_unur_isnan(area)) {
      /* When the PDF is extremely flat than NaN might occur. There are two possibilities: 
       * (1) We are far away from the center. Then the tail probabilty is almost 0.
       * (2) PDF(X) is still quite large. Then the tail probability cannot be estimated.
       * We assume (1).
       */
      _unur_warning(gen->genid,UNUR_ERR_NAN,"tail probability gives NaN --> assume 0.");
      return x;
    }

    if (sgn * df > 0.) {
      /* There is a maximum of the PDF in [fl,fr]. */
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not monotone at boundary");
      return x;
    }

    /* check accuracy of computation */
    if (fabs(area/crit-1.) < 1.e-4)
      return x;
    
    /* compute next point */
    xnew = x + fx/(lc*df) * ( pow(crit*fabs(df)*(lc+1.)/(fx*fx),lc/(lc+1.)) - 1.);

    /* check results */
    if (! _unur_isfinite(xnew)) {
      /* we cannot compute the next point */
      _unur_warning(gen->genid,UNUR_ERR_NAN,"numerical problems with cut-off point");
      return x;
    }
    
    if (sgn*dom < sgn*x) {
      /* boundary exceeded */
      return dom;
    }
    
    /* update point */
    x = xnew;
    
  }

  /* maximum number of iterations exceeded */
  return x;

} /* end of _unur_pinv_cut() */

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/

/* old version */

double
_unur_pinv_cut( struct unur_gen *gen, double dom, double w, double dw, double crit )
     /*----------------------------------------------------------------------*/
     /* [1c.] Calculate cut-off points for computational domain of           */
     /* distribution.                                                        */
     /* The area outside the cut-off point is given by 'crit'.               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   dom  ... boundary of domain / support of distribution              */
     /*   w    ... starting point for searching cut-off point                */
     /*   dw   ... length of working interval,                               */
     /*            sign of dw gives searching direction:                     */
     /*               dw < 0 ... left hand side cut-off point                */
     /*               dw > 0 ... right hand side cut-off point               */
     /*   crit ... u-error criterium for tail cut off                        */
     /*                                                                      */
     /* return:                                                              */
     /*   cut-off point                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double u, uplus;     /* tail probabilities */

  int found = FALSE;   /* indicates whether we have found an cut-off point */
  double dx;           /* step size for numeric differentiation */
  int j;               /* aux variable */

  double s = (dw>0) ? 1 : -1; /* searching direction */

  /* initial step size */
  dw /= 128.;

  /* step size for numeric differentiation */
  dx = dw/64.;

  /* search for cut-off point with tail probability less than 'crit'ical value */
  for (j=1; j<1000; j++) {
    
    /* check for boundary */
    if (s*dom < s*w) {
      /* boundary exceeded */
       return dom;
    }

    /* compute approximate tail probability at w */
    u = _unur_pinv_tailprob(gen, w, dw/64.);

    /* below threshold value ? */
    if (u < crit && u >= 0.) {
      found = TRUE;
      break;
    }

    /* else
       the tail probability is too large, or the approximation formula
       is not appropriate for the point.
       Hence we make a step towards + or - infinity ...
    */
    w += dw;

    /* ... and increase stepsize for next interation. */
    if (j>32) dw *= 1.5;
  }

  if(!found){
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cound not find valid cut point");
    return INFINITY;
  }


  if (_unur_iszero(u)) {
    /* the tail probability at w is approximately 0, i.e. _unur_pinv_tailprob() */
    /* assumes it is 0. Thus we only can return this point w.                   */
    return w;
  }

  /* now run secant method to find cut-off point with tail probability approx 'crit' */

  /* step size for numeric differentiation */
  dx=dw/64.;

  for (j=0; j<2048; j++) {

    /* check whether 'u' approx 'crit' */
    if (fabs(crit/u - 1.)<1.e-7) 
      return w;

    /* compute tail probability in the next point towards 'dx' */
    /* (required for numerical derivative) */
    uplus = _unur_pinv_tailprob(gen,w+dx,dx);
    if(uplus<0){
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"negative tail probability");
      return INFINITY;
    }
    if (_unur_iszero(uplus)) 
      return w;

    /* next iteration for secant method */
    w -= dx * (1./u-1./crit)/(1./uplus-1./u);

    /* check for boundary */
    if (s*dom < s*w)
      /* boundary exceeded */
       return dom;

    /* compute tail probability */
    u = _unur_pinv_tailprob(gen,w,dx);
    if(u<0){
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"negative tail probability");
      return INFINITY;
    }
    if (_unur_iszero(u))
      return w;
  }

  /* could not find point till now */
  _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"cound not find cut point");
  return INFINITY;

} /* end of _unur_pinv_cut() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_tailprob( struct unur_gen *gen, double x, double dx )
     /*----------------------------------------------------------------------*/
     /* [1c.] calculate approximate tail probability.                        */
     /* use formula                                                          */
     /*   area ~ f(x)^2 / ( (lc_f(x)+1) * abs(f'(x)) )                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... cut-off point                                              */
     /*   dx  ... step size for numeric differentiation                      */
     /*                                                                      */
     /* return:                                                              */
     /*   approximate tail probability                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return -INFINITY                                                   */
     /*----------------------------------------------------------------------*/
{
  double fx, fp, fm; /* value of PDF at x, x+dx, and x-dx */
  double df;         /* derivative of PDF                 */
  double lcplus1;    /* local concavity + 1               */
  double area;       /* tail probability                  */

  /* compute PDF */
  fx = PDF(x);
  fp = PDF(x+dx);
  fm = PDF(x-dx);

  /* We have a serious problem when PDF(X) == 0.                      */
  /* Thus our only chance is to assume that the support of the PDF is */
  /* connected. Hence we assume that we have found the cut point.     */
  /*   if ( _unur_iszero(fp) || _unur_iszero(fm) ) { */
  if ( _unur_iszero(fx) ) {
    /* TODO: we should start simple bisection to find the extremal point */
    /* where PDF(X) is non-zero.                                         */
    return 0.;
  } 

  /* check data */
  if( fm-2.*fx+fp < 0.|| fm<1.e-100 || fx<1.e-100 || fp<1.e-100) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,
		  "numerical problems with tail probability");
  }

  /* approximate local concavity + 1 */
  lcplus1 = fp/(fp-fx) + fm/(fm-fx);

  /* approximate derivative */
  df = (fp-fm)/(2.*dx);
  /** TODO: use dPDF if available **/
 
  /* approximate tail probability */
  area = (fx*fx) / (lcplus1 * fabs(df));

  /* check result */
  if (area < 0.) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"tail probability might be negative");
    /* Remark: 
       It is difficult to distinguish between an invalid PDF given by the user,
       a real numerical (round-off) error or simply a point too far from 
       the region of "computational irrelevant".
       Hence we just print a warning and let the calling routine handle the 
       situation.
     */
  }

  if (_unur_isnan(area)) {
    /* When the PDF is extremely flat than NaN might occur. There are two possibilities: 
     * (1) We are far away from the center. Then the tail probabilty is almost 0.
     * (2) PDF(X) is still quite large. Then we should return INFINITY.
     */
    _unur_warning(gen->genid,UNUR_ERR_NAN,"tail probability gives NaN --> assume 0.");
    return 0.;
  }

  /* return area below PDF in tail ("tail probability") */
  return area;

} /* end of _unur_pinv_tailprob() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
