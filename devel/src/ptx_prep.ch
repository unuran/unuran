/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ptx_prep.c                                                  *
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
_unur_ptx_preprocessing (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* 1. Find computational domain and compute PDF area.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  switch (gen->variant) {
  case PTX_VARIANT_PDF:

    /* 1a. Estimate computationally relevant domain (support) of PDF */
    if (_unur_ptx_relevant_support(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    
    /* 1b. Compute area below PDF over relevant domain approximately. */
    if (_unur_ptx_approx_pdfarea(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    
    /* 1c. Compute computational domain where inverse CDF is approximated */
    if (_unur_ptx_computational_domain(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;

    /* 1d. Compute area below PDF with requested accuracy and    */
    /*     store intermediate results from adaptive integration. */
    if (_unur_ptx_pdfarea(gen) != UNUR_SUCCESS) 
      return UNUR_FAILURE;
    break;

  case PTX_VARIANT_CDF:

    /* 1c. Compute computational domain where inverse CDF is approximated */
    if (_unur_ptx_computational_domain_CDF(gen) != UNUR_SUCCESS)
      return UNUR_FAILURE;
    break;

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_ptx_preprocessing() */


/*****************************************************************************/
/** Preprocessing when PDF is given                                         **/
/*****************************************************************************/

int
_unur_ptx_relevant_support ( struct unur_gen *gen )
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
    GEN->bleft = _unur_ptx_searchborder(gen, DISTR.center, GEN->bleft, 
					 &(GEN->dleft), &(GEN->sleft) );
    if (!_unur_isfinite(GEN->bleft)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get left boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }
  
  if(GEN->sright) {
    GEN->bright = _unur_ptx_searchborder(gen, DISTR.center, GEN->bright,
					  &(GEN->dright), &(GEN->sright) );
    if (!_unur_isfinite(GEN->bright)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot get right boundary of relevant domain.");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & PTX_DEBUG_SEARCHBD)
    _unur_ptx_debug_relevant_support(gen);
#endif

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_ptx_relevant_support() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_searchborder (struct unur_gen *gen, double x0, double bound,
			 double *dom, int *search)
     /*----------------------------------------------------------------------*/
     /* [1a.] Find left or right hand border of relevant domain.             */
     /*                                                                      */
     /* Calculate domain of computational relevant region.                   */
     /* Start at 'x0' and search towards 'bound'.                            */
     /* The boundary points of this domain are approximately given as        */
     /*      PDF(x0) * PTX_PDFLLIM                                          */
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
     /*   return UNUR_INFINITY                                               */
     /*----------------------------------------------------------------------*/
{
  double x;         /* current and previous searching point */
  double xs, xl;    /* point where PDF is less than and larger than threshold */
  double fx;        /* PDF at x */
  double fs, fl;    /* PDF at xs and xl */
  double fllim;     /* threshold value */
  double fulim;     /* threshold for detecting discontinuity at boundary */

  /* threshold value where we stop searching */
  fllim = PDF(x0) * PTX_PDFLLIM;
  fulim = 1.e4 * fllim;

  /* we already have checked PDF(center). but who knowns. */
  if (fllim <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(center) too small");
    return UNUR_INFINITY;
  }

  /* starting point */
  xl = x0; 
  fl = UNUR_INFINITY;
  x = _unur_arcmean(x0,bound);

  /* find points where PDF values bracket threshold: */
  /*   fs = PDF(xs) <= fllim <= PDF(xl) = fl           */
  while ( (fx=PDF(x)) > fllim ) {
    if (_unur_FP_same(x,bound))
      return bound;
    xl = x; fl = fx;
    x = _unur_arcmean(x,bound);
  }
  xs = x; fs = fx;

  /* check sign */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0");
    return UNUR_INFINITY;
  }

  /* decrease length of bracket if necessary */
  while (!_unur_FP_same(xs,xl)) {

    /* truncate domain if possible */
    if (_unur_iszero(fs)) {
      *dom = xs;
    }

    /* new point */
    x = xs/2. + xl/2.;
    fx = PDF(x);

    if (fx < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0");
      return UNUR_INFINITY;
    }

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
  
} /* end of _unur_ptx_searchborder() */

/*---------------------------------------------------------------------------*/

int
_unur_ptx_approx_pdfarea (struct unur_gen *gen )
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
     /*   However, then we probably need not method PTX                     */
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
    tol = PTX_UERROR_AREA_APPROX * GEN->area;

    /* check center of distribution */
    DISTR.center = _unur_max(DISTR.center, GEN->bleft);
    DISTR.center = _unur_min(DISTR.center, GEN->bright);

    /* compute area */
    GEN->area  =
      _unur_lobatto_adaptive(_unur_ptx_eval_PDF, gen,
			     GEN->bleft, DISTR.center - GEN->bleft, tol, NULL)
      + _unur_lobatto_adaptive(_unur_ptx_eval_PDF, gen,
			       DISTR.center, GEN->bright - DISTR.center, tol, NULL);

    /** FIXME: uerror **/

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
  /* write info into LOG file */
  if (gen->debug & PTX_DEBUG_SEARCHBD)
    _unur_ptx_debug_pdfarea(gen,TRUE);
#endif

  return res;

} /* end of _unur_ptx_approx_pdfarea() */

/*---------------------------------------------------------------------------*/

int
_unur_ptx_pdfarea (struct unur_gen *gen)
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
     /*   However, then we probably need not method PTX                     */
     /*----------------------------------------------------------------------*/
{
  double tol;   /* tolerated integration error */

  /* the tolerated U-error depends on the given u-resolution.        */
  /* however, we have to make some adjustments for integration error */
  /* in order to achieve the desired precison goal.                  */
  tol = GEN->u_resolution * GEN->area * PTX_UERROR_CORRECTION * PTX_UTOL_CORRECTION;

  /* check center of distribution */
  DISTR.center = _unur_max(DISTR.center, GEN->bleft);
  DISTR.center = _unur_min(DISTR.center, GEN->bright);

  /* create object that contains approximate CDF */
  GEN->aCDF = _unur_lobatto_create(_unur_ptx_eval_PDF, gen,
				   GEN->bleft, DISTR.center, GEN->bright,
				   tol, NULL, PTX_MAX_LOBATTO_IVS);

  /* create object that contains approximate CDF */
  GEN->aTRX = _unur_lobatto_create(_unur_ptx_eval_dTRX, gen,
				   GEN->bleft, DISTR.center, GEN->bright,
				   tol, _unur_ptx_uerror, PTX_MAX_LOBATTO_IVS);

  /* retrieve area below the PDF */
  GEN->area = _unur_lobatto_integral(GEN->aCDF);


  /* check estimated area */
  if ( !_unur_isfinite(GEN->area) || _unur_iszero(GEN->area) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot estimate area below PDF");
    return UNUR_FAILURE;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & PTX_DEBUG_SEARCHBD)
    _unur_ptx_debug_pdfarea(gen,FALSE);
#endif

  return UNUR_SUCCESS;

} /* end of _unur_ptx_pdfarea() */

/*---------------------------------------------------------------------------*/

int
_unur_ptx_computational_domain (struct unur_gen *gen)
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
  /* We use the given U-reslution * PTX_TAILCUTOFF_FACTOR * PDFarea.   */
  tailcut_error = GEN->u_resolution * PTX_TAILCUTOFF_FACTOR;
  tailcut_error = _unur_min( tailcut_error, PTX_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PTX_UERROR_CORRECTION;

  /* length of current working domain */
  range = GEN->bright-GEN->bleft;

  /* compute cut-off points for tails */
  if(GEN->sleft) {
    GEN->bleft = _unur_ptx_cut( gen, GEN->dleft, GEN->bleft, -range, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

  if(GEN->sright) {
    GEN->bright = _unur_ptx_cut( gen, GEN->dright, GEN->bright, range, tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & PTX_DEBUG_SEARCHBD)
    _unur_ptx_debug_computational_domain(gen);
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ptx_computational_domain() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_cut( struct unur_gen *gen, double dom, double w, double dw, double crit )
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
     /*   return UNUR_INFINITY                                               */
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
    if (_unur_iszero(lc)) {
      xnew = x + fx/df * log(crit*fabs(df)/(fx*fx));
    }
    else {
      xnew = x + fx/(lc*df) * ( pow(crit*fabs(df)*(lc+1.)/(fx*fx),lc/(lc+1.)) - 1.);
    }

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

} /* end of _unur_ptx_cut() */


/*****************************************************************************/
/** Preprocessing when CDF is given                                         **/
/*****************************************************************************/

int
_unur_ptx_computational_domain_CDF (struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* 1c. Compute computational domain where inverse CDF is approximated   */
     /*     (interval where we safely can compute coefficients of            */
     /*     interpolating polynomial).                                       */
     /*     Use CDF.                                                         */
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
  double fl, fr;

  /* first we have to chech the domain of the distribution */
  fl = CDF(DISTR.domain[0]);
  fr = CDF(DISTR.domain[1]);
  if (_unur_FP_approx(fl,fr)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain too narrow");
    return UNUR_FAILURE;
  }

  /* parameters for tail cut-off points: maximal area in tails          */
  /* We use the given U-reslution * PTX_TAILCUTOFF_FACTOR * PDFarea.   */
  tailcut_error = GEN->u_resolution * PTX_TAILCUTOFF_FACTOR;
  tailcut_error = _unur_min( tailcut_error, PTX_TAILCUTOFF_MAX );
  tailcut_error = _unur_max( tailcut_error, 2*DBL_EPSILON );
  tailcut_error *= GEN->area * PTX_UERROR_CORRECTION;

  /* length of current working domain */
  range = GEN->bright-GEN->bleft;

  /* compute cut-off points for tails */
  if(GEN->sleft) {
    GEN->bleft = _unur_ptx_cut_CDF( gen, GEN->dleft, DISTR.center, 0.5*tailcut_error, tailcut_error);
    if ( !_unur_isfinite(GEN->bleft) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find left boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

  if(GEN->sright) {
    GEN->bright = _unur_ptx_cut_CDF( gen, GEN->dright, DISTR.center, 1.-tailcut_error, 1.-0.5*tailcut_error);
    if ( !_unur_isfinite(GEN->bright) ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot find right boundary for computational domain");
      return UNUR_FAILURE;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & PTX_DEBUG_SEARCHBD)
    _unur_ptx_debug_computational_domain(gen);
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ptx_computational_domain_CDF() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_cut_CDF( struct unur_gen *gen, double dom, double x0, double ul, double uu )
     /*----------------------------------------------------------------------*/
     /* [1c.] Calculate cut-off points for computational domain of           */
     /* distribution using CDF.                                              */
     /* The area outside the cut-off point is within the interval ul and uu. */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   dom  ... boundary of domain / support of distribution              */
     /*   x0   ... starting point for searching cut-off point                */
     /*   ul   ... lower u-error criterium for tail cut off                  */
     /*   uu   ... upper u-error criterium for tail cut off                  */
     /*                                                                      */
     /* return:                                                              */
     /*   cut-off point                                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return UNUR_INFINITY                                               */
     /*----------------------------------------------------------------------*/
{
  double x;         /* current and previous searching point */
  double xs, xl;    /* point where CDF is less than and larger than threshold */
  double fx;        /* CDF at x */
  double fs, fl;    /* CDF at xs and xl */
  double f0, fdom;  /* CDF at x0 and dom */
  double dx;        /* step size for searching for relevant point */

  /* check length of interval */
  if (_unur_FP_same(x0,dom))
      return x0;

  /* check u-error (protect against u-values too close to 1) */
  if (1.-ul < 4*DBL_EPSILON) ul = 1. - 4*DBL_EPSILON;
  if (1.-uu < 2*DBL_EPSILON) ul = 1. - 2*DBL_EPSILON;

  /* starting point */
  x = x0;
  f0 = CDF(x0);
  fdom = CDF(dom);

  /* check for starting point */
  if (_unur_iszero(f0)) {
    for (dx=0.1; f0<ul; dx*=10.) {
      dom = x0; fdom = f0;
      x0 += dx;
      f0 = CDF(x0);
      if (!_unur_isfinite(x0))
	return UNUR_INFINITY;
    }
  }
  if (_unur_isone(f0)) {
    for (dx=0.1; f0>ul; dx*=10.) {
      dom = x0; fdom = f0;
      x0 -= dx;
      f0 = CDF(x0);
      if (!_unur_isfinite(x0))
	return UNUR_INFINITY;
    }
  }

  /* find points where CDF values bracket threshold critical u values */
  if ( (f0 < ul && fdom < ul) || (f0 > uu && fdom > uu) ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"CDF too small/large on given domain");
    return dom;
  }
  if (f0 >= ul && f0 <= uu) {
    /* we already have reached our goal */
    return x0;
  }
  if ( (x0 < dom && _unur_FP_greater(f0,fdom)) ||
       (x0 > dom && _unur_FP_less(f0,fdom)) ) {
    /* there is something wrong */
    return UNUR_INFINITY;
  }

  /* bracket */
  if (x0 > dom) {
    xs = dom; fs = fdom;
    xl = x0; fl = f0; }
  else {
    xs = x0; fs = f0;
    xl = dom; fl = fdom;
  }    
  x = x0;

  /* decrease length of bracket */
  while (!_unur_FP_same(xs,xl)) {

    /* new point */
    x = _unur_arcmean(xs,xl);
    fx = CDF(x);

    /* check result */
    if (fx >= ul && fx <= uu) {
      return x;
    }

    /* update bracket */
    if (fx < ul) {
      xs = x; fs = fx;
    }
    else {
      xl = x; fl = fx;
    }
  }

  return x;

} /* end of _unur_ptx_cut_CDF() */


/*****************************************************************************/

double
_unur_ptx_Udiff (struct unur_gen *gen, double x, double h, double utol)
     /*----------------------------------------------------------------------*/
     /* Compute difference CDF(x+h)-CDF(x) (approximately), where CDF is     */
     /* the integral of the given (quasi-) density.                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left boundary point of interval                            */
     /*   h   ... length of interval                                         */
     /*   tol ... tolerated ABSOLUTE error                                   */
     /*                                                                      */
     /* return:                                                              */
     /*    (approximate) difference CDF(x+h) - CDF(x)                        */
     /*----------------------------------------------------------------------*/
{
  switch (gen->variant) {
  case PTX_VARIANT_PDF:
/*     return _unur_ptx_Udiff_lobatto(gen,x,h,utol); */

    /** FIXME:PDF **/

  case PTX_VARIANT_CDF:
    return CDF(x+h) - CDF(x);

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_INFINITY;
  }

} /* end of _unur_ptx_Udiff() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_Tdiff (struct unur_gen *gen, double x, double h, double utol)
     /*----------------------------------------------------------------------*/
     /* Compute difference CDF(x+h)-CDF(x) (approximately), where CDF is     */
     /* the integral of the given (quasi-) density.                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left boundary point of interval                            */
     /*   h   ... length of interval                                         */
     /*   tol ... tolerated ABSOLUTE error                                   */
     /*                                                                      */
     /* return:                                                              */
     /*    (approximate) difference CDF(x+h) - CDF(x)                        */
     /*----------------------------------------------------------------------*/
{
  switch (gen->variant) {
  case PTX_VARIANT_PDF:

    /** FIXME:PDF **/

  case PTX_VARIANT_CDF:
    return (INVCDFIN(CDF(x+h)) - INVCDFIN(CDF(x)));

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_INFINITY;
  }

} /* end of _unur_ptx_Udiff() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_uerror (struct unur_gen *gen, double terror, double t)
     /*----------------------------------------------------------------------*/
     /* estimate u-error from t error at value t.                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   terror ... t-error                                                 */
     /*   t      ... point t                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*    (approximate) u-error                                             */
     /*----------------------------------------------------------------------*/
{
  double uerror;

  uerror = terror * PDFIN(t);

  printf("t=%g, terror=%g, PDF(t)=%g, uerror=%g\n",t,terror,PDFIN(t),uerror);

  return uerror;
  
} /* end of _unur_ptx_uerror() */

/*---------------------------------------------------------------------------*/
