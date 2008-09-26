/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv_newton.c                                                *
 *                                                                           *
 *   Routines for Newton interpolating polynomial.                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/** Compute coefficients for Newton  polynomial                             **/
/*****************************************************************************/

int
_unur_pinv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Create table for Newton interpolation.                               */
     /* The computational domain must be already computed / given.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double utol;               /* tolerated maximum U-error */
  double maxerror;           /* maximum U-error in a particular interval */ 
  double h;                  /* step size / length of intervals */
  int i;                     /* number of interval at work */
  int iter;                  /* number of iterations */
  int cont;                  /* whether we have to continue or loop */
  int k;                     /* auxiliary variable */

  int right_bd;              /* whether we have reach right boundary */
  double iarea;              /* area below PDF in subinterval */

  double chebyshev[MAX_ORDER+1]; /* Chebyshev points */
  double xval[MAX_ORDER+1];  /* x-values for construction points for Newton polynomial */

  int n_decr_h = 0;          /* number of steps where h is decreased */
  int n_incr_h = 0;          /* number of steps where h is increased */
  int n_incr_h_left = 0;     /* number of steps where h is increased at left boundary 
				to avoid round-off error */
  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_ERR_COOKIE);

  /* compute construction points for Chebyshev points */
  _unur_pinv_chebyshev_points(GEN->order,chebyshev);

  /* tolerated U-error */
  utol = GEN->u_resolution * GEN->area * PINV_UERROR_CORRECTION;

  /* initialize step size for subintervals */
  h = (GEN->bright-GEN->bleft)/128.;

  /* initialize array of interval: starting interval */
  if (_unur_pinv_interval( gen, 0, GEN->bleft, 0.) != UNUR_SUCCESS) 
    return UNUR_ERR_GEN_CONDITION;

  /* initialize counter and control variables */
  i = 0;        /* start at interval 0 ;-) */
  cont = TRUE;  /* there is at least one iteration */

  /* compute intervals and coefficients of Newton polynomials */
  for (iter=0; cont ; iter++) {

    /* check number of iterations */
    if (iter >= PINV_MAX_ITER_IVS) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "maximum number of iterations exceeded");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* right boundary reached ? */
    if(GEN->iv[i].xi+h > GEN->bright) {
      h = GEN->bright - GEN->iv[i].xi;
      cont = FALSE;    /* we probably can stop after this iteration */
      right_bd = TRUE; /* we have reached right boundary */
    }
    else {
      right_bd = FALSE;
    }

    /* construction points for Newton interpolation polynomial: */
    /* uses Chebyshev points for x values.                      */
    for(k=0; k<=GEN->order; k++)
      xval[k] = GEN->iv[i].xi + h * chebyshev[k];

    /* compute Newton interpolation polynomial */
    if (_unur_pinv_newton_create(gen,&(GEN->iv[i]),xval,utol) != UNUR_SUCCESS) {
      /* area below PDF = 0 or serious round-off errors */
      /* we have to deal with three cases: 
       *   1. left boundary
       *   2. interior
       *   3. right boundary
       */
      
      /* -- 1. left boundary (first interval) -- */
      if (i==0) { 
	/* check area: CDF(x) - CDF(left) */
	iarea = _unur_pinv_Udiff(gen, GEN->bleft, h, utol);

	if (iarea < 0.1*utol) {
	  /* cut right boundary */
	  GEN->bleft = GEN->iv[i].xi + h;
	  GEN->iv[i].xi = GEN->bleft;
	}
	else if (n_incr_h_left < 100) {
	  /* increase step size */
	  h *= 2.;
	  ++n_incr_h;
	  ++n_incr_h_left;
	}
	else {
	  /* emergency exit: since we were not successfull til now --> abort */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		      "PDF too close to 0 on relevant part of domain --> abort");
	  return UNUR_ERR_GEN_CONDITION;
	}
	continue;
      }
      
/*       /\* -- 3. right boudary -- *\/ */

/*       /\* check area: CDF(bright) - CDF(x) *\/ */
/*       iarea = _unur_pinv_Udiff(gen, GEN->iv[i].xi, GEN->bright - GEN->iv[i].xi, utol); */

/*       if (iarea < 0.1*utol) { */
/* 	/\* cut right boundary *\/ */
/* 	GEN->bright = GEN->iv[i].xi; */
/* 	cont = FALSE; */
/* 	break; */
/*       } */

/*       /\* -- 2. interior interval -- *\/ */
/*       if (!right_bd) { */
/* 	/\* try to replace interval [xi,xi+h] by [xi,bright]         *\/ */
/* 	/\* (avoid round-off errors by making h a little bit larger) *\/ */
/* 	h = 1.1*(GEN->bright - GEN->iv[i].xi); */
/* 	++n_incr_h; */
/* 	continue; */
/*       } */

/*       /\* -- no idea what to do now --> abort *\/ */
/*       _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION, */
/* 		  "PDF too close to 0 on relevant part of domain --> abort"); */
/*       return UNUR_ERR_GEN_CONDITION; */

      
      /* -- 2. interior interval -- */
      if (!right_bd) {
	/* try to replace interval [xi,xi+h] by [xi,bright]         */
	/* (avoid round-off errors by making h a little bit larger) */

	h = 1.1*(GEN->bright - GEN->iv[i].xi);
	++n_incr_h;
	continue;
      }
      
      /* -- 3. right boudary -- */

      /* check area: CDF(bright) - CDF(x) */
      iarea = _unur_pinv_Udiff(gen, GEN->iv[i].xi, GEN->bright - GEN->iv[i].xi, utol);

      if (iarea < 10*utol) {
	/* cut right boundary */
	if (iarea > utol) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,
			"PDF near right boundary too close to 0 --> cut domain");
	}
	GEN->bright = GEN->iv[i].xi;
	cont = FALSE;
	break;
      }

      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		    "PDF too close to 0 on relevant part of domain --> abort");
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    /* estimate error of Newton interpolation */
    maxerror = _unur_pinv_newton_maxerror(gen,&(GEN->iv[i]),xval,utol);

    if (!(maxerror <= utol)) {
      /* error too large: reduce step size */
      h *= (maxerror > 4.*utol) ? 0.81 : 0.9;
      cont = TRUE;  /* we need another iteration */
      ++n_decr_h;
    }

    else {
      /* create next interval */
      if ( _unur_pinv_interval( gen, i+1, GEN->iv[i].xi+h, 
				GEN->iv[i].cdfi +(GEN->iv)[i].ui[GEN->order-1])
	   /* cdfi holds CDF value at the left border of the interval,                  */
	   /* ui[order-1] holds area below PDF in interval, i.e. CDF(right) - CDF(left) */
      	   != UNUR_SUCCESS )
	return UNUR_ERR_GEN_CONDITION;

      /* increase step size for very small errors */
      if (maxerror < 0.3*utol) {
	h *= (maxerror < 0.1*utol) ? 2. : 1.2;
	++n_incr_h;
      }
      
      /* continue with next interval */
      i++;
    }
  }

  /* update size of array */
  GEN->iv = _unur_xrealloc( GEN->iv, (GEN->n_ivs+1) * sizeof(struct unur_pinv_interval) );
  
  /* set range for uniform random numbers */
  /* Umin = 0, Umax depends on area below PDF, tail cut-off points and round-off errors */
  GEN->Umax = GEN->iv[GEN->n_ivs].cdfi;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & PINV_DEBUG_SEARCHBD)
    _unur_pinv_debug_create_table(gen,iter,n_incr_h,n_decr_h);
#endif

  /* o.k. */
  return UNUR_SUCCESS;
}  /* end of _unur_pinv_create_table() */

/*---------------------------------------------------------------------------*/

int 
_unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx )
     /*----------------------------------------------------------------------*/
     /* make a new interval i with left boundary point x and CDF(x).         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   i    ... index (number) of interval                                */
     /*   x    ... left boundary point of new interval                       */
     /*   cdfx ... CDF at x                                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
     /* We use Newtons interpolating polynomial. The construction points are */
     /* created by using Chebyshev points for the x-coordinates.             */
     /*                                                                      */
     /* Alternative implementation would be:                                 */
     /* use Chebyshev points for u-coordinates.                              */
     /* advantage:                                                           */
     /*  - smaller tables as u-coordinates need not be stored.               */
     /*  - fewer intervals required (especially in the tails).               */
     /*  - routine _unur_pinv_newton_testpoints() not required               */
     /*    (test points are fixed Chebyshev points -- closed form).          */ 
     /* disadvantage:                                                        */
     /*  - numerically less stable and more expensive                        */
     /*    as it requires explicit inversion.                                */
     /*----------------------------------------------------------------------*/
{
  struct unur_pinv_interval *iv;
#ifdef PINV_USE_CDFTABLE
  struct unur_pinv_CDFtable *CDFtable = GEN->CDFtable;
#endif

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);

  /* check for free intervals */
  if (i >= GEN->max_ivs) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"maximum number of intervals exceeded");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* set values */
  iv = GEN->iv+i;     /* pointer to interval */
  iv->xi = x;         /* left boundary of interval */
  iv->cdfi = cdfx;    /* CDF at left boundary */
  COOKIE_SET(iv,CK_PINV_IV);

  /* allocate space for coefficients for Newton interpolation */
  iv->ui = _unur_xmalloc( GEN->order * sizeof(double) );
  iv->zi = _unur_xmalloc( GEN->order * sizeof(double) );

  /* update size of array (number of intervals) */
  GEN->n_ivs = i;


#ifdef PINV_USE_CDFTABLE
  /* set bookmark in table of CDF values */
  if (CDFtable != NULL) {
    /* search for first entry in interval. */
    /* we can continue from the position in the last interval. */
    /* otherwise, to restart from the first entry uncomment this line */
    /*    CDFtable->cur_iv = 0; */
    
    while (CDFtable->cur_iv < CDFtable->n_values &&
	   CDFtable->values[CDFtable->cur_iv].x < x) 
      ++(CDFtable->cur_iv);
  }
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_pinv_interval() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
			  double *xval, double utol)
     /*----------------------------------------------------------------------*/
     /* 2a. Compute coefficients for Newton interpolation within a           */
     /* subinterval of the domain of the distribution.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   h    ... length of interval                                        */
     /*   xval ... x-values for constructing polynomial                      */
     /*   utol ... maximal tolerated u-error                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double *ui = iv->ui;   /* u-values for Newton interpolation */
  double *zi = iv->zi;   /* coefficients of Newton interpolation */
  double xi, dxi;        /* boundary and length of i-th subinterval */
  double area;           /* integral of PDF over subinterval */
  int i,k;               /* auxiliary variables */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* We use a smaller tolerance for the Gauss-Lobatto integral */
  utol *= PINV_UTOL_CORRECTION;

  /* compute tuples (ui,zi) for constructing polynomials */
  for(i=0; i<GEN->order; i++) {

    /* left boundary and length of subinterval for integration */
    xi = xval[i];
    dxi = xval[i+1]-xval[i];

    /* compute integral of PDF in interval (xi,xi+dxi) */
    area = _unur_pinv_Udiff(gen, xi, dxi, utol);
    if (_unur_iszero(area)) return UNUR_ERR_SILENT;

    /* construction points of interpolation polynomial of CDF^{-1} */
    ui[i] = (i>0) ? (ui[i-1]+area) : area;
    /* rescaled corresponding values of CDF^{-1} */ 
    zi[i] = dxi/area;
  }
  /* Remark: ui[GEN->order-1] is the probability of the interval */

  /* compute coefficients of interpolation polynomial */
  for(k=1; k<GEN->order; k++) {
    for(i=GEN->order-1; i>k; i--) {
      zi[i] = (zi[i]-zi[i-1]) / (ui[i]-ui[i-(k+1)]);
    }
    zi[k] = (zi[k]-zi[k-1]) / ui[k];
  }

  /* check result */
  for (i=0; i<GEN->order; i++) {
    if (!_unur_isfinite(zi[i])) 
      return UNUR_ERR_SILENT;
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_create() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_chebyshev_points (int order, double *pt)
     /*----------------------------------------------------------------------*/
     /* [2a.] Compute Chebyshev points.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   order ... order of polynomial                                      */
     /*   pt    ... pointer to array of size (order+1) for storing points    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;
  double phi = M_PI*0.5/(order+1); /* constant for computing Chebyshev points */
  
  pt[0] = 0.;
  for(i=1; i<order; i++)
    pt[i] = sin(i*phi) * sin((i+1)*phi)/cos(phi);
  pt[order] = 1.;

  return UNUR_SUCCESS;
} /* end of _unur_pinv_chebyshev_points() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_eval ( double q, double ui[], double zi[], int order )
     /*----------------------------------------------------------------------*/
     /* 2b. evaluate Newton interpolation polynomial using Horner scheme.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   q     ... argument                                                 */
     /*   ui    ... coefficients of polynomial (increasing order)            */
     /*   zi    ... coefficients of polynomial (increasing order)            */
     /*   order ... order of polynomial                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   value of interpolating polynomial at u=q                           */
     /*----------------------------------------------------------------------*/
{
  int k;
  double chi;

  chi = zi[order-1];
  for (k=order-2; k>=0; k--)
    chi = chi*(q-ui[k])+zi[k];

  return (chi*q);
} /* end of _unur_pinv_newton_eval() */

/*---------------------------------------------------------------------------*/

double
_unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv,
			    double xval[], double utol)
     /*----------------------------------------------------------------------*/
     /* 2c. Estimate maximal error of Newton interpolation in subinterval.   */
     /*     In addition it makes a simple check for monotonicity of the      */
     /*     inverse CDF.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to current interval                               */
     /*   xval ... x-values for constructing polynomial                      */
     /*   utol ... maximal tolerated u-error                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   estimated maximal u-error, or                                      */
     /*   1000 whenever the inverse CDF is not monotone                      */
     /*----------------------------------------------------------------------*/
{
  double x0 = iv->xi;    /* left boundary point of interval */
  double *ui = iv->ui;   /* u-values for Newton interpolation  */
  double *zi = iv->zi;   /* coefficient of Newton interpolation  */

  double maxerror = 0.;  /* maximum error */
  double uerror;         /* error for given U value */
  double x;              /* x = CDF^{-1}(U) */
  double u;              /* u = CDF(x) */

  double testu[MAX_ORDER]; /* array of U values for testing */ 
  int i;                 /* aux variable */

  /* check arguments */
  COOKIE_CHECK(gen,CK_PINV_GEN,UNUR_FAILURE);
  COOKIE_CHECK(iv,CK_PINV_IV,UNUR_FAILURE);

  /* We use a smaller tolerance for the Gauss-Lobatto integral */
  utol *= PINV_UTOL_CORRECTION;

  /* get U values for test (points with maximal worst case error) */
  _unur_pinv_newton_testpoints(GEN->order,ui,testu);

  /* calculate the max u-error at the test points */
  for(i=0; i<GEN->order; i++){
    
    /* inverse CDF for U test point */
    x = _unur_pinv_newton_eval(testu[i], ui, zi, GEN->order);

    /* check for monotonicity */
    if (! (xval[i] <= x0+x && x0+x <=xval[i+1]) ) {
      /* not monotone */
      return 1000.;
    }

    /* estimate CDF for interpolated x value */
    if (i==0 || xval==NULL)
      u = _unur_pinv_Udiff(gen, x0, x, utol);
    else
      u = ui[i-1] + _unur_pinv_Udiff(gen, xval[i], x+x0-xval[i], utol);

    /* check u-value */
    if (!_unur_isfinite(u))
      return INFINITY;

    /* compute u-error */
    uerror = fabs(u - testu[i]);

    /* update maximal error */
    if (uerror>maxerror) maxerror = uerror;
  }

  return maxerror;
} /* end of _unur_pinv_newton_maxerror() */

/*---------------------------------------------------------------------------*/

int
_unur_pinv_newton_testpoints (int order, double ui[], double utest[])
     /*----------------------------------------------------------------------*/
     /* [2c.] calculates the local maxima of the polynomial.                 */
     /* used as control points for error estimate.                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*    order ... order of interpolation polynomial                       */
     /*    ui    ... u-values of interpolation                               */ 
     /*    utest ... pointer to array for storing control points             */
     /*                                                                      */
     /* return:                                                              */
     /*    u-values of control points in the array utest                     */
     /*----------------------------------------------------------------------*/
{
  int k,j,i;
  double sum, qsum,x;
  
  /* compute approximate maxima of error polynomial */
  for(k=0; k<order; k++) {

    /* first approximation: use mean of consecuting construction points */
    x = (k>0) ? 0.5*(ui[k-1]+ui[k]) : 0.5*ui[k];

    /* make two iterations for root finding */
    for(j=1; j<=2; j++) {
      sum = 1./x;
      qsum = sum*sum;
      for(i=0; i<order; i++){
	sum += 1./(x-ui[i]);
	qsum += 1./((x-ui[i])*(x-ui[i]));
      }
      x += sum/qsum;
    }

    /* store result in utest */
    utest[k] = x;
  }
  
  return UNUR_SUCCESS;
} /* end of _unur_pinv_newton_testpoints() */

/*---------------------------------------------------------------------------*/
