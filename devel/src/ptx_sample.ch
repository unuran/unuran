/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ptx_sample.c                                                *
 *                                                                           *
 *   Sampling routines.                                                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008 Wolfgang Hoermann and Josef Leydold                  *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************/
/**  Sampling routines                                                      **/
/*****************************************************************************/

double
_unur_ptx_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return UNUR_INFINITY                                               */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_PTX_GEN,UNUR_INFINITY);

  /** FIXME: error: this function is not available **/

  return UNUR_INFINITY;

} /* end of _unur_ptx_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_ptx_eval_approxinvcdf( const struct unur_gen *gen, double t )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
     /* (internal call)                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   t   ... argument for inverse CDF (0<=u<=1, no validation!)         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return UNUR_INFINITY                                               */
     /*----------------------------------------------------------------------*/
{
  int i;
  double x;

  int bsmin, bsmax;   /* boundaries for binary search */

  /* check arguments */
  CHECK_NULL(gen,UNUR_INFINITY);  COOKIE_CHECK(gen,CK_PTX_GEN,UNUR_INFINITY);

  /* look up in guide table and search for interval */
  /*   i = GEN->guide[(int)(u * GEN->guide_size)]; */
/*   i = 0; */
/*   while (GEN->iv[i+1].cdfi < t) */
/*     i++; */


  /* binary search for interval */
  bsmin = 0;
  bsmax = GEN->n_ivs;
  do {
    i = (bsmin + bsmax) / 2;
    if (GEN->iv[i].cdfi < t) { 
      bsmin = i; 
    }
    else {
      bsmax = i; 
    }
  } while (bsmin<bsmax-1);
  i = bsmin;

  /** FIXME: error handling **/




  /* rescale for range (0, uq(CDF(right))-uq(CDF(left)) for interval */
  t -= GEN->iv[i].cdfi;

  /* evaluate polynomial */
  x = _unur_ptx_newton_eval(t, GEN->iv[i].ti, GEN->iv[i].zi, GEN->order);

  /* return point (add left boundary point to x) */
  return (GEN->iv)[i].xi + x;

} /* end of _unur_ptx_eval_approxinvcdf() */

/*---------------------------------------------------------------------------*/

double
unur_ptx_eval_approxinvcdf( const struct unur_gen *gen, double t )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial interpolation of inverse CDF at u                */
     /* (user call)                                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   t   ... argument for inverse CDF (0<=u<=1)                         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (approximate inverse CDF)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return UNUR_INFINITY                                               */
     /*----------------------------------------------------------------------*/
{
  double x;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_INFINITY );
  if ( gen->method != UNUR_METH_PTX ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return UNUR_INFINITY;
  }
  COOKIE_CHECK(gen,CK_PTX_GEN,UNUR_INFINITY);

  
  if (t<GEN->Tmin || t>GEN->Tmax) {
    /** FIXME: this case is extremely unlikely **/
    printf("junk\n");
  }

  /*   if ( ! (u>0. && u<1.)) { */
  /*     if ( ! (u>=0. && u<=1.)) { */
  /*       _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]"); */
  /*     } */
  /*     if (u<=0.) return DISTR.domain[0]; */
  /*     if (u>=1.) return DISTR.domain[1]; */
  /*     return u;  /\* = NaN *\/ */
  /*   } */
  
  /* compute inverse CDF */
  x = _unur_ptx_eval_approxinvcdf(gen,t);

  /* validate range */
  if (x<DISTR.domain[0]) x = DISTR.domain[0];
  if (x>DISTR.domain[1]) x = DISTR.domain[1];

  return x;

} /* end of unur_ptx_eval_approxinvcdf() */

/*****************************************************************************/

int
unur_ptx_estimate_error( const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE )
     /*----------------------------------------------------------------------*/
     /* Estimate maximal u-error and mean absolute error (MAE) by means of   */
     /* Monte-Carlo simulation.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   samplesize ... sample size for Monte Carlo simulation              */
     /*   max_error  ... pointer to double for storing maximal u-error       */
     /*   MAE        ... pointer to double for storing MA u-error            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double score;

  /* check arguments */
  _unur_check_NULL(GENTYPE, gen, UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_PTX_GEN,UNUR_ERR_COOKIE);

  /* run test */
  score = unur_test_inverror(gen, max_error, MAE, 1.e-20, samplesize, 
			     FALSE, FALSE, FALSE, NULL);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_ptx_estimate_error() */

/*****************************************************************************/
