/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_sample.c                                                 *
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

/*****************************************************************************/
/**  Sampling routines                                                      **/
/*****************************************************************************/

double
_unur_tdr_gw_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (original variant by Gilks & Wild)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*   squeeze(x) = f(x0) * exp(sq * (x-x0))                              */
     /*                                                                      */
     /*   left hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                   */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = f(x1) * exp( (Tf)'(x1) *  (x-x1) )                  */
     /*   generation:                                                        */
     /*      X = x1 + 1/(Tf)'(x1) * \log( (Tf)'(x1)/f(x1) * U + 1 )          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   squeeze(x) = 1 / (Tf(x0) + sq * (x-x0))^2                          */
     /*                                                                      */
     /*   left hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                  */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = 1 / (Tf(x1) + (Tf)'(x1) * (x-x1))^2                 */
     /*   generation:                                                        */
     /*      X = x1 + (Tf(x1)^2 * U) / (1 - Tf(x1) * (Tf)'(x1) * U)          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv, *pt;
  double U, V, X;
  double sqx, hx, fx;
  double Tsqx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN.Umin + _unur_call_urng(urng) * (GEN.Umax - GEN.Umin);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* left or right side of hat */
    if (-U < iv->Ahatr) { /* right */
      pt = iv->next;
      /* u unchanged */
    }
    else {                /* left */
      pt = iv;
      U += iv->Ahat;
    }

    /* we have three different types of transformations */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      /* random variate */
      if (pt->dTfx == 0.)
	X = pt->x + U / pt->fx;
      else
	{
	  double t = pt->dTfx * U / pt->fx;
	  if (fabs(t) > 1.e-6)
	    X = pt->x + log(t + 1.) * U / (pt->fx * t);
	  /* x = pt->x + log(t + 1.) / pt->dTfx; is cheaper but numerical unstable */
	  else if (fabs(t) > 1.e-8)
	    /* use Taylor series */
	    X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
	  else
	    X = pt->x + U / pt->fx * (1 - t/2.);
	}

      /* accept or reject */
      hx = pt->fx * exp(pt->dTfx*(X - pt->x));    /* value of hat at x */   
      V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at x */
      
      /* below mininum of density in interval ? */
      if (V <= iv->fx && V <= iv->next->fx)
	return X;

      /* below squeeze ? */
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;     /* value of squeeze at x */
      if (V <= sqx)
	return X;

      break;

    case TDR_VAR_T_SQRT:
      /* random variate */
      if (pt->dTfx == 0.)
	X = pt->x + U /pt->fx;
      else {
	/* it would be less expensive to use:
	   X = pt->x + pt->Tfx/pt->dTfx * (1. - 1./(1. + pt->dTfx * pt->Tfx * U) )
	   however, this is unstable for small pt->dTfx */
	X = pt->x + (pt->Tfx*pt->Tfx*U) / (1.-pt->Tfx*pt->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }

      /* accept or reject */
      Thx = pt->Tfx + pt->dTfx * (X - pt->x);     /* transformed hat at x */ 
      hx = 1./(Thx*Thx);
      V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at x */

      /* below mininum of density in interval ? */
      if (V <= iv->fx && V <= iv->next->fx)
	return X;

      /* below squeeze ? */
      Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (X - iv->x)) : -INFINITY; /* transformed squeeze at x */ 
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      if (V <= sqx)
	return X;
      
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;

      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* value of PDF at x */
    fx = PDF(X);

    /* being above squeeze is bad. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_gw_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    if (V <= fx)
      /* between PDF and squeeze */
      return X;

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }
} /* end of _unur_tdr_gw_sample() */

/*****************************************************************************/

double
_unur_tdr_gw_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv, *pt;
  double U, V, X;
  double fx, sqx, hx;
  double Tfx, Tsqx, Thx;
  int error = 0;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN.Umin + _unur_call_urng(urng) * (GEN.Umax - GEN.Umin);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* left or right side of hat */
    if (-U < iv->Ahatr) { /* right */
      pt = iv->next;
      /* U unchanged */
    }
    else {                /* left */
      pt = iv;
      U += iv->Ahat;
    }

    /* random variate */
    if (pt->dTfx == 0.)
      X = pt->x + U/pt->fx;
    else
      switch( gen->variant & TDR_VARMASK_T ) {
      case TDR_VAR_T_LOG:
	{
	  double t = pt->dTfx * U / pt->fx;
	  if (fabs(t) > 1.e-6)
	    X = pt->x + log(t + 1.) * U / (pt->fx * t);
	  else if (fabs(t) > 1.e-8)
	    /* use Taylor series */
	    X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
	  else
	    X = pt->x + U / pt->fx * (1 - t/2.);
	}
	break;
      case TDR_VAR_T_SQRT:
	X = pt->x + (pt->Tfx*pt->Tfx*U) / (1.-pt->Tfx*pt->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 (in theory!) */
	break;
      case TDR_VAR_T_POW:
	/** TODO **/
	X = 1.;
	break;
      default:  /* this should not happen */
	X = 0.;
      }

    fx = PDF(X);                                /* value of PDF at X */
    Thx = pt->Tfx + pt->dTfx * (X - pt->x);     /* transformed hat at X */ 
    Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (X - iv->x)) : -INFINITY; /* transformed squeeze at X */ 

    switch( gen->variant & TDR_VARMASK_T ) {
    case TDR_VAR_T_LOG:
      Tfx = (fx>0.) ? log(fx) : -INFINITY;
      hx = pt->fx * exp(pt->dTfx*(X - pt->x));    /* value of hat at X */   
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;     /* value of squeeze at X */
      break;
    case TDR_VAR_T_SQRT:
      Tfx = (fx>0.) ? -1./sqrt(fx) : -INFINITY;
      hx = 1./(Thx*Thx);
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
      Tfx = 0.;
      hx = 0.;
      sqx = 0.;
      break;
    default:  /* this should not happen */
      Tfx = 0.; hx = 0.; sqx = 0.;
    }

    /* check result */
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(Tfx,Thx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
      error = 1;
    }
    if (_unur_FP_less(Tfx,Tsqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_gw_debug_sample( gen, iv, pt, X, fx, hx, sqx ); 
#endif

    /* accept or reject */
    V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at X */

    /* below mininum of density in interval ? */
    if (V <= iv->fx && V <= iv->next->fx)
      return X;

    /* below squeeze ? */
    if (V <= sqx)
      return X;

    /* being above squeeze is bad. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_gw_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    if (V <= fx)
      /* between PDF and squeeze */
      return X;

    /* reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }
} /* end of _unur_tdr_gw_sample_check() */

/*****************************************************************************/
/*****************************************************************************/

double
_unur_tdr_ps_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (proportional squeeze)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*                                                                      */
     /*   hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                        */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                       */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  double U, V, X;
  double fx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN.Umin + _unur_call_urng(urng) * (GEN.Umax - GEN.Umin);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum - iv->Ahatr;    /* result: U in (-A_hatl, A_hatr) */

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	/* it would be less expensive to use:
	   X = iv->x + iv->Tfx/iv->dTfx * (1. - 1./(1. + iv->dTfx * iv->Tfx * U) )
	   however, this is unstable for small iv->dTfx */
	X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;

      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* accept or reject */
    V = _unur_call_urng(urng);

    /* squeeze rejection */
    if (V <= iv->sq)
      	return X;

    /* evaluate hat at X:
       get uniform random number between 0 and hat(X) */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      V *= iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      V *= 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* main rejection */
    if (V <= fx)
      return X;

    /* evaluation of PDF is expensive. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_ps_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }

} /* end of _unur_tdr_ps_sample() */

/*****************************************************************************/

double
_unur_tdr_ps_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results (proportional squeeze)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  double U, V, X;
  double fx, hx, Thx, sqx;
  int squeeze_rejection = FALSE;
  int error = 0;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN.Umin + _unur_call_urng(urng) * (GEN.Umax - GEN.Umin);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum - iv->Ahatr;    /* result: U in (-A_hatl, A_hatr) */

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	/* it would be less expensive to use:
	   X = iv->x + iv->Tfx/iv->dTfx * (1. - 1./(1. + iv->dTfx * iv->Tfx * U) )
	   however, this is unstable for small iv->dTfx */
	X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;

      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* accept or reject */
    V = _unur_call_urng(urng);

    /* squeeze rejection */
    if (V <= iv->sq)
      squeeze_rejection = TRUE;

    /* evaluate hat at X:
       get uniform random number between 0 and hat(X) */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x));
      V *= iv->fx * exp(iv->dTfx*(X - iv->x));
      break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx);
      V *= 1./(Thx*Thx);
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
    default:
      return 0.;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* evaluate squeeze */
    sqx = iv->sq*hx;

    /* check result */
    if (_unur_FP_less(X, DISTR.BD_LEFT) || _unur_FP_greater(X, DISTR.BD_RIGHT) ) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(fx, hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
      error = 1;
    }
    if (_unur_FP_less(fx, sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_ps_debug_sample( gen, iv, X, fx, hx, sqx ); 
#endif

    /* squeeze rejection */
    if (squeeze_rejection)
      return X;

    /* main rejection */
    if (V <= fx)
      return X;

    /* evaluation of PDF is expensive. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_ps_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }
} /* end of _unur_tdr_ps_sample_check() */

/*****************************************************************************/
/*****************************************************************************/

double
_unur_tdr_ia_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (immediate acceptance)                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*                                                                      */
     /*   hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                        */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                       */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat,0) */

    /* check for region of immediate acceptance */
    if (U >= - iv->sq * iv->Ahat) {
      /* region of immediate acceptance */
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      /* rejection from region between hat and squeeze */
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    /* result: U in (-A_hat,0) */

    /* U in (-A_hatl, A_hatr) */
    U += iv->Ahatr;

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; /* avoid one multiplication */
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* immedate acceptance */
    if (use_ia)
      return X;

    /* evaluate hat at X */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* from now on we use the auxilliary generator
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

    /* rejection from region between hat and (proportional) squeeze */
    V = _unur_call_urng(urng);

    /* get uniform random number between squeeze(X) and hat(X) */
    V = (iv->sq + (1 - iv->sq) * V) * hx;

    /* evaluate PDF at X */
    fx = PDF(X);

    /* main rejection */
    if (V <= fx)
      return X;


    /* evaluation of PDF is expensive. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_ps_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    /* else reject and try again */
  }

} /* end of _unur_tdr_ia_sample() */

/*****************************************************************************/

double
_unur_tdr_ia_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results (immediate acceptance)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  int use_ia;
  double U, V, X;
  double fx, hx, Thx, sqx;
  int error = 0;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_TDR_GEN,0.);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    U = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    iv =  GEN.guide[(int) (U * GEN.guide_size)];
    U *= GEN.Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat,0) */

    /* check for region of immediate acceptance */
    if (U >= - iv->sq * iv->Ahat) {
      /* region of immediate acceptance */
      U /= iv->sq;
      use_ia = 1;
    }
    else {
      /* rejection from region between hat and squeeze */
      U = (U + iv->sq * iv->Ahat) / (1. - iv->sq);
      use_ia = 0;
    }
    /* result: U in (-A_hat,0) */

    /* U in (-A_hatl, A_hatr) */
    U += iv->Ahatr;

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (iv->dTfx == 0.)
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (iv->dTfx == 0.)
	X = iv->x + U /iv->fx;
      else {
	U *= iv->Tfx; /* avoid one multiplication */
	X = iv->x + (iv->Tfx * U) / (1. - iv->dTfx * U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/
      return 1.;
      break;

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 1.;

    } /* end switch */

    /* evaluate hat at X */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      hx = iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      hx = 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
    default:
      /** TODO **/
      return 1.;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* evaluate squeeze */
    sqx = iv->sq*hx;

    /* check result */
    if (_unur_FP_less(X, DISTR.BD_LEFT) || _unur_FP_greater(X, DISTR.BD_RIGHT) ) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(fx, hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
      error = 1;
    }
    if (_unur_FP_less(fx, sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_ps_debug_sample( gen, iv, X, fx, hx, sqx ); 
#endif

    /* immedate acceptance */
    if (use_ia)
      return X;

    /* from now on we use the auxilliary generator
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

    /* rejection from region between hat and (proportional) squeeze */
    V = _unur_call_urng(urng);

    /* get uniform random number between squeeze(X) and hat(X) */
    V = (iv->sq + (1 - iv->sq) * V) * hx;

    /* main rejection */
    if (V <= fx)
      return X;

    /* evaluation of PDF is expensive. improve the situation! */
    if (GEN.n_ivs < GEN.max_ivs) {
      if (GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
	if ( !_unur_tdr_ps_interval_split(gen, iv, X, fx) ) {
	  /* condition for PDF is violated! */
	  _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	  if (gen->variant & TDR_VARFLAG_PEDANTIC) {
	    /* replace sampling routine by dummy routine that just returns INFINITY */
	    SAMPLE = _unur_sample_cont_error;
	    return INFINITY;
	  }
	}
      }
      else {
	/* no more construction points (avoid to many second if statement above */
	GEN.max_ivs = GEN.n_ivs;
      }
    }

    /* else reject and try again */
  }

} /* end of _unur_tdr_ia_sample_check() */

/*****************************************************************************/
/*****************************************************************************/
