/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_poisson_gen.c                                              *
 *                                                                           *
 *   Special generators for Poisson distribution                             *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void poisson_pdtabl_init( struct unur_gen *gen );
inline static void poisson_pdac_init( struct unur_gen *gen );
inline static void poisson_pprsc_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.dstd        /* data for parameter object         */
#define GEN       gen->data.dstd        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

#define MAX_gen_params  100    /* maximal number of parameters for generator */
#define MAX_gen_iparams  100   /* maximal number of integer param. for gen.  */

/* parameters */
#define theta  (DISTR.params[0])    /* shape */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_poisson_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Poisson distribution            */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);  COOKIE_CHECK(par,CK_DSTD_PAR,0);

  switch (par->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Tabulated Inversion combined with Acceptance Complement */
    if (gen==NULL) return 1; /* test existence only  */
    if (theta < 10.) {
      /* CASE B: Tabulated Inversion */
      _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_poisson_pdtabl );
      poisson_pdtabl_init( gen );
    }
    else { /* theta >= 10. */
      /* CASE A: acceptance complement */
      _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_poisson_pdac );
      poisson_pdac_init( gen );
    }
    return 1;

  case 2:  /* Tabulated Inversion combined with Patchwork Rejection */
    if (gen==NULL) return 1; /* test existence only  */
    if (theta < 10.) {
      /* CASE: Tabulated Inversion --> same as case 1 !! */
      _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_poisson_pdtabl );
      poisson_pdtabl_init( gen );
    }
    else { /* theta >= 10. */
      /* CASE: Patchwork Rejection */
      _unur_dstd_set_sampling_routine( par,gen,unur_stdgen_sample_poisson_pprsc );
      poisson_pprsc_init( gen );
    }
    return 1;

  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default: /* no such generator */
    if (gen) _unur_warning(par->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  
} /* end of _unur_stdgen_poisson_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Poisson Distribution: Tabulated Inversion combined with                   *
 *                       Acceptance Complement                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Poisson distribution with  *
 *               parameter theta > 0.                                        *
 *               Tabulated Inversion for  theta < 10                         *
 *               Acceptance Complement for theta >= 10.                      *
 *                                                                           *
 * REFERENCE: - J.H. Ahrens, U. Dieter (1982): Computer generation of        * 
 *              Poisson deviates from modified normal distributions,         *
 *              ACM Trans. Math. Software 8, 163-179.                        *
 *                                                                           *
 * Implemented by R. Kremer, August 1990                                     *
 * Revised by E. Stadlober, April 1992                                       *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define m    (GEN.gen_iparam[0])
#define ll   (GEN.gen_iparam[1])

#define p0   (GEN.gen_param[0])
#define q    (GEN.gen_param[1])
#define p    (GEN.gen_param[2])
#define pp   ((GEN.gen_param)+3)  /* array of length 36 */
/*---------------------------------------------------------------------------*/

inline static void
poisson_pdtabl_init( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  /* check arguments */
  CHECK_NULL(gen,/*void*/); COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
    GEN.n_gen_iparam = MAX_gen_iparams;
    GEN.gen_iparam = _unur_malloc(GEN.n_gen_param * sizeof(int));
  }

  /* -X- setup code -X- */
  m = (theta > 1.) ? ((int) theta) : 1;
  ll = 0;
  p0 = q = p = exp(-theta);
  /* -X- end of setup code -X- */

} /* end of poisson_pdtabl_init() */

int
unur_stdgen_sample_poisson_pdtabl( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  /* -X- generator code -X- */
  double U;
  int K,i;


  /* check arguments */
  CHECK_NULL(gen,0.); COOKIE_CHECK(gen,CK_DSTD_GEN,0.);
  
  while (1) {
    U = uniform();              /* Step U. Uniform sample */
    K = 0;
    if (U <= p0) 
      return K;

    /* Step T. Table comparison */
    if (ll != 0) {               
      i = (U > 0.458) ? min(ll,m) : 1;
      for (K = i; K <=ll; K++)
	if (U <= pp[K])
	  return K;
      if (ll == 35) continue;
    }

    /* Step C. Creation of new prob. */
    for (K = ll +1; K <= 35; K++) {
      p *= theta / (double)K;
      q += p;
      pp[K] = q;
      if (U <= q) {
	ll = K;
	return K;
      }
    }
    ll = 35;
  }
  
  /* -X- end of generator code -X- */
  
} /* end of unur_stdgen_sample_poisson_pdtabl() */

/*---------------------------------------------------------------------------*/
#undef m 
#undef ll
#undef p0
#undef q 
#undef p 
#undef pp
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#define l     (GEN.gen_iparam[0])

#define s     (GEN.gen_param[0])
#define d     (GEN.gen_param[1])
#define omega (GEN.gen_param[2])
#define b1    (GEN.gen_param[3])
#define b2    (GEN.gen_param[4])
#define c     (GEN.gen_param[5])
#define c0    (GEN.gen_param[6])
#define c1    (GEN.gen_param[7])
#define c2    (GEN.gen_param[8])
#define c3    (GEN.gen_param[9])

#define NORMAL  GEN.gen_aux    /* pointer to normal variate generator        */
/*---------------------------------------------------------------------------*/

inline static void
poisson_pdac_init( struct unur_gen *gen )
     /* Theta >= 10: acceptance complement */
{
  /* check arguments */
  CHECK_NULL(gen,/*void*/); COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
    GEN.n_gen_iparam = MAX_gen_iparams;
    GEN.gen_iparam = _unur_malloc(GEN.n_gen_param * sizeof(int));
  }

  /* -X- setup code -X- */

  /* make a normal variate generator (use default special generator) */
  NORMAL = unur_init( unur_cstd_new( unur_distr_normal(NULL,0) ));
  /* need same uniform random number generator as slash generator */
  NORMAL->urng = gen->urng;

  s = sqrt(theta);
  d = 6. * theta * theta;
  l = (int)(theta - 1.1484);

  /* Step P. Preparations for steps Q and H */
  omega = 0.3989423 / s;
  b1 = 0.416666666667e-1 / theta;
  b2 = 0.3 * b1 * b1;
  c3 = 0.1428571 * b1 * b2;
  c2 = b2 - 15.0 * c3;
  c1 = b1 - 6.0 * b2 + 45.0 * c3;
  c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
  c = 0.1069 / theta;

  /* -X- end of setup code -X- */

} /* end of poisson_pdac_init() */

/*---------------------------------------------------------------------------*/
#define  a0  -0.5000000002
#define  a1   0.3333333343
#define  a2  -0.2499998565
#define  a3   0.1999997049
#define  a4  -0.1666848753
#define  a5   0.1428833286
#define  a6  -0.1241963125
#define  a7   0.1101687109
#define  a8  -0.1142650302
#define  a9   0.1055093006
/*---------------------------------------------------------------------------*/

int
unur_stdgen_sample_poisson_pdac( struct unur_gen *gen )
     /* Theta >= 10: acceptance complement */
{
  /* -X- generator code -X- */
  /* factorial for 0 <= k <= 9 */
  const static int fac[] = {1,1,2,6,24,120,720,5040,40320,362880};

  double t,g,theta_k;
  double gx,gy,px,py,x,xx,delta,v;
  int sign;

  double E, U;
  int K;

  /* check arguments */
  CHECK_NULL(gen,0.); COOKIE_CHECK(gen,CK_DSTD_GEN,0.);


  /* Step N. Normal sample */
  t = unur_sample_cont(NORMAL);
  g = theta + s * t;

  if (g >= 0.) {
    K = (int) g;
    /* Step I. Immediate acceptance */
    if (K >= l) 
      return K;
    /* Step S. Squeeze acceptance */
    U = uniform();
    theta_k = theta - K;
    if (d * U >= theta_k * theta_k * theta_k)
      return K;

    /* FUNCTION F */
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else {  /* k >= 10 */
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8 * delta*delta*delta * (1.-1./(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    /* end FUNCTION F */

    /* Step Q. Quotient acceptance */
    if (gy * (1.0 - U)  <= py * exp(px - gx))
      return K;
  }

  /* Step E. Double exponential sample */
  while (1) {
    do {
      E = - log(uniform());
      U = uniform();
      U = U + U - 1.;
      sign = (U < 0.) ? -1 : 1;
      t = 1.8 + E * sign;
    } while (t <= -0.6744);
    K = (int)(theta + s * t);
    theta_k = theta - K;

    /* FUNCTION F */
    if (K < 10) {
      px = -theta;
      py = exp(K * log(theta)) / fac[K];
    }
    else { /* k >= 10 */
      delta = 0.83333333333e-1 / (double)K;
      delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*K*K));
      v = (theta_k) / (double)K;
      if (fabs(v) > 0.25)
	px = K * log(1. + v) - theta_k - delta;
      else {
	px = K * v * v;
	px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
		  a4)*v+a3)*v+a2)*v+a1)*v+a0;
	px -= delta;
      }
      py = 0.3989422804 / sqrt((double)K);
    }
    x = (0.5 - theta_k) / s;
    xx = x * x;
    gx = -0.5 * xx;
    gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
    /* end FUNCTION F */

    /* Step H. Hat acceptance */
    if (c * sign * U <= py * exp(px + E) - gy * exp(gx + E)) 
      return K;
  }

  /* -X- end of generator code -X- */
  
} /* end of unur_stdgen_sample_poisson_pdac() */

/*---------------------------------------------------------------------------*/
#undef  a0
#undef  a1
#undef  a2
#undef  a3
#undef  a4
#undef  a5
#undef  a6
#undef  a7
#undef  a8
#undef  a9

#undef l
#undef s
#undef d
#undef omega
#undef b1
#undef b2
#undef c 
#undef c0
#undef c1
#undef c2
#undef c3

#undef NORMAL
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Poisson Distribution: Tabulated Inversion combined with                   *
 *                       Patchwork Rejection                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:  - samples a random number from the Poisson distribution with   *
 *              parameter theta > 0.                                         *
 *              Tabulated Inversion for  theta < 10                          *
 *              Patchwork Rejection for theta >= 10.                         *
 *                                                                           *
 * REFERENCE: - H. Zechner (1994): Efficient sampling from continuous and    *
 *              discrete unimodal distributions,                             *
 *              Pd.D. Thesis, 156 pp., Technical University Graz, Austria.   *
 *                                                                           *
 * Implemented by H. Zechner, January 1994                                   *
 * Revised by F. Niederl, July 1994                                          *
 *****************************************************************************
 *                                                                           *
 * Patchwork Rejection:                                                      *
 * The area below the histogram function f(x) is rearranged in its body by   *
 * certain point reflections. Within a large center interval variates are    *
 * sampled efficiently by rejection from uniform hats. Rectangular immediate *
 * acceptance regions speed up the generation. The remaining tails are       *
 * covered by exponential functions.                                         *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#if 0
#define m    (GEN.gen_iparam[0])
#define ll   (GEN.gen_iparam[1])

#define p0   (GEN.gen_param[0])
#define q    (GEN.gen_param[1])
#define p    (GEN.gen_param[2])
#define pp   ((GEN.gen_param)+3)  /* array of length 36 */
#endif
/*---------------------------------------------------------------------------*/

inline static void
poisson_pprsc_init( struct unur_gen *gen )
     /* theta < 10: Tabulated inversion */
{
  /* check arguments */
  CHECK_NULL(gen,/*void*/); COOKIE_CHECK(gen,CK_DSTD_GEN,/*void*/);

  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = MAX_gen_params;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
    GEN.n_gen_iparam = MAX_gen_iparams;
    GEN.gen_iparam = _unur_malloc(GEN.n_gen_param * sizeof(int));
  }

#if 0
  /* -X- setup code -X- */
  m = (theta > 1.) ? ((int) theta) : 1;
  ll = 0;
  p0 = q = p = exp(-theta);
  /* -X- end of setup code -X- */
#endif

} /* end of poisson_pprsc_init() */

// #define my  (DISTR.params[0])    /* shape */


static double f(long int k, double l_nu, double c_pm)
{
        return  exp(k * l_nu - _unur_factorialln(k) - c_pm);
}

int
unur_stdgen_sample_poisson_pprsc( struct unur_gen *gen )
     /* theta >= 10: Patchwork Rejection */
{
  /* -X- generator code -X- */

  static double        theta_last = -1.0;
  static long int      m,  k2, k4, k1, k5;
  static double        dl, dr, r1, r2, r4, r5, ll, lr, l_theta, c_pm,
    f1, f2, f4, f5, p1, p2, p3, p4, p5, p6;
  long int             Dk, X, Y;
  double               Ds, U, V, W;
  
  if (theta != theta_last)
			{                               /* set-up           */
				theta_last = theta;

 /* approximate deviation of reflection points k2, k4 from theta - 1/2      */
				Ds = sqrt(theta + 0.25);

 /* mode m, reflection points k2 and k4, and points k1 and k5, which     */
 /* delimit the centre region of h(x)                                    */
				m  = (long int) theta;
				k2 = (long int) ceil(theta - 0.5 - Ds);
				k4 = (long int)     (theta - 0.5 + Ds);
				k1 = k2 + k2 - m + 1L;
				k5 = k4 + k4 - m;

 /* range width of the critical left and right centre region             */
				dl = (double) (k2 - k1);
				dr = (double) (k5 - k4);

 /* recurrence constants r(k) = p(k)/p(k-1) at k = k1, k2, k4+1, k5+1    */
				r1 = theta / (double) k1;
				r2 = theta / (double) k2;
				r4 = theta / (double)(k4 + 1L);
				r5 = theta / (double)(k5 + 1L);

 /* reciprocal values of the scale parameters of expon. tail envelopes   */
				ll =  log(r1);                                   /* expon. tail left */
				lr = -log(r5);                                   /* expon. tail right*/

 /* Poisson constants, necessary for computing function values f(k)      */
				l_theta = log(theta);
				c_pm = m * l_theta - _unur_factorialln(m);

 /* function values f(k) = p(k)/p(m) at k = k2, k4, k1, k5               */
				f2 = f(k2, l_theta, c_pm);
				f4 = f(k4, l_theta, c_pm);
				f1 = f(k1, l_theta, c_pm);
				f5 = f(k5, l_theta, c_pm);

 /* area of the two centre and the two exponential tail regions          */
 /* area of the two immediate acceptance regions between k2, k4          */
				p1 = f2 * (dl + 1.0);                            /* immed. left      */
				p2 = f2 * dl         + p1;                       /* centre left      */
				p3 = f4 * (dr + 1.0) + p2;                       /* immed. right     */
				p4 = f4 * dr         + p3;                       /* centre right     */
				p5 = f1 / ll         + p4;                       /* expon. tail left */
				p6 = f5 / lr         + p5;                       /* expon. tail right*/
			}

		for (;;)
		{
 /* generate uniform number U -- U(0, p6)                                */
 /* case distinction corresponding to U                                  */
			if ((U = uniform() * p6) < p2)
			{     /* centre left      */

 /* immediate acceptance region R2 = [k2, m) *[0, f2),  X = k2, ... m -1 */
				if ((V = U - p1) < 0.0)  return(k2 + (long int)(U/f2));
 /* immediate acceptance region R1 = [k1, k2)*[0, f1),  X = k1, ... k2-1 */
				if ((W = V / dl) < f1 )  return(k1 + (long int)(V/f1));

 /* computation of candidate X < k2, and its counterpart Y > k2          */
 /* either squeeze-acceptance of X or acceptance-rejection of Y          */
				Dk = (long int)(dl * uniform()) + 1L;
				if (W <= f2 - Dk * (f2 - f2/r2))
				{             /* quick accept of  */
					return(k2 - Dk);                             /* X = k2 - Dk      */
				}
				if ((V = f2 + f2 - W) < 1.0)
				{                 /* quick reject of Y*/
					Y = k2 + Dk;
					if (V <= f2 + Dk * (1.0 - f2)/(dl + 1.0))
					{  /* quick accept of  */
						return(Y);                                 /* Y = k2 + Dk      */
					}
					if (V <= f(Y, l_theta, c_pm))  return(Y);       /* final accept of Y*/
				}
				X = k2 - Dk;
			}
			else if (U < p4)
			{                               /* centre right     */
	/*  immediate acceptance region R3 = [m, k4+1)*[0, f4), X = m, ... k4    */
				if ((V = U - p3) < 0.0)  return(k4 - (long int)((U - p2)/f4));
 /* immediate acceptance region R4 = [k4+1, k5+1)*[0, f5)                */
				if ((W = V / dr) < f5 )  return(k5 - (long int)(V/f5));

 /* computation of candidate X > k4, and its counterpart Y < k4          */
 /* either squeeze-acceptance of X or acceptance-rejection of Y          */
				Dk = (long int)(dr * uniform()) + 1L;
				if (W <= f4 - Dk * (f4 - f4*r4))
				{             /* quick accept of  */
					return(k4 + Dk);                             /* X = k4 + Dk      */
				}
				if ((V = f4 + f4 - W) < 1.0)
				{                 /* quick reject of Y*/
					Y = k4 - Dk;
					if (V <= f4 + Dk * (1.0 - f4)/ dr)
					{         /* quick accept of  */
						return(Y);                                 /* Y = k4 - Dk      */
					}
					if (V <= f(Y, l_theta, c_pm))  return(Y);       /* final accept of Y*/
				}
				X = k4 + Dk;
			}
			else
			{
				W = uniform();
				if (U < p5)
				{                                  /* expon. tail left */
					Dk = (long int)(1.0 - log(W)/ll);
					if ((X = k1 - Dk) < 0L)  continue;           /* 0 <= X <= k1 - 1 */
					W *= (U - p4) * ll;                          /* W -- U(0, h(x))  */
					if (W <= f1 - Dk * (f1 - f1/r1))  return(X); /* quick accept of X*/
				}
				else
				{                                         /* expon. tail right*/
					Dk = (long int)(1.0 - log(W)/lr);
					X  = k5 + Dk;                                /* X >= k5 + 1      */
					W *= (U - p5) * lr;                          /* W -- U(0, h(x))  */
					if (W <= f5 - Dk * (f5 - f5*r5))  return(X); /* quick accept of X*/
				}
			}

 /* acceptance-rejection test of candidate X from the original area      */
 /* test, whether  W <= f(k),    with  W = U*h(x)  and  U -- U(0, 1)     */
 /* log f(X) = (X - m)*log(theta) - log X! + log m!                         */
			if (log(W) <= X * l_theta - _unur_factorialln(X) - c_pm)  return(X);
		}


  /* -X- end of generator code -X- */
  
} /* end of unur_stdgen_sample_poisson_pprsc() */

