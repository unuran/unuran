/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      beta_gen.c                                                   *
 *                                                                           *
 *   Special generators for Beta distribution                                *
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

#include <unur_methods.h>
#include <unur_distribution.h>
#include <unur_distribution_lib.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static void beta_bc_init( struct unur_gen *gen );
inline static void beta_bb_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       par->data.cstd        /* data for parameter object         */
#define GEN       gen->data.cstd        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen) /* call for uniform prng             */

#define a     (DISTR.params[0])
#define b     (DISTR.params[1])
#define a_par (par->distr->data.cont.params[0])
#define b_par (par->distr->data.cont.params[1])
#define boundary_left  (DISTR.params[2])
#define boundary_right (DISTR.params[3])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_beta_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for gamma distribution                  */
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
  CHECK_NULL(par,0.);
  COOKIE_CHECK(par,CK_CSTD_PAR,0.);

  switch (par->variant) {
  case 0:  /* Rejection with log-logistic envelopes */  /* DEFAULT */
    if( (a_par > 1.) && (b_par > 1.) ) {
      _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_beta_bb );
      beta_bb_init( gen );
    }
    else {
      _unur_cstd_set_sampling_routine( par,gen,unur_stdgen_sample_beta_bc );
      beta_bc_init( gen );
    }
    return 1;
  case UNUR_STDGEN_INVERSION:   /* inversion method */
  default:
    return 0;
  }

} /* end of _unur_stdgen_beta_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Beta Distribution: Acceptance/Rejection from log-logistic hats for the    *
 *                    beta prime distribution                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Beta distribution with     *
 *               parameters a,b (a > 0, b > 0).                              *
 *               It combines algorithms bb (a > 1, b > 1) and                *
 *               bc (a<=1 or b<=1).                                          *
 *                                                                           *
 * REFERENCE : - R.C.H. Cheng (1978): Generating beta variates with          *
 *               nonintegral shape parameters,                               *
 *               Communications of the ACM 21, 317-322.                      *
 *                                                                           *
 * Implemented by E. Stadlober, R. Kremer, 1990                              *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

#define am      (GEN.gen_param[0])
#define bm      (GEN.gen_param[1])
#define al      (GEN.gen_param[2])
#define alnam   (GEN.gen_param[3])
#define be      (GEN.gen_param[4])
#define ga      (GEN.gen_param[5])
#define si      (GEN.gen_param[6])
#define rk1     (GEN.gen_param[7])
#define rk2     (GEN.gen_param[8])

inline static void
beta_bc_init( struct unur_gen *gen )
{
  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = 9;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  am = (a > b) ? a : b;
  bm = (a < b) ? a : b;
  al = am + bm;
  alnam = al * log(al/am) - 1.386294361;
  be = 1.0 / bm;
  si = 1.0 + am - bm;
  rk1 = si * (0.013888889 + 0.041666667 * bm) / (am * be - 0.77777778);
  rk2 = 0.25 + (0.5 + 0.25 / si) * bm;
  /* -X- end of setup code -X- */
} /* end of beta_bbbc_init() */

double 
unur_stdgen_sample_beta_bc(  struct unur_gen *gen )
     /* a <= 1. || b <= 1. */ 
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,w,y,z;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  while (1) {

    /* Step 1 */
    u1 = uniform();
    u2 = uniform();

    if (u1 < 0.5) {
      /* Step 2 */
      y = u1 * u2;
      z = u1 * y;

      if ((0.25 * u2 - y + z) >= rk1) 
	continue;  /* goto 1 */

      /* Step 5 */
      v = be * log(u1 / (1.0 - u1));
      if (v > 80.) {
	if (alnam < log(z)) continue;
	X = (am == a) ? 1. : 0.;
	break;
      }
      else {
	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z)) 
	  continue;  /* goto 1 */

	/* Step 6_a */
	X = (am != a) ? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
    else {
      /* Step 3 */
      z = u1 * u1 * u2;
      if (z < 0.25) {
	/* Step 5 */
	v = be * log(u1 / (1.0 - u1));
	if (v > 80.) {
	  X = (am == a) ? 1.0 : 0.0;
	  break;
	}

	w = am * exp(v);
	X = (am != a) ? bm / (bm + w) : w / (bm + w);
	break;
      }
      else {
	if (z >= rk2) continue;

	v = be * log(u1 / (1.0 - u1));
	if ( v > 80.) {
	  if (alnam < log(z)) continue;
	  X = (am == a) ? 1. : 0.;
	  break;
	}

	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z)) 
	  continue;  /* goto 1 */
	      
	/* Step 6_b */
	X = (am != a) ? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
  }
  /* -X- end of generator code -X- */
  
  return X;

} /* end of unur_stdgen_sample_beta_bc() */

/*---------------------------------------------------------------------------*/

inline static void
beta_bb_init( struct unur_gen *gen )
{
  if (GEN.gen_param == NULL) {
    GEN.n_gen_param = 9;
    GEN.gen_param = _unur_malloc(GEN.n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  am = (a < b) ? a : b;
  bm = (a > b) ? a : b;
  al = am + bm;
  be = sqrt((al - 2.0)/(2.0 * a * b - al));
  ga = am + 1.0 / be;
  /* -X- end of setup code -X- */
} /* end of beta_bbbc_init() */

double 
unur_stdgen_sample_beta_bb(  struct unur_gen *gen )
     /* (a > 1.0) && (b > 1.0) */ 
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,w,z,r,s,t;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  while (1) {
    /* Step 1 */
    u1 = uniform();
    u2 = uniform();
    v = be * log(u1 / (1.0 - u1));
    w = am * exp(v);
    z = u1 * u1 * u2;
    r = ga * v - 1.386294361;
    s = am + r - w;

    /* Step 2 */
    if (s + 2.609437912 < 5.0 * z) {
      /* Step 3 */
      t = log(z);
      if (s < t)
	/* Step 4 */
	if (r + al * log(al/(bm + w)) < t) 
	  continue;
    }

    /* Step 5 */
    X = (am == a) ? w / (bm + w) : bm / (bm + w);
    break;
  }
  /* -X- end of generator code -X- */

  return X;

} /* end of unur_stdgen_sample_beta_bb() */

#undef am
#undef bm
#undef al
#undef alnam
#undef be
#undef ga
#undef si
#undef rk1
#undef rk2

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#if 0

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
#include<math.h>

extern double drand (unsigned long *seed);

/******************************************************************
 *                                                                *
 * Beta Distribution - Stratified Rejection/Patchwork Rejection   *
 *                                                                *
 ******************************************************************
 * For parameters a < 1 , b < 1  and  a < 1 < b   or  b < 1 < a   *
 * the stratified rejection methods b00 and b01 of Sakasegawa are *
 * used. Both procedures employ suitable two-part power functions *
 * from which samples can be obtained by inversion.               *
 * If a > 1 , b > 1 (unimodal case) the patchwork rejection       *
 * method b1prs of Zechner/Stadlober is utilized:                 *
 * The area below the density function f(x) in its body is        *
 * rearranged by certain point reflections. Within a large center *
 * interval variates are sampled efficiently by rejection from    *
 * uniform hats. Rectangular immediate acceptance regions speed   *
 * up the generation. The remaining tails are covered by          *
 * exponential functions.                                         *
 * If (a-1)(b-1) = 0  sampling is done by inversion if either a   *
 * or b are not equal to one. If  a = b = 1  a uniform random     *
 * variate is delivered.                                          *
 *                                                                *
 ******************************************************************
 *                                                                *
 * FUNCTION :   - bsprc samples a random variate from the beta    *
 *                distribution with parameters  a > 0, b > 0.     *
 * REFERENCES : - H. Sakasegawa (1983): Stratified rejection and  *
 *                squeeze method for generating beta random       *
 *                numbers, Ann. Inst. Statist. Math. 35 B,        *
 *                291-302.                                        *
 *              - H. Zechner, E. Stadlober (1993): Generating     *
 *                beta variates via patchwork rejection,          *
 *                Computing 50, 1-18.                             *
 *                                                                *
 * SUBPROGRAMS: - drand(seed) ... (0,1)-Uniform generator with    *
 *                unsigned long integer *seed.                    *
 *              - b00(seed,a,b) ... Beta generator for a<1, b<1   *
 *              - b01(seed,a,b) ... Beta generator for a<1<b or   *
 *                                  b<1<a                         *
 *              - b1prs(seed,a,b) ... Beta generator for a>1, b>1 *
 *                with unsigned long integer *seed, double a, b.  *
 *                                                                *
 * Implemented by H. Zechner and F. Niederl, July 1994            *
 ******************************************************************/

static double b00(unsigned long *seed ,double a, double b)
{
	static double      a_last = 0.0,
				 b_last = 0.0;
	static double      a_, b_, c, t, fa, fb, p1, p2;
	double             U, V, X, Z;

	if (a != a_last || b != b_last) {
		a_last = a;
		b_last = b;

		a_ = a - 1.0;
		b_ = b - 1.0;
		c = (b * b_) / (a * a_);                            /* b(1-b) / a(1-a) */
		t = (c == 1.0) ? 0.5 : (1.0 - sqrt(c))/(1.0 - c);   /* t = t_opt       */
		fa = exp(a_ * log(t));
		fb = exp(b_ * log(1.0 - t));                        /* f(t) = fa * fb  */

		p1 = t/a;                                           /* 0 < X < t       */
		p2 = (1.0 - t)/b + p1;                              /* t < X < 1       */
	}

	for (;;) {
		if ((U = drand(seed) * p2) <= p1) {               /*  X < t  */
			Z = exp(log(U/p1) / a);  X = t*Z;
	/* squeeze accept:   L(x) = 1 + (1 - b)x                                 */
			if ((V = drand(seed) * fb) <= 1.0 - b_*X)  break;
	/* squeeze reject:   U(x) = 1 + ((1 - t)^(b-1) - 1)/t * x                */
			if (V <= 1.0 + (fb - 1.0)*Z) {
	/* quotient accept:  q(x) = (1 - x)^(b-1) / fb                           */
	if (log(V) <= b_ * log(1.0 - X))  break;
			}
		}
		else {                                                      /*  X > t  */
			Z = exp(log((U-p1)/(p2-p1)) / b);  X = 1.0 - (1.0 - t)*Z;
	/* squeeze accept:   L(x) = 1 + (1 - a)(1 - x)                           */
			if ((V = drand(seed) * fa) <= 1.0 - a_*(1.0 - X))  break;
	/* squeeze reject:   U(x) = 1 + (t^(a-1) - 1)/(1 - t) * (1 - x)          */
			if (V <= 1.0 + (fa - 1.0)*Z) {
	/* quotient accept:  q(x) = x^(a-1) / fa                                 */
	if (log(V) <= a_ * log(X))  break;
			}
		}
	}
	return(X);
}

static double b01(unsigned long *seed, double a, double b)
{
	static double      a_last = 0.0,
				 b_last = 0.0;
	static double      a_, b_, t, fa, fb, ml, mu, p1, p2;
	double             U, V, X, Z;

	if (a != a_last || b != b_last) {
		a_last = a;
		b_last = b;

		a_ = a - 1.0;
		b_ = b - 1.0;
		t = a_/(a - b);                   /* one step Newton * start value t   */
		fb = exp((b_ - 1.0) * log(1.0 - t));  fa = a - (a + b_)*t;
		t -= (t - (1.0 - fa) * (1.0 - t)*fb / b) / (1.0 - fa*fb);
		fa = exp(a_ * log(t));
		fb = exp(b_ * log(1.0 - t));                        /* f(t) = fa * fb  */
		if (b_ <= 1.0) {
			ml = (1.0 - fb) / t;                              /*   ml = -m1      */
			mu = b_ * t;                                      /*   mu = -m2 * t  */
		}
		else {
			ml = b_;
			mu = 1.0 - fb;
		}
		p1 = t/a;                                           /*  0 < X < t      */
		p2 = fb * (1.0 - t)/b + p1;                         /*  t < X < 1      */
	}

	for (;;) {
		if ((U = drand(seed) * p2) <= p1) {               /*  X < t  */
			Z = exp(log(U/p1) / a);  X = t*Z;
	/* squeeze accept:   L(x) = 1 + m1*x,  ml = -m1                          */
			if ((V = drand(seed)) <= 1.0 - ml*X)  break;
	/* squeeze reject:   U(x) = 1 + m2*x,  mu = -m2 * t                      */
			if (V <= 1.0 - mu*Z) {
	/* quotient accept:  q(x) = (1 - x)^(b-1)                                */
	if (log(V) <= b_ * log(1.0 - X))  break;
			}
		}
		else {                                                      /*  X > t  */
			Z = exp(log((U-p1)/(p2-p1)) / b);  X = 1.0 - (1.0 - t)*Z;
	/* squeeze accept:   L(x) = 1 + (1 - a)(1 - x)                           */
			if ((V = drand(seed) * fa) <= 1.0 - a_*(1.0 - X))  break;
	/* squeeze reject:   U(x) = 1 + (t^(a-1) - 1)/(1 - t) * (1 - x)          */
			if (V <= 1.0 + (fa - 1.0)*Z) {
	/* quotient accept:  q(x) = (x)^(a-1) / fa                                 */
	if (log(V) <= a_ * log(X))  break;
			}
		}
	}
	return(X);
}

static double f(double x, double a, double b, double m)
{
	return  exp(a*log(x/m) + b*log((1.0 - x)/(1.0 - m)));
}

double b1prs(unsigned long *seed, double p, double q)
{
	static double     p_last = 0.0,
				q_last = 0.0;
	static double     a, b, s, m, D, Dl, x1, x2, x4, x5, f1, f2, f4, f5,
				ll, lr, z2, z4, p1, p2, p3, p4;
	double            U, V, W, X, Y;

	if (p != p_last || q != q_last) {
		p_last = p;
		q_last = q;

		a = p - 1.0;
		b = q - 1.0;
		s = a + b;   m = a / s;
		if (a > 1.0 || b > 1.0)  D = sqrt(m * (1.0 - m) / (s - 1.0));

		if (a <= 1.0) {
			x2 = (Dl = m * 0.5);  x1 = z2 = 0.0;  f1 = ll = 0.0;
		}
		else {
			x2 = m - D;  x1 = x2 - D;
			z2 = x2 * (1.0 - (1.0 - x2)/(s * D));
			if (x1 <= 0.0 || (s - 6.0) * x2 - a + 3.0 > 0.0) {
	x1 = z2;  x2 = (x1 + m) * 0.5;
	Dl = m - x2;
			}
			else {
	Dl = D;
			}
			f1 = f(x1, a, b, m);
			ll = x1 * (1.0 - x1) / (s * (m - x1));            /* z1 = x1 - ll   */
		}
		f2 = f(x2, a, b, m);

		if (b <= 1.0) {
			x4 = 1.0 - (D = (1.0 - m) * 0.5);  x5 = z4 = 1.0;  f5 = lr = 0.0;
		}
		else {
			x4 = m + D;  x5 = x4 + D;
			z4 = x4 * (1.0 + (1.0 - x4)/(s * D));
			if (x5 >= 1.0 || (s - 6.0) * x4 - a + 3.0 < 0.0) {
	x5 = z4;  x4 = (m + x5) * 0.5;
	D = x4 - m;
			}
			f5 = f(x5, a, b, m);
			lr = x5 * (1.0 - x5) / (s * (x5 - m));            /* z5 = x5 + lr   */
		}
		f4 = f(x4, a, b, m);

		p1 = f2 * (Dl + Dl);                                /*  x1 < X < m    */
		p2 = f4 * (D  + D) + p1;                            /*  m  < X < x5   */
		p3 = f1 * ll       + p2;                            /*       X < x1   */
		p4 = f5 * lr       + p3;                            /*  x5 < X        */
	}

	for (;;) {
		if ((U = drand(seed) * p4) <= p1) {
	 /* immediate accept:  x2 < X < m, - f(x2) < W < 0                      */
			if ((W = U/Dl - f2) <= 0.0)  return(m - U/f2);
	 /* immediate accept:  x1 < X < x2, 0 < W < f(x1)                       */
			if (W <= f1)  return(x2 - W/f1 * Dl);
	 /* candidates for acceptance-rejection-test                            */
			V = Dl * (U = drand(seed));
			X = x2 - V;  Y = x2 + V;
	 /* squeeze accept:    L(x) = f(x2) (x - z2) / (x2 - z2)                */
			if (W * (x2 - z2) <= f2 * (X - z2))  return(X);
			if ((V = f2 + f2 - W) < 1.0) {
	 /* squeeze accept:    L(x) = f(x2) + (1 - f(x2))(x - x2)/(m - x2)      */
	if (V <= f2 + (1.0 - f2) * U)  return(Y);
	 /* quotient accept:   x2 < Y < m,   W >= 2f2 - f(Y)                    */
	if (V <= f(Y, a, b, m))  return(Y);
			}
		}
		else if (U <= p2) {
			U -= p1;
	 /* immediate accept:  m < X < x4, - f(x4) < W < 0                      */
			if ((W = U/D - f4) <= 0.0)  return(m + U/f4);
	 /* immediate accept:  x4 < X < x5, 0 < W < f(x5)                       */
			if (W <= f5)  return(x4 + W/f5 * D);
	 /* candidates for acceptance-rejection-test                            */
			V = D * (U = drand(seed));
			X = x4 + V;  Y = x4 - V;
	 /* squeeze accept:    L(x) = f(x4) (z4 - x) / (z4 - x4)                */
			if (W * (z4 - x4) <= f4 * (z4 - X))  return(X);
			if ((V = f4 + f4 - W) < 1.0) {
	 /* squeeze accept:    L(x) = f(x4) + (1 - f(x4))(x4 - x)/(x4 - m)      */
	if (V <= f4 + (1.0 - f4) * U)  return(Y);
		/* quotient accept:   m < Y < x4,   W >= 2f4 - f(Y)                   */
	if (V <= f(Y, a, b, m))  return(Y);
			}
		}
		else if (U <= p3) {                                   /*      X < x1  */
			Y = log(U = (U - p2)/(p3 - p2));
			if ((X = x1 + ll * Y) <= 0.0)  continue;            /*      X > 0!! */
			W = drand(seed) * U;
	 /* squeeze accept:    L(x) = f(x1) (x - z1) / (x1 - z1)                */
	 /*                    z1 = x1 - ll,   W <= 1 + (X - x1)/ll             */
			if (W <= 1.0 + Y)  return(X);
			W *= f1;
		}
		else {                                                /* x5 < X       */
			Y = log(U = (U - p3)/(p4 - p3));
			if ((X = x5 - lr * Y) >= 1.0)  continue;            /*      X < 1!! */
			W = drand(seed) * U;
	 /* squeeze accept:    L(x) = f(x5) (z5 - x) / (z5 - x5)                */
	 /*                    z5 = x5 + lr,   W <= 1 + (x5 - X)/lr             */
			if (W <= 1.0 + Y)  return(X);
			W *= f5;
		}
	 /* density accept:  f(x) = (x/m)^a ((1 - x)/(1 - m))^b                 */
		if (log(W) <= a*log(X/m) + b*log((1.0 - X)/(1.0 - m)))  return(X);
	}
}

double bsprc(unsigned long *seed, double a, double b)
{
	if (a  > 1.0) {
		if (b  > 1.0)  return(b1prs(seed, a, b));
		if (b  < 1.0)  return(1.0 - b01(seed, b, a));
		if (b == 1.0)
		{
			return(exp(log( drand(seed)) / a));
		}
	}
	if (a  < 1.0) {
		if (b  > 1.0)  return(b01(seed, a, b));
		if (b  < 1.0)  return(b00(seed, a, b));
		if (b == 1.0)
		{
			return(exp(log( drand(seed)) / a));
		}
	}
	if (a == 1.0) {
		if (b != 1.0)  return(1.0 - exp(log(drand(seed)) / b));
		if (b == 1.0)  return(drand(seed));
	}
	return 0.0;
}




#endif
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#undef a
#undef b
#undef a_par
#undef b_par
#undef boundary_left
#undef boundary_right
/*---------------------------------------------------------------------------*/
