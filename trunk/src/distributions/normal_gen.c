/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      normal_gen.c                                                 *
 *                                                                           *
 *   Special generators for Normal distribution                              *
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
/* Prototypes for special generators                                         */

inline static double nbm(UNUR_URNG_TYPE urng);
/* Box Muller method                                                         */

inline static double npol(UNUR_URNG_TYPE urng);
/* Polarmethod with rejection                                                */

inline static double nnquo(UNUR_URNG_TYPE urng);
/* "Naive" ratio-of-uniforms method                                          */

inline static double nquo(UNUR_URNG_TYPE urng);
/* Ratio-of-uniforms method with squeeze                                     */

inline static double nleva(UNUR_URNG_TYPE urng);
/* Ratio-of-uniforms method  with quadratic bounding curves                  */

inline static double nkr(UNUR_URNG_TYPE urng);
/* Kindermann-Ramage method                                                  */

inline static double nacr(UNUR_URNG_TYPE urng);
/* Acceptance-complement ratio                                               */

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define GEN        gen->data.cstd
#define uniform()  (_unur_call_urng_prt(urng))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  get special sampling routine for distribution                          **/
/**                                                                         **/
/*****************************************************************************/

_UNUR_SAMPLING_ROUTINE_CONT *
_unur_stdgen_normal_get_routine(unsigned variant)
     /*----------------------------------------------------------------------*/
     /* get pointer to sampling routine                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   variant ... variant of special generator                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to sampling routine                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  switch (variant) {
  case 0: /* Default */
    return unur_stdgen_sample_normal_bm;    /* Box-Muller method */
  case 1:
    return unur_stdgen_sample_normal_pol;   /* Polarmethod with rejection */
  case 2:
    return unur_stdgen_sample_normal_nquo;  /* "Naive" ratio-of-uniforms */
  case 3:
    return unur_stdgen_sample_normal_quo;   /* Ratio-of-uniforms with squeeze */
  case 4:
    return unur_stdgen_sample_normal_leva;  /* Ratio-of-uniforms with quadratic bounding curves */
  case 5:
    return unur_stdgen_sample_normal_kr;    /* Kindermann-Ramage method */
  case 6:
    return unur_stdgen_sample_normal_acr;   /* Acceptance-complement ratio */
  case UNUR_STDGEN_INVERSION:
  default:
    return NULL;
  }

} /* end of _unur_stdgen_normal_get_routine() */

/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO

const char *
_unur_stdgen_normal_routinename(void *routine)
     /*----------------------------------------------------------------------*/
     /* get name of sampling routine                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   routine ... pointer to sampling routine                            */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to name of sampling routine                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define routinename(rn) if (routine==(rn)) return #rn

  routinename( unur_stdgen_sample_normal_bm );
  routinename( unur_stdgen_sample_normal_pol );
  routinename( unur_stdgen_sample_normal_nquo );
  routinename( unur_stdgen_sample_normal_quo );
  routinename( unur_stdgen_sample_normal_leva );
  routinename( unur_stdgen_sample_normal_kr );
  routinename( unur_stdgen_sample_normal_acr );

  return NULL;

} /* end of _unur_stdgen_normal_routinename() */

#endif

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Wrapper for special generators                                         **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_bm( struct unur_gen *gen )
     /* Box-Muller method                                                    */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nbm(gen->urng));

} /* end of unur_stdgen_sample_normal_bm() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_pol( struct unur_gen *gen )
     /* Polarmethod with rejection                                           */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * npol(gen->urng));

} /* end of unur_stdgen_sample_normal_pol() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_nquo( struct unur_gen *gen )
     /* "Naive" ratio-of-uniforms method                                     */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nnquo(gen->urng));

} /* end of unur_stdgen_sample_normal_nquo() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_quo( struct unur_gen *gen )
     /* Ratio-of-uniforms method with squeeze                                */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nquo(gen->urng));

} /* end of unur_stdgen_sample_normal_quo() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_leva( struct unur_gen *gen )
     /* Ratio-of-uniforms method  with quadratic bounding curves             */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nleva(gen->urng));

} /* end of unur_stdgen_sample_normal_leva() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_kr( struct unur_gen *gen )
     /* Kindermann-Ramage method                                             */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nkr(gen->urng));

} /* end of unur_stdgen_sample_normal_kr() */

/*---------------------------------------------------------------------------*/

double unur_stdgen_sample_normal_acr( struct unur_gen *gen )
     /* Acceptance-complement ratio                                          */
{
  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_CSTD_GEN,0.);

  return (GEN.pdf_param[0] + GEN.pdf_param[1] * nacr(gen->urng));

} /* end of unur_stdgen_sample_normal_acr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

inline static double
nbm(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Sinus-Cosinus or Box/Muller Method                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * This method is based on the transformation                                *
 * x = sqrt(-2 ln(u)) cos(2pi*v), y = sqrt(-2 ln(u)) sin(2pi*v)              *
 * which converts two independent (0,1)-Uniforms u and v to two              *
 * independent standard Normal variates x and y.                             *
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - G.E.P. Box, M.E. Muller (1958):                             *
 *               A note on the generation of random normal deviates,         *
 *               Annals Math. Statist. 29, 610-611.                          *
 *                                                                           *
 * Implemented by W. Hoermann, April 1992                                    *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/
{
  double u,v,s;
  static double x2;
  static double two_pi = 2.0 * M_PI;
  static char f = 1;

  f = -f;
  if (f > 0) return(x2);
  u = uniform();
  v = uniform();
  s = sqrt(-2.0 * log(u));
  x2 = s * sin(two_pi * v);
  return(s * cos(two_pi * v));
} /* end of nbm() */

/*---------------------------------------------------------------------------*/

inline static double
npol(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Polarmethod with rejection                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the standard                   *
 *               Normal distribution  N(0,1).                                *
 *                                                                           *
 * REFERENCE:  - G. Marsaglia (1962): Improving the Polar Method for         *
 *               Generating a Pair of Random Variables,                      *
 *               Boeing Sci. Res. Lab., Seattle, Washington.                 *
 *                                                                           *
 * Implemented by W. Hoermann, April 1992                                    *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/
{
  double s,x,y,tmp;
  static double store;
  static char f = 1;

  f = -f;
  if (f > 0) return( store );

  while(1) {
    x = 2. * uniform() - 1.;
    y = 2. * uniform() - 1.;
    s = x*x + y*y;
    if( s < 1. ) {
      tmp = sqrt( -2. * log(s) / s );
      store = y*tmp;
      return (x*tmp);
    }
  }
} /* end of npol() */

/*---------------------------------------------------------------------------*/

inline static double
nnquo(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: "Naive" Ratio of uniforms                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - A.J. Kindermann, F.J.Monahan (1977): Computing generation   *
 *               of random variables using the ratio of uniform deviates,    *
 *               ACM TOMS 3(3), 257-260                                      *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/
{
  double u,v,x;

  while (1) {
    u = uniform();
    if (u==0.) u = 1.;
    v = (uniform() - 0.5) * 0.857763885 * 2;
    x = v/u;
    if (x*x <= -4. * log(u)) 
      return x;
  }
} /* end of nnquo() */

/*------------------------------------------------------------------*/

inline static double
nquo(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Ratio of uniforms with squeeze                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - L. Barabesi (1993): Random variate generation               *
 *               by using the ratio-of-uniforms method, p. 133               *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/
{
  double r,rnormal,w;

  while (1) {
    r = uniform();
    rnormal = (2.101083837941101 * uniform() - 1.050541918970551) / sqrt(r);
    w = rnormal * rnormal;
    if (4. - 4.186837275258269 * r < w) {
      if (1.5/r - 0.920558458320164 < w) 
	continue;
      if (-3.*log(r) < w )
	continue;
    }
    return rnormal;
  }
} /* end of nquo() */

/*------------------------------------------------------------------*/

inline static double
nleva(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Ratio of uniforms with quadratic bounding curves     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - J.L. Leva (1992):                                           *
 *               Algorithm 712; a normal random number generator,            *
 *               ACM TOMS 18(4), 454-455                                     *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/
{
#define S    0.449871
#define T   -0.386595
#define A    0.19600
#define B    0.25472
#define RA   0.27597
#define RB   0.27846

  double u,v,x,y,q;

  while (1) {
    u = uniform();
    v = uniform();
    v = 1.7156 * (v - 0.5);
    x = u - S;
    y = fabs(v) - T;
    q = x * x + y * (A * y - B * x);
    if( q < RA ) return(v/u);
    if( q > RB ) continue;
    if (v*v > -4.*log(u)*u*u) continue;
    return(v/u);
  }

#undef S
#undef T  
#undef A  
#undef B  
#undef RA 
#undef RB 
} /* end of nleva() */

/*------------------------------------------------------------------*/

inline static double
nkr(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Kindermann-Ramage (patchwork) method                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - Kinderman A.J., Ramage J.G. (1976):                         *
 *               Computer Generation of Normal Random Variables,             *
 *               J. Am. Stat. Assoc. 71(356), 893 - 898.                     *
 *                                                                           *
 * Implemented by:  M. Lehner April 1992                                     *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/
{
#define XI 2.216035867166471
#define PIhochK 0.3989422804

   double t, u, v, w, x, z;

   u = uniform();

   if (u < 0.884070402298758) {
     v = uniform();
     x = XI * (1.131131635444180 * u + v - 1.);
   }

   else if (u >= 0.973310954173898) {
     do {
       v = uniform();
       w = uniform();
       if (w==0.) { t=0.; continue; }    /** TODO: is this necessary ?? **/
       t = XI * XI/2. - log(w);
     } while ( (v*v*t) > (XI*XI/2.) );
     x = (u < 0.986655477086949) ? pow(2*t,0.5) : -pow(2*t,0.5);
   }

   else if (u>=0.958720824790463) {
     do {
       v = uniform();
       w = uniform();
       z = v - w;
       t = XI - 0.630834801921960 * min(v,w);
     } while (max(v,w) > 0.755591531667601 &&
	      0.034240503750111 * fabs(z) > (PIhochK * exp(t*t/(-2.)) - 0.180025191068563*(XI-fabs(t))) );
     x = (z<0) ? t : -t;
   }

   else if (u>=0.911312780288703) {
     do {
       v = uniform();
       w = uniform();
       z = v - w;
       t = 0.479727404222441 + 1.105473661022070 * min(v,w);
     } while (max(v,w) > 0.872834976671790 &&
	      0.049264496373128*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
     x = (z<0) ? t : -t;
   }

   else {
     do {
       v = uniform();
       w = uniform(); 
       z = v - w;
       t = 0.479727404222441 - 0.595507138015940 * min(v,w);
     } while (max(v,w)>0.805777924423817 &&
	      0.053377549506886*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
     x = (z<0) ? t : -t;
   }

   return x;

#undef XI
#undef PIhochK 
#undef min
#undef max
} /* end of nkr() */

/*------------------------------------------------------------------*/

inline static double
nacr(UNUR_URNG_TYPE urng)
/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Acceptance-complement ratio                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - W. Hoermann and G. Derflinger (1990):                       *
 *               The ACR Methodfor generating normal random variables,       *
 *               OR Spektrum 12 (1990), 181-185.                             *
 *                                                                           *
 * Implemented by:  M. Lehner April 1992                                     *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/
{
  double  rn,x,y,z;

  static double c1=1.448242853,c2=3.307147487,c3=1.46754004,
    d1=1.036467755,d2=5.295844968,d3=3.631288474,
    hm=0.483941449,zm=0.107981933,hp=4.132731354,
    zp=18.52161694,phln=0.4515827053,
    hm1=0.516058551,hp1=3.132731354,hzm=0.375959516,
    hzmp=0.591923442/*,zhm=0.967882898*/;

  static double
    as=0.8853395638,bs=0.2452635696,cs=0.2770276848,b=0.5029324303,
    x0=0.4571828819,ym=0.187308492,s=0.7270572718,t=0.03895759111;

  y = uniform();

  if (y>hm1) 
    return (hp*y-hp1);

  else if (y<zm) {  
    rn = zp*y-1;
    return ( (rn>0) ? (1+rn) : (-1+rn));
  } 

  else if (y<hm) {  
    rn = uniform();
    rn = rn-1+rn;
    z = (rn>0) ? 2-rn : -2-rn;
    if ((c1-y)*(c3+fabs(z))<c2) return z;
    else {  
      x = rn*rn;
      if ((y+d1)*(d3+x)<d2)
	return rn;
      else if (hzmp-y<exp(-(z*z+phln)/2)) 
	return z;
      else if (y+hzm<exp(-(x+phln)/2))
	return rn;
    }
  }

  while (1) {
    x = uniform();
    y = ym * uniform();
    z = x0 - s*x - y;
    if (z>0) 
      rn = 2+y/x;
    else {
      x = 1-x;
      y = ym-y;
      rn = -(2+y/x);
    }
    if ((y-as+x)*(cs+x)+bs<0) 
      return rn;
    else if (y<x+t)
      if (rn*rn<4*(b-log(x)))
	return rn;
  }

} /* end of nacr() */

/*---------------------------------------------------------------------------*/

