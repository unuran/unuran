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
  }

  /* -X- setup code -X- */

  /* make a normal variate generator (use default special generator) */
  NORMAL = unur_init( unur_cstd_new( unur_distr_normal(NULL,0) ));
  /* need same uniform random number generator as slash generator */
  NORMAL->urng = gen->urng;

  /* -X- end of setup code -X- */

} /* end of poisson_pdac_init() */


// #define my  (DISTR.params[0])    /* shape */

int
unur_stdgen_sample_poisson_pdac( struct unur_gen *gen )
     /* Theta >= 10: acceptance complement */
{
  /* -X- generator code -X- */
 static double theta_prev = -1.0, theta_alt = -1.0,
                         a0 =-0.5000000002, a1 = 0.3333333343, a2 =-0.2499998565,
                         a3 = 0.1999997049, a4 =-0.1666848753, a5 = 0.1428833286,
                         a6 =-0.1241963125, a7 = 0.1101687109, a8 =-0.1142650302,
                         a9 = 0.1055093006;
 static long fac[] = {
                                1,1,2,6,24,120,720,5040,40320,362880
                          };
 static long l;
 static double s,d;
 double t,g,theta_k;

 static double omega,b1,b2,c0,c1,c2,c3,c;
 double gx,gy,px,py,e,x,xx,delta,v;
 long sign;

 double u;
 long k;


  /* check arguments */
  CHECK_NULL(gen,0.); COOKIE_CHECK(gen,CK_DSTD_GEN,0.);


          if (theta_prev != theta)
                 {
        theta_prev = theta;
        s = sqrt(theta);
        d = 6.0 * theta * theta;
        l = (long int)(theta - 1.1484);
                 }
          t = unur_sample_cont(NORMAL);             /* Step N. Normal sample */
          g = theta + s * t;
          if (g >= 0.0)
                 {
        k = (long int)(g);
        if (k >= l) return(k);     /* Step I. Immediate acceptance */
        u = uniform();           /* Step S. Squeeze acceptance */
        theta_k = theta - k;
        if (d * u >= theta_k * theta_k * theta_k) return(k);
                 }
          if (theta_alt != theta)    /* Step P. Preparations for steps Q and H */
                 {
        theta_alt = theta;
        omega = 0.3989423 / s;
        b1 = 0.416666666667e-1 / theta;
        b2 = 0.3 * b1 * b1;
        c3 = 0.1428571 * b1 * b2;
        c2 = b2 - 15.0 * c3;
        c1 = b1 - 6.0 * b2 + 45.0 * c3;
        c0 = 1.0 - b1 + 3.0 * b2 - 15.0 * c3;
        c = 0.1069 / theta;
                 }
          if (g >= 0.0)
                 {
        /* FUNCTION F */
        if (k < 10)
          {
                px = -theta;
                py = exp(k * log(theta)) / fac[k];
          }
        else  /* k >= 10 */
          {
                delta = 0.83333333333e-1 / k;
                delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*k*k));
                v = (theta_k) / (double)k;
                if (fabs(v) > 0.25)
                  {
                        px = k * log(1.0 + v) - theta_k - delta;
                  }
                else
                  {
                        px = k * v * v;
                        px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
                                 a4)*v+a3)*v+a2)*v+a1)*v+a0;
                        px -= delta;
                  }
                py = 0.3989422804 / sqrt((double)k);
          }
        x = (0.5 - theta_k) / s;
        xx = x * x;
        gx = -0.5 * xx;
        gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
                                        /* Step Q. Quotient acceptance */
        if (gy * (1.0 - u)  <= py * exp(px - gx)) return(k);
                 }
          for(;;)
                 {
        do
          {
                e = - log(uniform()); /* Step E. Double exponential sample */
                u = uniform();
                u = u + u - 1.0;
                sign = (u < 0.0)? -1 : 1;
                t = 1.8 + e * sign;
          }
        while (t <= -0.6744);
        k = (long int)(theta + s * t);
        theta_k = theta-k;
        /* FUNCTION F */
        if (k < 10)
          {
                px = -theta;
                py = exp(k * log(theta)) / fac[k];
          }
        else  /* k >= 10 */
          {
                delta = 0.83333333333e-1 / (double)k;
                delta = delta - 4.8*delta*delta*delta*(1.0-1.0/(3.5*k*k));
                v = (theta_k) / (double)k;
                if (fabs(v) > 0.25)
                  {
                        px = k * log(1.0 + v) - theta_k - delta;
                  }
                else
                  {
                        px = k * v * v;
                                  px *= ((((((((a9*v+a8)*v+a7)*v+a6)*v+a5)*v+
                                 a4)*v+a3)*v+a2)*v+a1)*v+a0;
                        px -= delta;
                  }
                py = 0.3989422804 / sqrt((double)k);
          }
        x = (0.5 - theta_k) / s;
        xx = x * x;
        gx = -0.5 * xx;
        gy = omega * (((c3 * xx + c2) * xx + c1) * xx + c0);
                 /* Step H. Hat acceptance */
        if (c * sign * u <= py * exp(px + e) - gy * exp(gx + e)) return(k);
                 }



  /* -X- end of generator code -X- */
  
} /* end of unur_stdgen_sample_poisson_pdac() */

/*---------------------------------------------------------------------------*/
#undef NORMAL
/*---------------------------------------------------------------------------*/





