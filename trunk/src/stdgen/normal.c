/*****************************************************************************
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_stdgen.h>
#include <unur_math.h>

/*---------------------------------------------------------------------------*/

/***************************************************************
 * NORMAL Distribution - NBM (Sinus-Cosinus or Box/Muller      *
 *                       Method)                               *
 ***************************************************************
 * FUNCTION: -    nsc samples a random number from the         *
 *                standard normal distribution                 *
 * REFERENCE: -   G.E.P. Box and M.E. Muller: A note on the    *
 *                generation of random normal deviates.        *
 *                Annals of Mathem. Stat. 29, (1958), 610-611. *
 ***************************************************************/

double nbm(UNUR_URNG_TYPE urng)
{
  double u,v,s;
  static double x2;
  static double two_pi = 2.0 * M_PI;
  static char f = 1;

  f = -f;
  if (f > 0) return(x2);
  u = _unur_call_urng_prt(urng);
  v = _unur_call_urng_prt(urng);
  s = sqrt(-2.0 * log(u));
  x2 = s * sin(two_pi * v);
  return(s * cos(two_pi * v));
} /* end of nbm() */

/*------------------------------------------------------------------*/

/*****************************************************************************

 NOT MODIFIED BELOW THIS LINE !!!!!									      

 *****************************************************************************/

#if 0

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *****************************************************************************/
/*******************************************************************
  
  "naive" ratio-of-uniforms, MBR, no squeezes
  e.g.\ Fishman (1996), Monte Carlo. Concepts, algorithms, and applications,
  p.~181, algorithm ROU

  ******************************************************************/

double nnquo(void)
{
  double u,v,x;

  while (1) {
    u=uniform();
    if (u==0.0) u = 1.0;
    v=(uniform()-0.5)*0.857763885*2;
    x=v/u;
    if (x*x<=-4*log(u)) return(x);
  }
} /* end of nnquo() */

/*------------------------------------------------------------------*/

/*******************************************************************
  
  ratio-of-uniforms, MBR, optimized by using squeezes
  Barabesi (1993), Random variate generation 
  by using the ratio-of-uniforms method, p. 133

  ******************************************************************/

double nquo(void)
{
  double r,rnormal,w;

  while( 1 ) {
    r=uniform();
    rnormal = (2.101083837941101*uniform()-1.050541918970551)/sqrt(r);
    w=rnormal*rnormal;
    if(4. - 4.186837275258269*r < w) {
      if(1.5/r - 0.920558458320164 < w) continue;
      if(-3.*log(r) < w ) continue;
    }
    return( rnormal );
  }
} /* end of nquo() */

/*------------------------------------------------------------------*/
    
/*******************************************************************
  
  ratio-of-uniforms, quadratic bounding curves
  Leva (1992), Algorithm 712; a normal random number generator,
  ACM TOMS 18(4), pp.~454--455

  ******************************************************************/

#define abs(x)   (((x)<0) ? (-(x)) : (x))

#define S    0.449871
#define T   -0.386595
#define A    0.19600
#define B    0.25472
#define RA   0.27597
#define RB   0.27846

double nleva(void)
{
  double u,v,x,y,q;

  while(1) {
    u = uniform();
    v = uniform();
    v = 1.7156 * (v - 0.5);
    x = u - S;
    y = abs(v) - T;
    q = x * x + y * (A * y - B * x);
    if( q < RA ) return(v/u);
    if( q > RB ) continue;
    if( v*v > -4.*log(u)*u*u) continue;
    return(v/u);
  }
} /* end of nleva() */

#undef S
#undef T  
#undef A  
#undef B  
#undef RA 
#undef RB 

#undef abs

/*------------------------------------------------------------------*/

/*******************************************************************

  Polarmethod with rejection
  Marsaglia (1962)

  ******************************************************************/

double npol(void)
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

/*------------------------------------------------------------------*/

/***************************************************************
 * NORMAL Distribution- KR (not subject to C-rand library)     *
 ***************************************************************
 * FUNCTION:  KR    samples a random number from the standard  *
 *                  normal distribution                        *
 * REFERENCE: -     Kinderman A.J., Ramage J.G.:               *
 *                  Computer Generation of Normal Random       *
 *                  Variables, in: Journal of the American     *
 *                  Statistical Association, Dec. 1976,        *
 *                  Volume 71, Number 356, pp. 893 - 898.      *
 *                                                             *
 * Implemented by:  M. Lehner April 1992                       *
 ***************************************************************/

#define XI 2.216035867166471
#define PIhochK 0.3989422804

#define min(v,w) ( ((v)<(w)) ? (v) : (w) )
#define max(v,w) ( ((v)>(w)) ? (v) : (w) )

double nkr(void)
{
   double t, u, v, w, x, z;

   u=  uniform();
   if(u<0.884070402298758)
      {v=  uniform();
       x=XI * (1.131131635444180*u+v-1);
      }

   else
   if(u>=0.973310954173898)
      {do{v=  uniform();
          w=  uniform();
          if(w==0) {w+=0.000001;} /* ?????? printf("%lf ",w);} */
          t=XI*XI/2-log(w);
         }while((v*v*t) > (XI*XI/2));
       if(u<0.986655477086949)
          x=pow(2*t,0.5);
       else
	 x= -pow(2*t,0.5);
      }

   else
   if(u>=0.958720824790463)
      {do{v=  uniform();
          w=  uniform();
          z=v-w;
          t=XI-0.630834801921960*min(v,w);
         }while(max(v,w)>0.755591531667601 &&
                0.034240503750111*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
       if(z<0) x=t;
       else x=(-t);
      }

   else
   if(u>=0.911312780288703)
      {do{v=  uniform();
          w=  uniform();
          z=v-w;
          t=0.479727404222441+1.105473661022070*min(v,w);
         }while(max(v,w)>0.872834976671790 &&
                0.049264496373128*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
       if(z<0) x=t;
       else x=(-t);
      }

   else
      {do{v=  uniform();
          w=  uniform(); z=v-w;
          t=0.479727404222441-0.595507138015940*min(v,w);
         }while(max(v,w)>0.805777924423817 &&
                0.053377549506886*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
       if(z<0) x=t;
       else x=(-t);
      }
   return x;
  } /* end of nkr() */

#undef XI
#undef PIhochK 
#undef min
#undef max

/*------------------------------------------------------------------*/

/*******************************************************************

 * FUNCTION:  KR    samples a random number from the standard  *
 *                  normal distribution                        *
 * REFERENCE: -     Kinderman A.J., Ramage J.G.:               *
 *                  Computer Generation of Normal Random       *
 *                  Variables, in: Journal of the American     *
 *                  Statistical Association, Dec. 1976,        *
 *                  Volume 71, Number 356, pp. 893 - 898.      *
 *                                                             *
 Zweite Version ...

  ******************************************************************/

#define abs(x)   ((x)<0) ? (-(x)) : (x)

double f(double t)
{
  return( (t>=0.) ? 0.39894228040*exp(-t*t/2.0)-0.180025191068563*(2.216035867166471-abs(t)) : 0. );
}

#undef abs

/*---------------------------------------------------------------------------*/

#define XI  2.216035867166471
#define XI2 4.91081496456825

double nkira(void)
{
  double u,v,w,maxvw,minvw,z,t;

  /*s1*/
  u = uniform();
  if (u<0.884070402298758) {
    v = uniform();
    return( XI*(1.131131635444180*u+v-1) );
  }
  /*s2*/
  if (u<0.973310954173898) goto s4;
  /*s3*/
 s3:
  v = uniform();
  w = uniform();
  t=XI2/2-log(w);
  if (v*v*t>XI2/2) goto s3;
  else {
    if (u<0.986655477086949) return( sqrt(2*t) );
    else return( -sqrt(2*t) );
  } /* endif */
  /*s4*/
 s4:
  if (u<0.958720824790463) goto s6;
  /*s5*/
 s5:
  v = uniform();
  w = uniform();
  if (v<=w) {
    minvw=v;
    maxvw=w;
  }
  else {
    minvw=w;
    maxvw=v;
  } /* endif */
  z=v-w;
  t=XI-0.630834801921960*minvw;
  if (maxvw<=0.755591531667601) goto s9;
  if (0.034240503750111*fabs(z)<=f(t)) goto s9;
  else goto s5;
  /*s6*/
 s6:
  if (u<=0.911312780288703) goto s8;
  /*s7*/
 s7:
  v = uniform();
  w = uniform();
  if (v<=w) {
    minvw=v;
    maxvw=w;
  }
  else {
    minvw=w;
    maxvw=v;
  } /* endif */
  z=v-w;
  t=0.479727404222441+1.105473661022070*minvw;
  if (maxvw<=0.872834976671790) goto s9;
  if (0.049264496373128*fabs(z)<=f(t)) goto s9;
  else goto s7;
  /*s8*/
 s8:
  v = uniform();
  w = uniform();
  if (v<=w) {
    minvw=v;
    maxvw=w;
  }
  else {
    minvw=w;
    maxvw=v;
  } /* endif */
  z=v-w;
  t=0.479727404222441-0.595507138015940*minvw;
  if (maxvw<=0.805577924423817) goto s9;
  if (0.053377549506886*fabs(z)<=f(t)) goto s9;
  else goto s8;
  /*s9*/
 s9:
  if (z<0) return(t);
  else return(-t);

} /* end of nkira() */

#undef XI  
#undef XI2 

/*---------------------------------------------------------------------------*/

/***************************************************************
 * NORMAL Distribution-NACR (acceptance complement ratio o.u.) *
 ***************************************************************
 * FUNCTION:  nacr  samples a random number from the standard  *
 *                  normal distribution                        *
 * REFERENCE: -     W. Hoermann and G. Derflinger:             *
 *                  The ACR Methodfor generating normal        *
 *                  random variables.                          *
 *                  OR Spektrum 12 (1990), 181-185.            *
 * SUBPROGRAMS: -   uniform() ... (0,1)--uniform generator,    *
 *                  contained in cranduni.c                    *
 *                                                             *
 ***************************************************************/

/*------------------------------------------------------------------*/

static double nt(void)
{  
  double rn,x,y,g;
  static double
    as=0.8853395638,bs=0.2452635696,cs=0.2770276848,b=0.5029324303,
    x0=0.4571828819,ym=0.187308492,s=0.7270572718,t=0.03895759111;

  while (1) {
    x = uniform();
    y = ym * uniform();
    g = x0-s*x-y;
    if (g>0) rn=2+y/x;
    else {
      x=1-x;
      y=ym-y;
      rn= -(2+y/x);
    }
    
    if ((y-as+x)*(cs+x)+bs<0) return(rn);
    else if (y<x+t) if (rn*rn<4*(b-log(x))) return(rn);
  }
}

/*---------------------------------------------------------------------------*/

double nacr(void)
{
  double  rn,xx,y,z;
  static double c1=1.448242853,c2=3.307147487,c3=1.46754004,
    d1=1.036467755,d2=5.295844968,d3=3.631288474,
    hm=0.483941449,zm=0.107981933,hp=4.132731354,
    zp=18.52161694,phln=0.4515827053,
    hm1=0.516058551,hp1=3.132731354,hzm=0.375959516,
    hzmp=0.591923442/*,zhm=0.967882898*/;

  y = uniform();
  if (y>hm1) 
    return (hp*y-hp1);
  else if (y<zm) {  
    rn=zp*y-1;
    if (rn>0) return(1+rn);
    else return(-1+rn);
  } 
  else if (y<hm) {  
    rn = uniform();
    rn=rn-1+rn;
    if (rn>0) z=2-rn;
    else z= -2-rn;
    if((c1-y)*(c3+fabs(z))<c2) return( z );
    else {  
      xx=rn*rn;
      if ((y+d1)*(d3+xx)<d2)              return( rn );
      else if (hzmp-y<exp(-(z*z+phln)/2)) return( z );
      else if (y+hzm<exp(-(xx+phln)/2))   return( rn );
      else                                return( nt() );
    }					
  }
  else return( nt() );
}

/*---------------------------------------------------------------------------*/

/*******************************************************************
  
  ziggurat method
  @Article{Marsaglia;Tsang:1984a,
  author =       {Marsaglia, G. and
                  Tsang, W. W.},
  year =         1984,
  title =        {A fast, easily implemented method for sampling from
                  decreasing or symmetric unimodal density functions},
  journal =      {SIAM J. Sci. Statist. Comput.},
  volume =       5,
  pages =        {349--359},
  }

  DO NOT USE!!!

  ******************************************************************/

#define abs(x)   (((x)<0) ? (-(x)) : (x))

/*---------------------------------------------------------------------------*/

double nzig(void)
{
  static const double
    aa = 12.37586, b = .4878992, c = 12.67706, rmax = .4656613e-9,
    c1 = .9689279, c2 = 1.301198, pc = .1958303e-1, xn = 2.776994;
  static const double v[65] = { .3409450, .4573146, .5397792, .6062427, .6631690,
			 .7136974, .7596124, .8020356, .8417227, .8792102, .9148948,
			 .9490791, .9820005, 1.013848, 1.044780, 1.074924, 1.104391,
			 1.133273, 1.161653, 1.189601, 1.217181, 1.244452, 1.271463,
			 1.298265, 1.324901, 1.351412, 1.377839, 1.404221, 1.430593,
			 1.456991, 1.483452, 1.510012, 1.536706, 1.563571, 1.590645,
			 1.617968, 1.645579, 1.673525, 1.701850, 1.730604, 1.759842,
			 1.789622, 1.820009, 1.851076, 1.882904, 1.915583, 1.949216,
			 1.983924, 2.019842, 2.057135, 2.095992, 2.136644, 2.179371,
			 2.224517, 2.272518, 2.323934, 2.379500, 2.440222, 2.507511,
			 2.583466, 2.671391, 2.776994, 2.776994, 2.776994, 2.776994 };

  int j;
  double u,rnor,x,y,s;

  /* fast part */
  u = uniform();
/*    j = ((int) (64*uniform())) % 64; */
  j = ((int) (64*u)) % 64;
  rnor = (2*u-1) * v[j+1];
  if( abs(rnor) <= v[j] ) return rnor;

  /* slow part */
  x = ( abs(rnor) - v[j] ) / (v[j+1] - v[j] );
  y = uniform();
  s = x+y;
  if( s > c2 ) goto s11;
  if( s <= c1 ) return rnor;
  if( y > c - aa*exp( -0.5*(b-b*x)*(b-b*x)) ) goto s11;
  if( (exp( -0.5*v[j+1]*v[j+1]) + y*pc/v[j+1]) <= exp( -0.5*rnor*rnor ) ) return rnor;
  
  /* tail part */
 s22:
  x = .3601016 * log(uniform());
  if( -2. * log(uniform()) <= x*x ) goto s22;
 s33:
  return ((rnor>0) ? xn-x : x-xn);
 s11:
  return ((rnor>0) ? b-b*x : b*x-b);

} /* end of nzig() */

#undef abs

/*------------------------------------------------------------------*/

#endif


