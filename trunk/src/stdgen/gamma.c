
/*****************************************************************************
 *****************************************************************************/
/*---------------------------------------------------------------------------*/

#include <unur_stdgen.h>
#include <unur_math.h>

/*---------------------------------------------------------------------------*/

/**********************************************************************/
/*Gamma rng using Acceptance Rejection with log-logistic envelopes*/

double gammarand(double a, UNUR_URNG_TYPE urng)
{
 static double aa,bb,cc,a_in = -1.0;
 double u1,u2,v,r,z,gl;

 if (a != a_in)
   {
    a_in = a;
    aa = (a > 1.0)? sqrt(a + a - 1.0) : a;
    bb = a - 1.386294361;
    cc = a + aa;
   }
 for(;;)
   {
    u1 = _unur_call_urng_prt(urng);
    u2 = _unur_call_urng_prt(urng);
    v = log(u1 / (1.0 - u1)) / aa;
    gl = a * exp(v);
    r = bb + cc * v - gl;
    z = u1 * u1 * u2;
    if (r + 2.504077397 >= 4.5 * z) break;
    if (r >= log(z)) break;
   }
 return(gl);
}
/************************************************************/
