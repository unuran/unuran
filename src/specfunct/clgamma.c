/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Logarithm of complex gamma function                                       */
/*                                                                           */
/*---------------------------------------------------------------------------*/

// Complex log-gamma function
// It returns log(|w|)=Re(log(w)) for w=gamma(z) with complex argument z=x+iy
// Taken from C++ code 
//  cgamma.cpp -- Complex gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
// It is the same as the fortran code of the R package 'fAsianOptions' 


#include <utils/unur_fp_source.h>
#include <utils/unur_math_source.h>
#include <utils/umath.h>

/*---------------------------------------------------------------------------*/

double
_unur_clgamma (double x, double y)
/*---------------------------------------------------------------------------*/
/* Logarithm of complex gamma function: log(Gamma(|z|)) = Re(log(Gamma(z))   */
/*                                                                           */
/* parameters:                                                               */
/*   x ... real part of argument z                                           */
/*   y ... imaginary part of argument z                                      */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* references:                                                               */
/*   S. Zhang & J. Jin: "Computation of Special Functions" (Wiley, 1996).    */
/*   online: http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html      */
/*                                                                           */

/* // Complex log-gamma function */
/* // It returns log(|w|)=Re(log(w)) for w=gamma(z) with complex argument z=x+iy */
/* // Taken from C++ code  */
/* //  cgamma.cpp -- Complex gamma function. */
/* //      Algorithms and coefficient values from "Computation of Special */
/* //      Functions", Zhang and Jin, John Wiley and Sons, 1996. */
/* // It is the same as the fortran code of the R package 'fAsianOptions'  */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Original FORTRAN code by S. Zhang and J. Jin.                             */
/* http://iris-lee3.ece.uiuc.edu/~jjin/routines/routines.html                */
/*                                                                           */
/* C All the programs and subroutines contained in this archive are          */
/* C copyrighted. However, we give permission to the user who downloads      */
/* C these routines to incorporate any of these routines into his or         */
/* C her programs provided that the copyright is acknowledged.               */
/*                                                                           */
/* C Contact Information                                                     */
/* C Email: j-jin1@uiuc.edu                                                  */
/* C Phone: (217) 244-0756                                                   */
/* C Fax: (217) 333-5962                                                     */
/* C Professor Jianming Jin                                                  */
/* C Department of Electrical and Computer Engineering                       */
/* C University of Illinois at Urbana-Champaign                              */
/* C 461 William L Everitt Laboratory                                        */
/* C 1406 West Green Street                                                  */
/* C Urbana, IL 61801-2991                                                   */
/*                                                                           */
/* Translated into C code by ????                                            */
/*                                                                           */
/* Modified by Josef Leydold on Tue Nov  1 13:22:09 CET 2011                 */
/*                                                                           */
/*---------------------------------------------------------------------------*/
{
  double x0,q1,q2,th,th1,th2,gr1,gi1;
  double t,x1,y1,sr,si;
  int j,k,na;
  double logq1;

  static const double a[] = {
    8.333333333333333e-02, -2.777777777777778e-03,
    7.936507936507937e-04, -5.952380952380952e-04,
    8.417508417508418e-04, -1.917526917526918e-03,
    6.410256410256410e-03, -2.955065359477124e-02,
    1.796443723688307e-01, -1.39243221690590 };

  /* real part and imaginary part of Gamma(z) */
  double gr, gi;

  y1 = 0.0;
  x1 = 0.0;
  na = 0;

  if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
    return UNUR_INFINITY;
  else if (x < 0.0) {
    x1 = x;
    y1 = y;
    x = -x;
    y = -y;
  }

  x0 = x;
  if (x <= 7.0) {
    na = (int)(7.0-x);
    x0 = x+na;
  }
  
  q1 = hypot(x0,y);    /* = sqrt(x0*x0+y*y) */
  th = atan(y/x0);
  logq1 = log(q1);
  gr = (x0-0.5)*logq1-th*y-x0+0.5*log(2.0*M_PI);
  gi = th*(x0-0.5)+y*logq1-y;

  for (k=0; k<10; k++){
    t = pow(q1,-1.0-2.0*k);
    gr += (a[k]*t*cos((2.0*k+1.0)*th));
    gi -= (a[k]*t*sin((2.0*k+1.0)*th));
  }

  if (x <= 7.0) {
    gr1 = 0.0;
    gi1 = 0.0;
    for (j=0; j<na; j++) {
      gr1 += (0.5*log((x+j)*(x+j)+y*y));
      gi1 += atan(y/(x+j));
    }
    gr -= gr1;
    gi -= gi1;
  }
  
  if (x1 < 0.0) {
    q1 = hypot(x,y);    /* = sqrt(x*x+y*y) */
    th1 = atan(y/x);
    sr = -sin(M_PI*x)*cosh(M_PI*y);
    si = -cos(M_PI*x)*sinh(M_PI*y);
    q2 = hypot(sr,si);  /* = sqrt(sr*sr+si*si) */
    th2 = atan(si/sr);
    if (sr < 0.0) th2 += M_PI;
    gr = log(M_PI/(q1*q2))-gr;
    gi = -th1-th2-gi;
    x = x1;
    y = y1;
  }

  /* For the case that we once need Gamma(z): */
  /* if (!islog) {                            */
  /*   double g0 = exp(gr);                   */
  /*   gr = g0 * cos(gi);                     */
  /*   gi = g0 * sin(gi);                     */
  /* }					      */

  /* return result */
  return gr;
  
} /* end of _unur_clgamma() */

/*---------------------------------------------------------------------------*/
