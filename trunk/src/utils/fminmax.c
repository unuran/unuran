/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      fminmax.c                                                    *
 *                                                                           *
 *   Find maximum or minimum of a function using Brent's algorithm.          *
 *                                                                           *
 *****************************************************************************
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

#include <unur_source.h>
#include "fminmax_source.h"

static double _unur_function_minmax(struct UNUR_FUNCT_GENERIC fs, int minmax, 
                             double a, double b, double c, double tol);  
/* calculation routine for minimum or maximum of a continous function */

/*---------------------------------------------------------------------------*/


/* Wrapper function for the max calculation routine */
double
_unur_function_max(fs, a, b, c, tol)  /* An estimate to the max location     */
     struct UNUR_FUNCT_GENERIC fs;    /* Function struct                     */
     double a;                          /* Left border | of the range	     */
     double b;                          /* Right border| the min is seeked   */
     double c;                          /* first guess for the max           */
     double tol;                        /* Acceptable tolerance              */
{
   return _unur_function_minmax(fs, -1, a, b, c, tol);
} /* end of _unur_function_max() */

/*---------------------------------------------------------------------------*/

/* Wrapper function for the min calculation routine */
double
_unur_function_min(fs, a, b, c, tol)  /* An estimate to the min location   */
     struct UNUR_FUNCT_GENERIC fs;    /* Function struct                     */
     double a;                          /* Left border | of the range	     */
     double b;                          /* Right border| the min is seeked   */
     double c;                          /* first guess for the min           */
     double tol;                        /* Acceptable tolerance              */
{
   return _unur_function_minmax(fs, +1,  a, b, c, tol);
} /* end of _unur_function_min() */

/*---------------------------------------------------------------------------*/

/*
 *****************************************************************************
 *	    		    C math library                                   *
 * function FMINBR - one-dimensional search for a function minimum           *
 *			  over the given range                               *
 *                                                                           *
 * Author: Oleg Keselyov.                                                    *
 *                                                                           *
 * modified by Josef Leydold (documentation unchanged)                       *
 *                                                                           *
 * Input                                                                     *
 *	double fminbr(f, a,b,c,tol)                                          *
 *	double a; 			Minimum will be seeked for over      *
 *	double b;  			a range [a,b], a being < b.          *
 *      double c;                       c within (a,b) is first guess        *
 *	double (*f)(double x);		Name of the function whose minimum   *
 *					will be seeked for                   *
 *	double tol;			Acceptable tolerance for the minimum *
 *					location. It have to be positive     *
 *					(e.g. may be specified as EPSILON)   *
 *                                                                           *
 * Output                                                                    *
 *	Fminbr returns an estimate for the minimum location with accuracy    *
 *	3*SQRT_EPSILON*abs(x) + tol.                                         *
 *	The function always obtains a local minimum which coincides with     *
 *	the global one only if a function under investigation being          *
 *	unimodular.                                                          *
 *	If a function being examined possesses no local minimum within       *
 *	the given range, Fminbr returns 'a' (if f(a) < f(b)), otherwise      *
 *	it returns the right range boundary value b.                         *
 *                                                                           *
 * Algorithm                                                                 *
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical    *
 *	computations. M., Mir, 1980, p.202 of the Russian edition            *
 *                                                                           *
 *	The function makes use of the "gold section" procedure combined with *
 *	the parabolic interpolation.                                         *
 *	At every step program operates three abscissae - x,v, and w.         *
 *	x - the last and the best approximation to the minimum location,     *
 *	    i.e. f(x) <= f(a) or/and f(x) <= f(b)                            *
 *	    (if the function f has a local minimum in (a,b), then the both   *
 *	    conditions are fulfiled after one or two steps).                 *
 *	v,w are previous approximations to the minimum location. They may    *
 *	coincide with a, b, or x (although the algorithm tries to make all   *
 *	u, v, and w distinct). Points x, v, and w are used to construct      *
 *	interpolating parabola whose minimum will be treated as a new        *
 *	approximation to the minimum location if the former falls within     *
 *	[a,b] and reduces the range enveloping minimum more efficient than   *
 *	the gold section procedure.                                          *
 *	When f(x) has a second derivative positive at the minimum location   *
 *	(not coinciding with a or b) the procedure converges superlinearly   *
 *	at a rate order about 1.324                                          *
 *                                                                           *
 *****************************************************************************/

/* in case of any error INFINITY is returned */

double
_unur_function_minmax(fs, minmax, a, b, c, tol)  
                                   /* An estimate to the min or max location */
     struct UNUR_FUNCT_GENERIC fs; /* Function struct                        */
     int minmax ; 	           /* -1 for maximum, +1 for minimum         */
     double a;                     /* Left border | of the range	     */
     double b;                     /* Right border| the min is seeked        */
     double c;                     /* first guess for the min/max            */
     double tol;                   /* Acceptable tolerance                   */
{
#define f(x) ( (minmax) * ((fs.f)(x, fs.params)) )          
#define SQRT_EPSILON  (1.e-7)           /* tolerance for relative error      */
#define MAXIT         (1000)            /* maximum number of iterations      */


  int i;

  double x,v,w;                         /* Abscissae, descr. see above       */
  double fx;                            /* f(x)                              */
  double fv;                            /* f(v)                              */
  double fw;                            /* f(w)                              */
  const double r = (3.-sqrt(5.0))/2;    /* Gold section ratio                */

  /* check arguments */
  CHECK_NULL(fs.f,INFINITY);
  if ( tol < 0. || b <= a || c <= a || b <= c) {
    _unur_error("CMAX",UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;
  }


  /* Origially the third point was computed by golden section. In the        */
  /* modified version it is given as point `c' by the calling function.      */
  /* v = a + r*(b-a);  fv = f(v);        First step - always gold section    */

  v = c;  fv = f(v);                    /* First step */
  x = v;  w = v;
  fx=fv;  fw=fv;

  for(i=0; i < MAXIT; i++) {            /* Main iteration loop	*/
    double range = b-a;                 /* Range over which the minimum is   */
                                        /* seeked for		             */
    double middle_range = (a+b)/2;
    double tol_act =                    /* Actual tolerance                  */
		SQRT_EPSILON*fabs(x) + tol/3;
    double new_step;                    /* Step at this iteration            */
       
    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;                         /* Acceptable approx. is found       */

                                        /* Obtain the gold section step      */
    new_step = r * ( x<middle_range ? b-x : a-x );


                                 /* Decide if the interpolation can be tried */
    if( fabs(x-w) >= tol_act ) {        /* If x and w are distinct           */
                                        /* interpolatiom may be tried	     */
      register double p;                /* Interpolation step is calculated  */
      register double q;                /* as p/q; division operation        */
                                        /* is delayed until last moment      */
      register double t;

      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);

      if( q>(double)0 )                 /* q was calculated with the         */
	p = -p;                         /* opposite sign; make q positive    */
      else                              /* and assign possible minus to	p    */
	q = -q;

      if( fabs(p) < fabs(new_step*q) && /* If x + p/q falls in [a,b] not too */
	  p > q*(a-x+2*tol_act)      && /* close to a and b, and isn't too   */
	  p < q*(b-x-2*tol_act)  )      /* large, it is accepted.            */
	new_step = p/q;                 /* If p/q is too large then the	     */
                                        /* gold section procedure can reduce */
                                        /* [a,b] range to more extent.       */
    }

    if( fabs(new_step) < tol_act ) {    /* Adjust the step to be not less    */
      if( new_step > (double)0 )        /* than tolerance.                   */
	new_step = tol_act;
      else
	new_step = -tol_act;
    }
                                 /* Obtain the next approximation to min     */
    {                            /* and reduce the enveloping range.         */
      register double t = x + new_step;	/* Tentative point for the min       */
      register double ft = f(t);

      if( ft <= fx ) {                  /* t is a better approximation	     */
	if( t < x )   
	  b = x;                        /* Reduce the range so that          */
	else                            /* t would fall within it            */
	  a = x;
      
	v = w;  w = x;  x = t;          /* Assign the best approx to x	     */
	fv=fw;  fw=fx;  fx=ft;
      }
      else {                            /* x remains the better approx       */
	if( t < x )
	  a = t;                        /* Reduce the range enclosing x      */
	else
	  b = t;
      
        if( ft <= fw || w==x ) {
	  v = w;  w = t;
	  fv=fw;  fw=ft;
        }
        else if( ft<=fv || v==x || v==w ) {
	  v = t;
	  fv=ft;
        }
      }
    }                   /* ----- end-of-block ----- */
  }                /* ===== End of for loop ===== */

  /* maximal number of interations exceeded */
  return INFINITY;

#undef f
#undef MAXIT
#undef SQRT_EPSILON
} /* end of _unur_function_min() */


/*---------------------------------------------------------------------------*/
