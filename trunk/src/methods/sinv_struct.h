/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: sinv_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method SINV                               *
 *         (Spline approximation for INVerse of CDF)                         *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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
/* Information for constructing the generator                                */

struct unur_sinv_par { 
  int order;               /* order of interpolation polynomial              */
  double u_resolution;     /* maximal error in u                             */
  double  guide_factor;    /* relative size of guide table                   */
  double  bleft;           /* left border of the computational domain        */
  double  bright;          /* right border of the computational domain       */
};

/*---------------------------------------------------------------------------*/
/* store information about splines                                           */

#define UNUR_SINV_MAX_ORDER   (5)

struct unur_sinv_interval {
  double spline[UNUR_SINV_MAX_ORDER+1];   /* coefficients of spline           */
  double p;                /* left design point (node) in interval            */  
  double u;                /* CDF at node p (u=CDF(p))                        */
  double f;                /* PDF at node p (u=CDF(p))                        */
  double df;               /* derivative of PDF at node p (u=CDF(p))          */

  struct unur_sinv_interval *next;  /* pointer to next element in list        */

#ifdef UNUR_COOKIES
  unsigned cookie;         /* magic cookie                                    */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_sinv_gen { 
  int order;               /* order of interpolation polynomial              */

  int N;                   /* total number of division points = #intervals+1 */
  double *intervals;       /* pointer to array for storing data for intervals 
			      in blocks of size order+2:
			      [0] ... u_{i-1} = CDF at left design point
			      [1] ... p_{i-1} = left design point = spline[0]
			      [2]-[order+1] ... spline[1] - spline[order] 
			      size of the array = N * (2+order)              */

  int    *guide;           /* pointer to guide table                         */ 
  int     guide_size;      /* size of guide table                            */
  double  guide_factor;    /* relative size of guide table                   */

  double  Umin, Umax;      /* bounds for iid random variable in respect to
			      the given (truncated) domain of the distr.     */
  double  CDFmin, CDFmax;  /* CDF-bounds of domain                           */
  double  u_resolution;    /* maximal error in u                             */
  double  bleft;           /* left border of the computational domain        */
  double  bright;          /* right border of the computational domain       */

  struct unur_sinv_interval *iv; /* linked list of splines (only used in setup) */
};

/*---------------------------------------------------------------------------*/
























