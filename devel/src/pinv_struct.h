/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: pinv_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method PINV                               *
 *         (Polynomial interpolation based INVersion of CDF)                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

struct unur_pinv_par { 
  int order;               /* order of interpolating polynomial              */
  double u_resolution;     /* maximal error in u                             */
  double bleft;            /* left boundary of the computational domain      */
  double bright;           /* right boundary of the computational domain     */
  int sleft;               /* whether to search for left boundary point      */
  int sright;              /* whether to search for right boundary point     */
/*   const double *stp;       /\* pointer to array of starting points            *\/ */
/*   int     n_stp;           /\* number of construction points at start         *\/ */
};

/*---------------------------------------------------------------------------*/
/* store information about splines                                           */


struct siv{
  double *ui;  //[g+1];    
  double *zi;  //[g+1];
  double xi;
  double cdfi;
/* #ifdef UNUR_COOKIES */
/*   unsigned cookie;         /\* magic cookie                                    *\/ */
/* #endif */
};

/* #define UNUR_PINV_MAX_ORDER   (5) */

/* struct unur_pinv_interval { */
/*   double spline[UNUR_PINV_MAX_ORDER+1];   /\* coefficients of spline           *\/ */
/*   double p;                /\* left design point (node) in interval            *\/   */
/*   double u;                /\* CDF at node p (u=CDF(p))                        *\/ */
/*   double f;                /\* PDF at node p (u=CDF(p))                        *\/ */
/*   double df;               /\* derivative of PDF at node p (u=CDF(p))          *\/ */

/*   struct unur_pinv_interval *next;  /\* pointer to next element in list        *\/ */

/* #ifdef UNUR_COOKIES */
/*   unsigned cookie;         /\* magic cookie                                    *\/ */
/* #endif */
/* }; */

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_pinv_gen { 
  int order;               /* order of interpolating polynomial              */

/*   int N;                   /\* total number of division points = #intervals+1 *\/ */
/*   double *intervals;       /\* pointer to array for storing data for intervals */
/* 			      in blocks of size order+2: */
/* 			      [0] ... u_{i-1} = CDF at left design point */
/* 			      [1] ... p_{i-1} = left design point = spline[0] */
/* 			      [2]-[order+1] ... spline[1] - spline[order]  */
/* 			      size of the array = N * (2+order)              *\/ */

  int    *guide;           /* pointer to guide table                         */ 
  int     guide_size;      /* size of guide table                            */

  double  Umax;            /* upper bound for uniform random variable U 
			      [ Umin = 0. ]                                  */

/*   double  CDFmin, CDFmax;  /\* CDF-bounds of domain                           *\/ */

  double  u_resolution;    /* maximal error in u                             */
  double  bleft;           /* left border of the computational domain        */
  double  bright;          /* right border of the computational domain       */

/*   struct unur_pinv_interval *iv; /\* linked list of splines (only used in setup) *\/ */
/*   double  tailcutoff_left; /\* cut point for left hand tail (u-value)         *\/  */
/*   double  tailcutoff_right;/\* cut point for right hand tail (u-value)        *\/  */
/*   const double *stp;       /\* pointer to array of starting points            *\/ */
/*   int     n_stp;           /\* number of construction points at start         *\/ */

  double  bleft_par;       /* border of the computational domain as ...      */
  double  bright_par;      /* ... given by user                              */

  int sleft;               /* whether to search for left boundary point      */
  int sright;              /* whether to search for right boundary point     */




  /* from pinvwh: */

  struct siv *iv;//[maxint+1] for setup; for sampling [ni+1]
  int ni;//number of sub intervals

};

/*---------------------------------------------------------------------------*/
