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
  double  u_resolution;      /* maximal error in u                           */
  double  bleft;             /* left border of the domain                    */
  double  bright;            /* right border of the domain                   */
};

/*---------------------------------------------------------------------------*/
/* store information about splines                                           */

#define UNUR_SINV_SPLINE_ORDER   (3)

struct unur_sinv_interval {

  double spline[UNUR_SINV_SPLINE_ORDER+1];   /* coefficients of spline       */
  double  bl, br;                   /* boundary of spline interval           */

  struct unur_sinv_interval *next;  /* pointer to next element in list       */

#ifdef UNUR_COOKIES
  unsigned cookie;      /* magic cookie                                      */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_sinv_gen { 
  double  u_resolution;      /* maximal error in u                           */
  double  bleft;             /* left border of the domain                    */
  double  bright;            /* right border of the domain                   */

  struct unur_sinv_interval *splines; /* pointer to splines                  */

  double  Umin, Umax;        /* bounds for iid random variable in respect to
                                the given (truncated) domain of the distr.   */
  double  CDFmin, CDFmax;    /* CDF-bounds of domain                         */
};

/*---------------------------------------------------------------------------*/
























