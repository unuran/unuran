/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ptx_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method PTX                               *
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

struct unur_ptx_par { 
  int order;               /* order of interpolating polynomial              */
  double u_resolution;     /* maximal error in u                             */
  double bleft;            /* left boundary of the computational domain      */
  double bright;           /* right boundary of the computational domain     */
  int sleft;               /* whether to search for left boundary point      */
  int sright;              /* whether to search for right boundary point     */
  int max_ivs;             /* maximum number of subintervals                 */

  const struct unur_distr *din;
};

/*---------------------------------------------------------------------------*/
/* store information about splines                                           */

struct unur_ptx_interval {
  double *ti;  /* points for constructing Newton interpolation */
  double *zi;  /* values of inverse CDF at these T values */
  double xi;   /* left point of interval */
  double cdfi; /* CDF at left point of interval */

#ifdef UNUR_COOKIES
  unsigned cookie;         /* magic cookie                                    */
#endif
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_ptx_gen { 
  int order;               /* order of interpolating polynomial              */

/*   int    *guide;           /\* pointer to guide table                         *\/  */
/*   int     guide_size;      /\* size of guide table                            *\/ */


  double  Tmin;
  double  Tmax;

  double  Umax;            /* upper bound for uniform random variable U 
			      [ Umin = 0. ]                                  */

  double  u_resolution;    /* maximal error in u                             */
  double  bleft;           /* left border of the computational domain        */
  double  bright;          /* right border of the computational domain       */

  struct unur_ptx_interval *iv; /* list of intervals                        */
  int n_ivs;               /* number of subintervals                         */
  int max_ivs;             /* maximum number of subintervals                 */


  double  bleft_par;       /* border of the computational domain as ...      */
  double  bright_par;      /* ... given by user                              */

  double  dleft;           /* left and right boundary of domain / support    */
  double  dright;          /* of distribution                                */

  int sleft;               /* whether to search for left boundary point      */
  int sright;              /* whether to search for right boundary point     */

  double area;             /* approximate area below PDF                     */ 
  double logPDFconstant;   /* rescaling constant for logPDF                  */


  struct unur_lobatto_table *aCDF; /* polynomial approximation of CDF        */

  struct unur_lobatto_table *aTRX; /* polynomial approximation of FIXME      */


};

/*---------------------------------------------------------------------------*/
