/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mrg31k3p.c                                                        *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         uniform random number generator provided by UNURAN                *
 *         random number generators inside UNURAN.                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *         Combined multiple recursive generator.                            *
 *         The two components of the generator are                           *
 *                                                                           *
 *   x_{1,k}  = (2^{22} x_{1,k-2} + (2^7 +1)x_{1,k-3}) mod (2^{31}-1)        *
 *   x_{2,k}  = (2^{15} x_{2,k-1} + (2^{15} +1)x_{2,k-3} mod (2^{31}-21069)  *
 *                                                                           *
 *         x_{1,k} and x_{2,k} are combined together and the result is       *
 *         multiplied by 1/(2^{31}-1) to have a number between 0 and 1.      *
 *                                                                           *
 *   REFERENCE:                                                              *
 *   L'Ecuyer, P. and R. Touzin (2000): Fast Combined Multiple Recursive     *
 *      Generators with Multipliers of the Form a = ±2^q±2^r.                *
 *      in: J.A. Jones, R.R. Barton, K. Kang, and P.A. Fishwick (eds.),      *
 *      Proc. 2000 Winter Simulation Conference, 683-689.                    *  
 *                                                                           *
 *   Copyright for the original code by Renee Touzin.                        * 
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
#if UNUR_URNG_TYPE == UNUR_URNG_SIMPLE || UNUR_URNG_TYPE == UNUR_URNG_GENERIC
/*---------------------------------------------------------------------------*/

double 
unur_urng_MRG31k3p (void)
{

# define m1      2147483647
# define m2      2147462579
# define norm    4.656612873077393e-10
# define mask11  511
# define mask12  16777215
# define mask20  65535
 
  /* seed */
  static unsigned long x10 = 12345;
  static unsigned long x11 = 23456;
  static unsigned long x12 = 34567;
  static unsigned long x20 = 45678;
  static unsigned long x21 = 56789;
  static unsigned long x22 = 67890;

  register unsigned long y1, y2;  /* For intermediate results */
  
  /* First component */
  y1 = ( (((x11 & mask11) << 22) + (x11 >> 9))
	 + (((x12 & mask12) << 7)  + (x12 >> 24)) );
  if (y1 > m1) y1 -= m1;
  y1 += x12;
  if (y1 > m1) y1 -= m1;
  x12 = x11;  x11 = x10;  x10 = y1;
 
  /* Second component */
  y1 = ((x20 & mask20) << 15) + 21069 * (x20 >> 16);
  if (y1 > m2) y1 -= m2;
  y2 = ((x22 & mask20) << 15) + 21069 * (x22 >> 16);
  if (y2 > m2) y2 -= m2;
  y2 += x22;
  if (y2 > m2) y2 -= m2;
  y2 += y1;
  if (y2 > m2) y2 -= m2;
  x22 = x21;  x21 = x20;  x20 = y2;

  /* Combination */
  if (x10 <= x20)
    return ((x10 - x20 + m1) * norm);
  else 
    return ((x10 - x20) * norm);

} /* end of unur_urng_MRG31k3p() */
 
/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

