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
 *      Generators with Multipliers of the Form a = �2^q�2^r.                *
 *      in: J.A. Jones, R.R. Barton, K. Kang, and P.A. Fishwick (eds.),      *
 *      Proc. 2000 Winter Simulation Conference, 683-689.                    *  
 *                                                                           *
 *   Copyright for generator code by Renee Touzin.                           * 
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
#include <unur_source.h>
#include "unur_uniform.h"
/*---------------------------------------------------------------------------*/

/* seed (must not be 0!) */
#define SEED10  (12345L)
#define SEED11  (23456L)
#define SEED12  (34067L)
#define SEED20  (45678L)
#define SEED21  (56789L)
#define SEED22  (67890L)

/* status variable */
static unsigned long x10 = SEED10;
static unsigned long x11 = SEED11;
static unsigned long x12 = SEED12;
static unsigned long x20 = SEED20;
static unsigned long x21 = SEED21;
static unsigned long x22 = SEED22;

/* seed of last stream */
static unsigned long x10_start = SEED10;
static unsigned long x11_start = SEED11;
static unsigned long x12_start = SEED12;
static unsigned long x20_start = SEED20;
static unsigned long x21_start = SEED21;
static unsigned long x22_start = SEED22;

/*---------------------------------------------------------------------------*/

double
unur_urng_MRG31k3p (void)
     /* Combined multiple recursive generator.                               */
     /* Copyright (c) 2002 Renee Touzin.                                     */
{

# define m1      2147483647
# define m2      2147462579
# define norm    4.656612873077393e-10
# define mask11  511
# define mask12  16777215
# define mask20  65535
 

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

int
unur_urng_MRG31k3p_seed (long seed)
{
  if (seed==0) {
    _unur_error("URNG.fish",UNUR_ERR_GENERIC,"seed = 0");
    return 0;
  }
  
  /* the following is not really optimal */
  x10 = x10_start = seed; 
  x11 = x11_start = seed; 
  x12 = x12_start = seed; 
  x20 = x20_start = seed; 
  x21 = x21_start = seed; 
  x22 = x22_start = seed; 

  return 1;
} /* end of unur_urng_MRG31k3p_seed() */

/*---------------------------------------------------------------------------*/

int 
unur_urng_MRG31k3p_reset (void)
{
  x10 = x10_start;
  x11 = x11_start;
  x12 = x12_start;
  x20 = x20_start;
  x21 = x21_start;
  x22 = x22_start;

  return 1;
} /* end of unur_urng_MRG31k3p_reset() */

/*---------------------------------------------------------------------------*/
