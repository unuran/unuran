/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: fish.c                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         uniform random number generator provided by UNURAN                *
 *         random number generators inside UNURAN.                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *         Linear congruential generator                                     *
 *         x_(n+1) = 16807 * x_n mod (2^32 - 1)  ["Minimal Standard"]        *
 *                                                                           *
 *   WARNING:                                                                *
 *         Not state-of-the-art. SHOULD NOT BE USED ANY MORE.                *
 *         In UNURAN only as auxilliary second stream.                       *
 *         Should be replaced in future releases.                            *
 *                                                                           *
 *   REFERENCE:                                                              *
 *   Park, S. K. and Miller, K. W. (1988):                                   *
 *      Random number generators: good ones are hard to find,                *
 *      Comm. ACM 31, 1192--1201.                                            *
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
unur_urng_mstd (void)
{

# define a 16807       /* multiplicator */
# define m 2147483647  /* modulus */
# define q 127773      /* m / a */
# define r 2836        /* m % a */

  /* seed (must not be 0!) */
  static unsigned long x = 1804289383;

  int hi, lo, test;   /* intermediate results */

  hi = x / q;
  lo = x % q;
  test = a * lo - r * hi;
  x = (test > 0 ) ? test : test + m;
  return (x * 4.656612875245796924105750827e-10);

} /* end of unur_urng_mstd() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
