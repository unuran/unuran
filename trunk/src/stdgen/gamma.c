/**
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
