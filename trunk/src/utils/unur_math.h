/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_math.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros, constants, structures, function prototypes, etc. *
 *         for using mathematics in UNURAN.                                  *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file.                                             *
 *         included in all source files that need mathematical functions.    *
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
#ifndef __UNUR_MATH_H_SEEN
#define __UNUR_MATH_H_SEEN
/*---------------------------------------------------------------------------*/

#include <math.h>
#include <unur_defs.h>

/*---------------------------------------------------------------------------*/
/* Defining infinity (just to avoid writing UNUR_INIFINITY)                  */

#define INFINITY        UNUR_INFINITY    /* (see unur_defs.h)                */

/*---------------------------------------------------------------------------*/
/* mathematical constants                                                    */

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328
#endif

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_MATH_H_SEEN */
/*---------------------------------------------------------------------------*/






