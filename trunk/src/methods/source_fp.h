/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_fp.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares macros for floating point arithmetic.                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file.                                             *
 *         only included in source_unuran.h                                  *
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
#ifndef __SOURCE_FP_H_SEEN
#define __SOURCE_FP_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Comparisons                                                               */

/* a == b (except precision bit) */
#define _unur_FP_same(a,b) \
 ((a)==(b) || \
 fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b)) * DBL_EPSILON)

/* a == b */
#define _unur_FP_equal(a,b) \
 ((a)==(b) || \
 fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b)) * UNUR_EPSILON)

/* a is approximately equal to b */
#define _unur_FP_approx(a,b) \
 ((a)==(b) || \
 fabs((a)-(b)) <= ((fabs(a)<fabs(b))?fabs(a):fabs(b)) * FLT_EPSILON)

/* a < b */
#define _unur_FP_less(a,b) \
 (!_unur_FP_equal((a),(b)) && ((a) < (b)))

/* a > b */
#define _unur_FP_greater(a,b) \
 (!_unur_FP_equal((a),(b)) && ((a) > (b)))

/*---------------------------------------------------------------------------*/
/* Infinity                                                                  */

/* Defining infinity (just to avoid writing UNUR_INIFINITY)                  */

#define INFINITY  UNUR_INFINITY  /* This must be already defined in x_math.h */

/* check for infinity */

/* +oo */
#define _unur_FP_is_infinity(a)  ((a) >= INFINITY)

/* -oo */
#define _unur_FP_is_minus_infinity(a)  ((a) <= -INFINITY)


/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_FP_H_SEEN */
/*---------------------------------------------------------------------------*/






