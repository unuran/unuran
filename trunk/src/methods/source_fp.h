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

/* a == b */
#define _FP_equal(a,b) \
 (fabs((a)-(b)) <= ((fabs(a)<fabs(b)) ? fabs(a) : fabs(b))*UNUR_EPSILON)

/* a < b */
#define _FP_less(a,b) \
 (!_FP_equal((a),(b)) && ((a) < (b)))

/* a > b */
#define _FP_greater(a,b) \
 (!_FP_equal((a),(b)) && ((a) > (b)))

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_FP_H_SEEN */
/*---------------------------------------------------------------------------*/






