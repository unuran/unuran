/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_misc.h                                                       *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for miscelleanous          *
 *         mathematical routines                                             *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
#ifndef X_MATH_H_SEEN
#define X_MATH_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/* 
   =NODE Math Mathematics

   =UP Misc [0]

   =DESCRIPTION
      The following macros have been defined

      @ftable @code
      @item UNUR_INFINITY
      indicates infinity for floating point numbers (of type @code{double}).
      Internally @code{HUGE_VAL} is used.

      @item INT_MAX
      @itemx INT_MIN
      indicate infinity and minus infinity, resp., for integers
      (defined by ISO C standard).

      @item  TRUE
      @itemx FALSE
      bolean expression for return values of @code{set} functions.
      @end ftable

   =END
*/

/*---------------------------------------------------------------------------*/
/* Defining infinity                                                         */
/* (we use the largest possible value to indicate infinity)                  */
#include <math.h>

#define UNUR_INFINITY      HUGE_VAL     

/*---------------------------------------------------------------------------*/
/* True and false                                                            */

#ifndef TRUE
#define TRUE   (1)
#endif

#ifndef FALSE
#define FALSE  (0)
#endif

/*---------------------------------------------------------------------------*/
#endif  /* X_MATH_H_SEEN */
/*---------------------------------------------------------------------------*/






