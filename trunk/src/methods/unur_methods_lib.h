/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_methods_lib.h                                                *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and declares structures and function prototypes    *
 *         for all UNURAN methods                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in all methods source files.                        *
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
#ifndef __UNUR_METHODS_LIB_H_SEEN
#define __UNUR_METHODS_LIB_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>

/*---------------------------------------------------------------------------*/
/* typedefs                                                                  */

/*---------------------------------------------------------------------------*/
/* List of methods                                                           */

#define UNUR_MASK_VARIANT  0x00000fffu   /* indicate variant (see the corresponding .c files) */
#define UNUR_MASK_METHOD   0xfff00000u   /* indicate method                   */
#define UNUR_MASK_TYPE     0xf0000000u   /* indicate type of method           */

/* bits 13-20 are used for flags common to all generators */
#define UNUR_MASK_SCHECK   0x00000001u   /* turns check sampling on/off       */


#define UNUR_MASK_MODE     0x00002000u   /* use mode                          */

/* discrete, univariate */
#define UNUR_METH_DISCR    0x10000000u

#define UNUR_METH_DAU      0x10100000u
#define UNUR_METH_DIS      0x10200000u

/* continuous, univariate */
#define UNUR_METH_CONT     0x20000000u

#define UNUR_METH_AROU     0x20300000u
#define UNUR_METH_TABL     0x20400000u
#define UNUR_METH_TDR      0x20500000u
#define UNUR_METH_UNIF     0x20600000u
#define UNUR_METH_UTDR     0x20700000u

/* continuous, multivariate */
#define UNUR_METH_VEC      0x40000000u

#define UNUR_METH_RECT     0x40700000u

/* generators for standard distributions                                     */
/* for definitions of methods for standard distributions see "stand.c"       */
#define UNUR_MASK_DISTR    0x000ffff0u   /* indicate distribution           */

#define UNUR_METH_CSTD     0x2f000000u   /* is of type UNUR_METH_CONT !! */

/* to indicate unkown type */
#define UNUR_METH_UNKNOWN  0xf0000000u


/*---------------------------------------------------------------------------*/
/* check object                                                              */

/* check if parameter object is of correct type, return 0 otherwise       */
#define _unur_check_par_object( type ) \
  if ( par->method != UNUR_METH_##type ) { \
    _unur_warning(GENTYPE,UNUR_ERR_PAR_INVALID,""); \
    return 0; } \
  COOKIE_CHECK(par,CK_##type##_PAR,0)

/*---------------------------------------------------------------------------*/
/* we cannot load the next files until all definitions are done              */

/*---------------------------------------------------------------------------*/
/* misc prototype                                                            */

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_METHODS_LIB_H_SEEN */
/*---------------------------------------------------------------------------*/


