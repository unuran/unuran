/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_methods.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines bitmasks to identify used method in generator objects     *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_masks.h                                   *
 *                                                                           *
 *                                                                           *
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
#ifndef __SOURCE_METHODS_H_SEEN
#define __SOURCE_METHODS_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Bitmask to indicate methods                                            **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* bitmasks                                                                  */

#define UNUR_MASK_TYPE     0xff000000u   /* indicate type of method           */

/* discrete distributions */
#define UNUR_METH_DISCR    0x10000000u

#define UNUR_METH_DAU      0x10000001u
#define UNUR_METH_DIS      0x10000002u

/* continuous distributions */
#define UNUR_METH_CONT     0x20000000u

#define UNUR_METH_AROU     0x20000100u
#define UNUR_METH_SROU     0x20000200u
#define UNUR_METH_STDR     0x20000300u
#define UNUR_METH_TABL     0x20000400u
#define UNUR_METH_TDR      0x20000500u
#define UNUR_METH_UNIF     0x20000600u
#define UNUR_METH_UTDR     0x20000700u

/* multivariate continuous distributions */
#define UNUR_METH_VEC      0x40000000u

#define UNUR_METH_RECT     0x40010000u

/* generators for standard distributions */
#define UNUR_METH_CSTD     0x2000ff00u   /* is of type UNUR_METH_CONT !! */

/* to indicate unkown type */
#define UNUR_METH_UNKNOWN  0xff000000u

/*****************************************************************************/
/**  Macros                                                                 **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* check if parameter object is of correct type, return 0 otherwise       */

#define _unur_check_par_object( type ) \
  if ( par->method != UNUR_METH_##type ) { \
    _unur_warning(#type,UNUR_ERR_PAR_INVALID,""); \
    return 0; } \
  COOKIE_CHECK(par,CK_##type##_PAR,0)

/*---------------------------------------------------------------------------*/
#endif  /* end __SOURCE_METHODS_H_SEEN */
/*---------------------------------------------------------------------------*/
