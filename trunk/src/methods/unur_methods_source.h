/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_methods_source.h                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines bitmasks to identify used method in generator objects     *
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
#ifndef UNUR_METHODS_SOURCE_H_SEEN
#define UNUR_METHODS_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Bitmask to indicate methods                                            **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* bitmasks                                                                  */

#define UNUR_MASK_TYPE     0xff000000u   /* indicate type of method           */

/* discrete distributions */
#define UNUR_METH_DISCR    0x10000000u

#define UNUR_METH_DARI     0x10000001u
#define UNUR_METH_DAU      0x10000002u
#define UNUR_METH_DGT      0x10000003u
#define UNUR_METH_DSROU    0x10000004u
#define UNUR_METH_DSS      0x10000008u

/* continuous distributions */
#define UNUR_METH_CONT     0x20000000u

#define UNUR_METH_AROU     0x20000100u
#define UNUR_METH_HINV     0x20000200u
#define UNUR_METH_HRB      0x20000a00u
#define UNUR_METH_HRD      0x20000b00u
#define UNUR_METH_HRI      0x20000c00u
#define UNUR_METH_NINV     0x20000300u
#define UNUR_METH_SROU     0x20000400u
#define UNUR_METH_SSR      0x20000500u
#define UNUR_METH_TABL     0x20000600u
#define UNUR_METH_TDR      0x20000700u
#define UNUR_METH_UNIF     0x20000800u
#define UNUR_METH_UTDR     0x20000900u

/* univariate continuous empirical distributions */
#define UNUR_METH_CEMP     0x40000000u

#define UNUR_METH_EMPK     0x40001100u
#define UNUR_METH_EMPL     0x40001200u

/* multivariate continuous distributions */
#define UNUR_METH_VEC      0x80000000u

#define UNUR_METH_VMT      0x80010000u
#define UNUR_METH_VEMPK    0x80020000u

/* generators for standard distributions */
#define UNUR_METH_CSTD     0x2000f100u   /* is of type UNUR_METH_CONT !!  */
#define UNUR_METH_DSTD     0x1000f200u   /* is of type UNUR_METH_DISCR !! */

/* automatically selected generator */
#define UNUR_METH_AUTO     0x00a00000u   /* can be any type of distribution */

/* to indicate unkown type */
#define UNUR_METH_UNKNOWN  0xff000000u

/*****************************************************************************/
/**  Macros                                                                 **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* check if parameter object is of correct type, return 0 otherwise       */

#define _unur_check_par_object( par,type ) \
  if ( (par)->method != UNUR_METH_##type ) { \
    _unur_warning(#type,UNUR_ERR_PAR_INVALID,""); \
    return (UNUR_ERR_PAR_INVALID); } \
  COOKIE_CHECK((par),CK_##type##_PAR,UNUR_ERR_COOKIE)

/*---------------------------------------------------------------------------*/
/* check if generator object is of correct type, return 0 otherwise          */

#define _unur_check_gen_object( gen,type,rval ) \
  if ( (gen)->method != UNUR_METH_##type ) { \
    _unur_warning((gen)->genid,UNUR_ERR_GEN_INVALID,""); \
    return rval; } \
  COOKIE_CHECK((gen),CK_##type##_GEN,UNUR_ERR_COOKIE)

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_METHODS_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
