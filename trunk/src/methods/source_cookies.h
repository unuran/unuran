/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_cookies.h                                                  *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines magic cookies.                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file                                              *
 *          included only in source_unuran.h                                 *
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
#ifndef __SOURCE_COOKIES_H_SEEN
#define __SOURCE_COOKIES_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unuran_config.h>

/*---------------------------------------------------------------------------*/
#ifdef UNUR_COOKIES    /* use magic cookies                                  */
/*---------------------------------------------------------------------------*/
/* name of cookies                                                           */

/* generators discrete distributions */
#define CK_DIS_PAR       0x00000010u
#define CK_DIS_GEN       0x00000011u
#define CK_DAU_PAR       0x00000020u
#define CK_DAU_GEN       0x00000021u

/* generators continuous distributions */
#define CK_AROU_PAR      0x00100010u
#define CK_AROU_GEN      0x00100011u
#define CK_AROU_SEG      0x00100012u
#define CK_SROU_PAR      0x00100060u
#define CK_SROU_GEN      0x00100061u
#define CK_STDR_PAR      0x00100070u
#define CK_STDR_GEN      0x00100071u
#define CK_TABL_PAR      0x00100020u
#define CK_TABL_GEN      0x00100021u
#define CK_TABL_IV       0x00100022u
#define CK_TDR_PAR       0x00100030u
#define CK_TDR_GEN       0x00100031u
#define CK_TDR_IV        0x00100032u
#define CK_UNIF_PAR      0x00100040u
#define CK_UNIF_GEN      0x00100041u
#define CK_UTDR_PAR      0x00100050u
#define CK_UTDR_GEN      0x00100051u

/* generators multivariate continuous distributions */
#define CK_RECT_PAR      0x00200010u
#define CK_RECT_GEN      0x00200011u

#define CK_CSTD_PAR      0x10000010u
#define CK_CSTD_GEN      0x10000011u

#define CK_MBLOCK        0xf0000001u

/* distribution objects */
#define CK_DISTR         0xe0000000u
#define CK_DISTR_CONT    0xe0000001u
#define CK_DISTR_DISCR   0xe0000002u

#define CK_SPECIALGEN_CONT 0xd00001u

/*---------------------------------------------------------------------------*/
/* macros for dealing with magic cookies                                     */

#define COOKIE_SET(ptr,ck)           (ptr)->cookie=(ck)

#define COOKIE_SET_ARRAY(ptr,ck,n)   {                                          \
  register int i;                                                               \
  for (i=0;i<(n);i++)                                                           \
    ((ptr)+(i))->cookie=(ck);                                                   \
}

#define COOKIE_CHECK(ptr,ck,rval)                                               \
  if((ptr)->cookie!=(ck)) {                                                     \
    _unur_stream_printf(NULL,__FILE__,__LINE__,                                 \
			"warning: %s (observed = %#lx, expected = %#lx)",       \
                        unur_get_strerror(UNUR_ERR_COOKIE),                     \
                        (ptr)->cookie, (ck));                                   \
    return rval;                                                                \
  }

/*---------------------------------------------------------------------------*/
#else                               /* do not use magic cookies              */
/*---------------------------------------------------------------------------*/

#define COOKIE_SET(ptr,ck) 
#define COOKIE_SET_ARRAY(ptr,ck,n)
#define COOKIE_CHECK(ptr,ck,rval) 

/*---------------------------------------------------------------------------*/
#endif 
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* __SOURCE_COOKIES_H_SEEN */
/*---------------------------------------------------------------------------*/






