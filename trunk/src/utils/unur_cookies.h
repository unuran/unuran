/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_cookies.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines magic cookies.                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         internal header file.                                             *
 *         included in all source files that use magis cookies.              *
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
#ifndef __UNUR_COOKIES_H_SEEN
#define __UNUR_COOKIES_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>

/*---------------------------------------------------------------------------*/
#if UNUR_DEBUG & UNUR_DB_COOKIES    /* use magic cookies                     */
/*---------------------------------------------------------------------------*/
/* name of cookies                                                           */

#define CK_DIS_PAR       0x00000010UL
#define CK_DIS_GEN       0x00000011UL
#define CK_DAU_PAR       0x00000020UL
#define CK_DAU_GEN       0x00000021UL

#define CK_AROU_PAR      0x00100010UL
#define CK_AROU_GEN      0x00100011UL
#define CK_AROU_SEG      0x00100012UL
#define CK_TABL_PAR      0x00100020UL
#define CK_TABL_GEN      0x00100021UL
#define CK_TABL_IV       0x00100022UL
#define CK_TDR_PAR       0x00100030UL
#define CK_TDR_GEN       0x00100031UL
#define CK_TDR_IV        0x00100032UL
#define CK_UNIF_PAR      0x00100040UL
#define CK_UNIF_GEN      0x00100041UL
#define CK_UTDR_PAR      0x00100050UL
#define CK_UTDR_GEN      0x00100051UL

#define CK_RECT_PAR      0x00200010UL
#define CK_RECT_GEN      0x00200011UL

#define CK_CSTD_PAR      0x10000010UL
#define CK_CSTD_GEN      0x10000011UL

#define CK_MBLOCK        0xf0000001UL

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
#endif  /* __UNUR_COOKIES_H_SEEN */
/*---------------------------------------------------------------------------*/






