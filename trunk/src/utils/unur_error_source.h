/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_error_source.h                                               *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes error messages             *
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
#ifndef UNUR_ERROR_SOURCE_H_SEEN
#define UNUR_ERROR_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_WARNINGS_ON    /* warnings enabled */
/*---------------------------------------------------------------------------*/
#define _unur_error(genid,errortype,str) \
   do { \
      unur_errno = (errortype); \
      _unur_stream_printf((genid),__FILE__,__LINE__,"error: %s: %s", \
                          unur_get_strerror(errortype), (str) ); \
   } while (0)

#define _unur_warning(genid,errortype,str) \
   do { \
      unur_errno = (errortype); \
      _unur_stream_printf((genid),__FILE__,__LINE__,"warning: %s: %s", \
                          unur_get_strerror(errortype), (str) ); \
   } while (0)
/*---------------------------------------------------------------------------*/
#else   /* warnings disabled */
/*---------------------------------------------------------------------------*/
#define _unur_error(genid,errortype,str)      do { unur_errno = (errortype); } while(0)
#define _unur_warning(genid,errortype,str)    do { unur_errno = (errortype); } while(0)
/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_WARNINGS_ON */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif  /* UNUR_ERROR_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
