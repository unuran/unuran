/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_debug.h                                                    *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and function prototypes for debugging routines.    *
 *                                                                           *
 *   USAGE:                                                                  *
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
/* set generator id                                                          */

#ifdef UNUR_ENABLE_GENID
char *_unur_make_genid( const char *gentype );
#define _unur_set_genid(gentype) _unur_make_genid(gentype)
#define _unur_free_genid(gen)    free((gen)->genid)
#else
#define _unur_set_genid(gentype) (gentype)
#define _unur_free_genid(gen)    do { } while(0)
#endif

/*---------------------------------------------------------------------------*/
/* default debugging flag for generator                                      */

extern unsigned _unur_default_debugflag;     /* default debugging flags      */

/*---------------------------------------------------------------------------*/
/* warnings and error messages                                               */

/* Function prototypes                                                       */
void _unur_stream_printf( const char *genid, char *filename, int line, const char *format, ... );

/* an abbreviation */
#define _unur_print_if_default(par,flag)   if(!((par)->set & (flag))) fprintf(log,"  [default]")

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
/* Check for NULL pointer                                                    */

#ifdef UNUR_ENABLE_CHECKNULL

#define CHECK_NULL(ptr,rval)             \
  if (!(ptr)) {                          \
    _unur_error(NULL,UNUR_ERR_NULL,"");  \
    return rval;                         \
  }

#else               /* do not check (be carefull) */

#define CHECK_NULL(ptr,rval)  do {} while(0)

#endif

/* the second macro cannot be switched off by a compiler switch */
#define _unur_check_NULL(gid,ptr,rval)    \
  if (!(ptr)) {                           \
    _unur_error((gid),UNUR_ERR_NULL,"");  \
    return rval;                          \
  }

/*---------------------------------------------------------------------------*/
