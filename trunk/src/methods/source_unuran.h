/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: source_unuran.h                                                   *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         defines macros and declares structures and function prototypes    *
 *         for all UNURAN source files                                       *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source files.                                    *
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
#ifndef __UNURAN_SOURCE_H_SEEN
#define __UNURAN_SOURCE_H_SEEN
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* include main header files                                                  */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <config.h>
#include <unuran.h>

#include <unuran_errno.h>
#include <source_cookies.h>
#include <source_distr.h>
#include <source_math.h>
#include <source_methods.h>

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
/* Function prototypes for allocating memory blocks                          */
void *_unur_malloc(size_t size);
void  _unur_add_mblocks( struct unur_mblock **mblocks, void *ptr );
void  _unur_free_mblocks( struct unur_mblock *mblocks );

/*---------------------------------------------------------------------------*/

#define TOLERANCE 1.e-10  /* Provisorium !!!! */   /** TODO !! **/

/*---------------------------------------------------------------------------*/
/* True and false                                                            */

#ifndef TRUE
#define TRUE   (1)
#endif

#ifndef FALSE
#define FALSE  (0)
#endif

/*---------------------------------------------------------------------------*/
/* warnings and error messages                                               */

/* Function prototypes                                                       */
void _unur_stream_printf( const char *genid, char *filename, int line, const char *format, ... );

extern unsigned unur_errno;  /* global variable used to record errors        */

/*---------------------------------------------------------------------------*/
#ifdef UNUR_WARNINGS_ON    /* warnings enabled */
/*---------------------------------------------------------------------------*/

#define _unur_error(genid,errortype,str) \
   do { \
      unur_errno = (errortype); \
      _unur_stream_printf((genid),__FILE__,__LINE__,"error: %s %s", \
                          unur_get_strerror(errortype), (str) ); \
   } while (0)

#define _unur_warning(genid,errortype,str) \
   do { \
      unur_errno = (errortype); \
      _unur_stream_printf((genid),__FILE__,__LINE__,"warning: %s %s", \
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

#ifdef UNUR_CHECKNULL

#define CHECK_NULL(ptr,rval)                 \
  if (!(ptr)) {                              \
    _unur_error(NULL,UNUR_ERR_NULL,"");      \
    return rval;                             \
  }

#else               /* do not check (be carefull) */

#define CHECK_NULL(ptr,rval)  do {} while(0)

#endif

/*---------------------------------------------------------------------------*/
/* write infos into log file                                                 */

/* an abbreviation */
#define _unur_print_if_default(par,flag)   if(!((par)->set & (flag))) fprintf(log,"  [default]")


/*---------------------------------------------------------------------------*/

double _unur_arcmean( double x0, double x1 );

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

#define min(x,y)   (((x)<(y)) ? (x) : (y))
#define max(x,y)   (((x)>(y)) ? (x) : (y))


/*---------------------------------------------------------------------------*/
#endif  /* end __UNURAN_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
