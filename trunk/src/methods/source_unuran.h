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
/* include main header file                                                  */
#include <unuran.h>


#include <source_cookies.h>
#include <source_distr.h>
#include <source_errno.h>
#include <source_math.h>
#include <source_methods.h>

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
char *_unur_make_genid( const char *gentype );

/*---------------------------------------------------------------------------*/
/* Macros                                                                    */

#define min(x,y)   (((x)<(y)) ? (x) : (y))
#define max(x,y)   (((x)>(y)) ? (x) : (y))


/*---------------------------------------------------------------------------*/
#endif  /* end __UNURAN_SOURCE_H_SEEN */
/*---------------------------------------------------------------------------*/
