/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tools.h                                                      *
 *                                                                           *
 *   auxilliary subroutines                                                  *
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
#ifndef __UNUR_UTILS_H_SEEN
#define __UNUR_UTILS_H_SEEN
/*---------------------------------------------------------------------------*/

#include <unur_defs.h>

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

#if UNUR_DEBUG & UNUR_DB_CHECKNULL

#define CHECK_NULL(ptr,rval)                 \
  if (!(ptr)) {                              \
    _unur_error(NULL,UNUR_ERR_NULL,"");      \
    return rval;                             \
  }

#else               /* do not check (be carefull) */

#define CHECK_NULL(ptr,rval)                 \
        do {} while(0)

#endif

/*---------------------------------------------------------------------------*/

double _unur_arcmean( double x0, double x1 );

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_UTILS_H_SEEN */
/*---------------------------------------------------------------------------*/






