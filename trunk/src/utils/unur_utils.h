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
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Mon Oct 18 20:46:18 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
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

#define CHECK_NULL(ptr,rval)

#endif

/*---------------------------------------------------------------------------*/
/* Check number of arguments                                                 */

#if UNUR_DEBUG & UNUR_DB_CHECKARGS

#define CHECK_N_PARAMS(obsed,expted,rval)    \
  if (obsed != expted ) {                    \
     _unur_error(NULL,UNUR_ERR_NPARAM,"");   \
     return rval;                            \
  }

#else               /* do not check (be carefull) */

#define CHECK_N_PARAMS(obsed,expted,rval)

#endif 

/* an abbreviation */
#define CHECKARGS (UNUR_DEBUG & UNUR_DB_CHECKARGS)

/*---------------------------------------------------------------------------*/

double _unur_arcmean( double x0, double x1 );

/*---------------------------------------------------------------------------*/
#endif  /* __UNUR_UTILS_H_SEEN */
/*---------------------------------------------------------------------------*/






