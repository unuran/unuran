/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: unur_fp.c                                                         *
 *                                                                           *
 *   miscelleanous routines for floating point arithmetic                    *
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

#include <unur_source.h>

/*---------------------------------------------------------------------------*/

int
_unur_isfinite (const double x)
     /*----------------------------------------------------------------------*/
     /* Check whether x is a finite number.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x is finite                                              */
     /*   FALSE  otherwise (i.e., if x is +/- infinity or NaN)               */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISFINITE
  return (isfinite(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  if (x < INFINITY && x > -INFINITY)
    return TRUE;
  else
    return FALSE;
#else
# error
# error +--------------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! ....... !
# error ! Please contact <unuran@statistik.wu-wien.ac.at>. !
# error +--------------------------------------------------+
# error
#endif
} /* end of _unur_isfinite() */

/*---------------------------------------------------------------------------*/

int
_unur_isnan (const double x)
     /*----------------------------------------------------------------------*/
     /*  Check whether x is a NaN (not a number) value.                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   TRUE   if x is NaN                                                 */
     /*   FALSE  otherwise                                                   */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISNAN
  return (isnan(x) ? TRUE : FALSE);
#elif HAVE_IEEE_COMPARISONS
  return ((x!=x) ? TRUE : FALSE);
#else
# error
# error +--------------------------------------------------+
# error ! Sorry, Cannot handle NaN correctly! ............ !
# error ! Please contact <unuran@statistik.wu-wien.ac.at>. !
# error +--------------------------------------------------+
# error
#endif
} /* end of _unur_isnan() */

/*---------------------------------------------------------------------------*/

int
_unur_isinf (const double x)
     /*----------------------------------------------------------------------*/
     /*  Check whether x is infinity.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x ... floating-point number                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   -1  ... if x is -infinity                                          */
     /*    1  ... if x is +infinity                                          */
     /*    0  ... otherwise                                                  */
     /*----------------------------------------------------------------------*/
{
#if HAVE_DECL_ISINF
  if (isinf(x)==0)
    return 0;
  else if (x>0) 
    return 1;
  else
    return -1;
#elif HAVE_IEEE_COMPARISONS
  if (x>=INFINITY)
    return 1;
  else if (x<=-INFINITY)
    return -1;
  else
    return 0;
#else
# error
# error +--------------------------------------------------+
# error ! Sorry, Cannot handle INFINITY correctly! ....... !
# error ! Please contact <unuran@statistik.wu-wien.ac.at>. !
# error +--------------------------------------------------+
# error
#endif
} /* end of _unur_isinf() */

/*---------------------------------------------------------------------------*/
