/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: x_misc.c                                                          *
 *                                                                           *
 *   miscellaneous routines                                                  *
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

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
#ifdef HAVE_LIBMD    /* have CEPHES library                                  */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* The CEPHES library writes error messages (over/underflow) to stdout       */
/* using the command mtherr().                                               */
/* This behavior is annoying. So we override mtherr() with our own version   */
/* that does nothing. It should work as long as libunuran is linked          */
/* _before_ libmd (which must be done in any way to link all necessary       */
/* routines from the library).                                               */
/*---------------------------------------------------------------------------*/
int
mtherr( char *name, int code )
{
  return 0;
} /* end of mtherr() */

/*---------------------------------------------------------------------------*/
#endif  /* HAVE_LIBMD */
/*---------------------------------------------------------------------------*/



