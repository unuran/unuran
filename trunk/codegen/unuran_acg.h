/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   file:      unuran_acg                                                   *
 *                                                                           *
 *   declarations for Automatic Code Generator                               *
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
#ifndef __UNURAN_ACG_H_SEEN
#define __UNURAN_ACG_H_SEEN
/*---------------------------------------------------------------------------*/

int unur_acg_C( struct unur_gen *gen, FILE *out, const char *distr_name );
/*---------------------------------------------------------------------------*/
/* Automatic code generator (C version)                                      */
/*---------------------------------------------------------------------------*/

int unur_acg_FORTRAN( struct unur_gen *gen, FILE *out, const char *distr_name );
/*---------------------------------------------------------------------------*/
/* Automatic code generator (FORTRAN version)                                */
/*---------------------------------------------------------------------------*/

int unur_acg_JAVA( struct unur_gen *gen, FILE *out, const char *distr_name );
/*---------------------------------------------------------------------------*/
/* Automatic code generator (JAVA version)                                   */
/*---------------------------------------------------------------------------*/

int unur_acg_UNURAN( struct unur_gen *gen, FILE *out, const char *distr_name );
/*---------------------------------------------------------------------------*/
/* Automatic code generator (UNURAN version)                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif   /* __UNURAN_ACG_H_SEEN */ 
/*---------------------------------------------------------------------------*/




