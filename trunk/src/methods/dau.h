/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dau.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DAU                                *
 *         ((Discrete) Alias-Urn)                                            *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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

/* 
   =METHOD  DAU  (Discrete) Alias-Urn method

   DAU samples from arbitrary but finite probability vectors of length
   N. The algorithmus is based on an ingeneous method by A.J. Walker
   and requires a table of size (at least) N and needs only one
   comparison for each generated random variate.
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dau_new( UNUR_DISTR *distribution );
/* Get default parameters for generator.                                     */

/*...........................................................................*/

int unur_dau_set_urnfactor( UNUR_PAR *parameters, double factor );
/* 
   Set size of urn table relative to length of probability vector.  It
   must not be less than 1. Larger tables result in (slightly) faster
   generation times but require a more expensive setup. However sizes
   larger than 2 are not recommended; default is 1.
*/

/* =END */

/*---------------------------------------------------------------------------*/










