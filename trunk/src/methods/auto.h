/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: auto.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method AUTO                               *
 *         (selects a method for a given distribution object AUTOmatically)  *
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
   =METHOD  AUTO  Select method automatically

   =UP Methods [05]

   =DESCRIPTION
      AUTO selects a an appropriate method for the given distribution
      object automatically. There are no parameters for this method,
      yet. But it is planned to give some parameter to describe the
      task for which the random variate generator is used for and thus
      make the choice of the generating method more appropriate.
      Notice that the required sampling routine for the generator
      object depends on the type of the given distribution object.

      IMPORTANT: This is an experimental version and the method chosen
      may change in future releases of UNURAN.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_auto_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/


/* =END */

/*---------------------------------------------------------------------------*/












