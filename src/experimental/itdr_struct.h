/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: itdr_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method ITDR                               *
 *         (Inverse Transformed Density Rejection)                           *
 *                                                                           *
 *****************************************************************************
     $Id: nrou_struct.h 2588 2005-11-05 11:47:58Z leydold $
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
/* Information for constructing the generator                                */

struct unur_itdr_par { 
  double bx;                 /* splitting point between pole and tail region */
  double cp, ct;             /* c-value for pole and tail region, resp.      */
};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_itdr_gen { 
  double bx;                 /* splitting point between pole and tail region */
  double cp, ct;             /* c-value for pole and tail region, resp.      */
  double xp, xt;             /* desig point in pole and tail region, resp.   */
  double alphap, betap;      /* parameters for hat in pole region            */
  double Tfxt, dTfxt;        /* parameters for hat in tail region            */
  double by;                 /* hat of pole region at bx                     */
  double Ap, Ac, At;         /* areas in upper pole, center, and tail region */     
  double Atot;               /* total area below hat                         */     
};

/*---------------------------------------------------------------------------*/
