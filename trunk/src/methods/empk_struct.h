/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: empk_struct.h                                                     *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         declares structures for method EMPK                               *
 *         (EMPirical distribution with Kernel smoothing)                    *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in source_struct.h                                  *
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
  n -- > n_observ

*/


/*---------------------------------------------------------------------------*/
/* Information for constructing the generator                                */

struct unur_empk_par {
  /* the observed sample is stored in the distribution object */

  int     positive;    /* 1 ... only positive variates (mirroring)
			  0 ... negative variates possible                   */
  int     varcor;      /* 1 ... use variance correction
			  0 ... no variance correction                       */

  double (*kernrvg)(UNUR_GEN *gen);  /* random variate generator for kernel  */

  double  alfa;        /* alfa is used to compute the optimal bandwidth from
			  the point of view of minimizing the mean integrated
			  square error (MISE).
			  alfa depends on the type of kernel being used.     */

  double  kernvar;     /* variance of used kernel, only used if varcor == 1  */

};

/*---------------------------------------------------------------------------*/
/* The generator object                                                      */

struct unur_empk_gen {
  double *observ;      /* pointer to the array of the observations           */
  int     n_observ;    /* number of observations                             */

  int     positive;    /* 1 ... only positive variates (mirroring)
			  0 ... negative variates possible                   */
  int     varcor;      /* 1 ... use variance correction
			  0 ... no variance correction                       */

  double (*kernrvg)(UNUR_GEN *gen);  /* random variate generator for kernel  */

  double  bwidth;      /* bandwidth for kernel density estimation            */
  double  xbar;        /* sample mean                                        */
  double  stdev;       /* sample standard deviation; not used in sample!     */
  double  sconst;      /* constant used for variance corrected version 
			  of kernel method                                   */
};

/*---------------------------------------------------------------------------*/
