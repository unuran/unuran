/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: empk.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method EMPK                                *
 *         (EMPirical distribution with Kernel smoothing)                    *
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
   =METHOD  EMPK   EMPirical distribution with Kernel smoothing

   .....

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_empk_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_empk_set_kernel( UNUR_PAR *parameters, UNUR_DISTR *kernel);
/* 
   Set generator for the kernel used the density estimation.
   Default is gaussian distribution.
*/

int unur_empk_set_positive(  UNUR_PAR *parameters, int positive );
/*
  If @code{positive} is TRUE then only nonnegative random variates are
  generated. This is done by means of mirroring technique.
  Default is FALSE.
*/

int unur_empk_set_varcor(  UNUR_PAR *parameters, int varcor );
/*
  Set whether the variance corrected version of the density estimation
  is used. If @code{varcor} is TRUE then the variance of the used
  density estimation is the same as the sample variance. However this 
  increases the efficiency of the estimation a little bit.
  Default is FALSE.
*/

int unur_empk_set_kernvar(  UNUR_PAR *parameters, double kernvar );
/*
  Variance of the used kernel. It is only required for the variance
  reduced version of the density estimation.
  Default is 1.
*/

int unur_empk_set_alfa(  UNUR_PAR *parameters, double alfa );
/*
  alfa is used to compute the optimal bandwidth from
  the point of view of minimizing the mean integrated
  square error (MISE).
  alfa depends on the type of kernel being used.     
  DEFAULT is ??.


  kannst du hier noch ein paar hinweise geben, was der benutzer
  eingeben soll??

  oder reicht ein 
  There is no need to change this factor it almost all situations.
*/

