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
   Set kernel distribution.
   Default is a Gaussian kernel.
*/

int unur_empk_set_kernelgen( UNUR_PAR *parameters, UNUR_GEN *kernelgen);
/* 
   Set generator for the kernel used the density estimation.
   It is not necessary to set both the kernel and the kernel generator.
   If no kernel generator is provided one is automatically created
   from the distribution object given by the unur_empk_set_kernel()
   call. On the other hand the unur_empk_set_kernelgen() overwrites
   the unur_empk_set_kernel() call.
   Notice that the uniform random number generator of the kernel
   generator is overwritten during the unur_init() and at each
   unur_chg_urng() call with generator for the empirical
   distribution.
   Default is a Gaussian kernel.
*/

int unur_empk_set_alpha( UNUR_PAR *parameters, double alpha );
/*
  alpha is used to compute the optimal bandwidth from the point of
  view of minimizing the mean integrated square error (MISE).
  alpha depends on the type of kernel K being used and is given by 
     alpha(K) = Var(K)^(-2/5){ \int K(t)^2 dt}^(1/5)
  For standard kernels (see ????) alpha is computed by the algorithm.
  For all other kernels, it must be given. 
*/
/*
  Otherwise the default (Gaussian) kernel is used.
*/

int unur_empk_set_beta( UNUR_PAR *parameters, double beta );
/*
  beta is used to compute the optimal bandwidth from the point of
  view of minimizing the mean integrated square error (MISE).
  beta depends on the (unknown) distribution of the sampled data
  points. Thus its value contains some guess on this distribution.
  By default Gaussian distribution is assumed for the sample
  (beta = 1.3637439). There is no requirement to set beta.
*/

int unur_empk_set_smoothing( UNUR_PAR *parameters, double smoothing );
/*
  The smoothing factor controlles how "smooth" the resulting density
  estimation will be. A smoothing factor equal to 0 results in naive
  resampling. A very large smoothing factor (together with the
  variance correction) results in a density which is approximately
  equal to the kernel.
  Default is 1.
*/

int unur_empk_set_varcor( UNUR_PAR *parameters, int varcor );
/*
  Set whether the variance corrected version of the density estimation
  is used. If @code{varcor} is TRUE then the variance of the used
  density estimation is the same as the sample variance. However this 
  increases the efficiency of the estimation a little bit.
  Default is TRUE.
*/

int unur_empk_set_kernelvar( UNUR_PAR *parameters, double kernvar );
/*
  Variance of the used kernel. It is only required for the variance
  reduced version of the density estimation (which is used by default).
  For standard kernels (see ????) kernvar is computed by the algorithm.
  Default is 1.
*/

int unur_empk_set_positive( UNUR_PAR *parameters, int positive );
/*
  If @code{positive} is TRUE then only nonnegative random variates are
  generated. This is done by means of mirroring technique.
  Default is FALSE.
*/

/*---------------------------------------------------------------------------*/
