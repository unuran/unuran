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

   EMPK generates random variates from an empirical distribution that is
   given by an observed sample. The idea is that simply choosing a random
   point from the sample and to return it with some added noise results
   in a method that has very nice properties, as it can be seen as sampling
   from a kernel density estimate.
   Clearly we have to decide about the density of the noise (called kernel)
   and about the standard deviation of the noise.
   The mathematical theory of kernel density estimation shows us that we
   are comparatively free in choosing the kernel. It also supplies us with
   a simple formula to compute the optimal standarddeviation of the noise,
   called bandwidth (or window width) of the kernel.
   For most applications it is perfectly ok to use the default values offered.
   Unless you have some knowledge on density estimation we do not recommend
   to change anything. Only exception is the case that you are especially
   interested in a fast sampling algorithm. Then use the call

   unur_empk_set_kernel( par, UNUR_DISTR_BOXCAR);

   to change the used noise distribution from the default Gaussian
   distribution to the uniform distribution.

   All other parameters are only necessary for people knowing the theory
   of kernel density estimation.

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_empk_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_empk_set_kernel( UNUR_PAR *parameters, unsigned kernel);
/* 
   Select one of the supported kernel distributions. Currently the following 
   kernels are supported:

     UNUR_DISTR_GAUSSIAN  ... Gaussian (normal) kernel
     UNUR_DISTR_BOXCAR    ... Boxcar (uniform, rectangular) kernel
     UNUR_DISTR_STUDENT   ... t3 kernel (Student's distribution with 3 degrees of freedom)
     UNUR_DISTR_LOGISTIC  ... logistic kernel

   For other kernels (including kernels with Student's distribution other with
   3 degrees of freedom) use the unur_empk_set_kernelgen() call.
   However then unur_empk_set_alpha() and (if variance correction is used)
   unur_empk_set_kernelvar() has to be called.

   It is not possible to call unur_empk_set_kernel() twice.

   Default is a Gaussian kernel.
*/

int unur_empk_set_kernelgen( UNUR_PAR *parameters, UNUR_GEN *kernelgen);
/* 
   Set generator for the kernel used for density estimation.
   It is not possible to call unur_empk_set_kernelgen() after a standard kernel
   has been selected by a unur_empk_set_kernel() call.
   After calling this function unur_empk_set_alpha() and (if variance
   correction is used) unur_empk_set_kernelvar() has to be called!

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
     alpha(K) = Var(K)^(-2/5)@{ \int K(t)^2 dt@}^(1/5)
  For standard kernels (see above) alpha is computed by the algorithm.
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
  Default is 1 which results in a smoothing parameter minimising
  the MISE (mean integrated squared error) if the data are not too
  far away from normal.
*/

int unur_empk_chg_smoothing( UNUR_GEN *generator, double smoothing );
/* 
   Change smoothing factor in generator.
*/

int unur_empk_set_varcor( UNUR_PAR *parameters, int varcor );
/*
  Set whether the variance corrected version of the density estimation
  is used. If @code{varcor} is TRUE then the variance of the used
  density estimation is the same as the sample variance. However this 
  increases the MISE of the estimation a little bit.
  Default is TRUE.
*/

int unur_empk_chg_varcor( UNUR_GEN *generator, int varcor );
/* 
   Switch variance correction in generator on/off.
   Default is FALSE.
*/

int unur_empk_set_kernelvar( UNUR_PAR *parameters, double kernvar );
/*
  Variance of the used kernel. It is only required for the variance
  reduced version of the density estimation (which is used by default).
  For standard kernels (see above) kernvar is computed by the algorithm.
  Default is 1.
*/

int unur_empk_set_positive( UNUR_PAR *parameters, int positive );
/*
  If @code{positive} is TRUE then only nonnegative random variates are
  generated. This is done by means of mirroring technique.
  Default is FALSE.
*/

/* =END */

/*---------------------------------------------------------------------------*/
