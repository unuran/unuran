/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vempk.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VEMPK                              *
 *         ((Vector) EMPirical distribution with Kernel smoothing)           *
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
   =METHOD  VEMPK   (Vector) EMPirical distribution with Kernel smoothing

*/

#if 0
   VEMPK generates random variates from an empirical distribution that is
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

 double smooth;/*should be mostly between 0 and 1, controls the smoothing*/ 
 /*1. smooth estimate optimal for normal distribution, 0. no noise added
   for bimodal distributions a value smaller than 1 should be used to
   avoid oversmoothing;
   if smooth is chosen very big (eg. 1000) together with varcor=1.
   this results approximately in fitting the kernel distribution (multi
   normal distribution) to the data*/



#endif

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vempk_new( UNUR_DISTR *distribution );
/* Get default parameters for generator                                      */

/*...........................................................................*/

int unur_vempk_set_smoothing( UNUR_PAR *parameters, double smoothing );
/*
  The smoothing factor controlles how "smooth" the resulting density
  estimation will be. A smoothing factor equal to 0 results in naive
  resampling. A very large smoothing factor (together with
  variance correction) results in a density which is approximately
  equal to the kernel.
  Default is 1 which results in a smoothing parameter minimising
  the MISE (mean integrated squared error) if the data are not too
  far away from normal.
*/

int unur_vempk_chg_smoothing( UNUR_GEN *generator, double smoothing );
/* 
   Change smoothing factor in generator.
*/

int unur_vempk_set_varcor( UNUR_PAR *parameters, int varcor );
/*
  Set whether the variance corrected version of the density estimation
  is used. If @code{varcor} is TRUE then the variance of the used
  density estimation is the same as the sample variance. However this 
  increases the MISE of the estimation a little bit.
  Default is TRUE.
*/

int unur_vempk_chg_varcor( UNUR_GEN *generator, int varcor );
/* 
   Switch variance correction in generator on/off.
   Default is FALSE.
*/

/* =END */

/*---------------------------------------------------------------------------*/
