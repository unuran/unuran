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

   =UP  Methods_for_CVEMP

   =DESCRIPTION
      VEMPK generates random variates from a multivariate empirical
      distribution that is given by an observed sample. The idea is
      that simply choosing a random point from the sample and to
      return it with some added noise results in a method that has
      very nice properties, as it can be seen as sampling from a
      kernel density estimate. 
      Clearly we have to decide about the density of the noise (called kernel)
      and about the standard deviation of the noise.
      The mathematical theory of kernel density estimation shows us that we
      are comparatively free in choosing the kernel. 
      It also supplies us with a simple formula to compute the optimal
      standarddeviation of the noise, called bandwidth (or window
      width) of the kernel.

      Currently only a Gaussian kernel with the same covariance matrix
      as the given sample is implemented.
      However it is possible to choose between a variance corrected
      version or those with optimal MISE.
      Additionally a smoothing factor can be set.

   =END

*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vempk_new( UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_vempk_set_smoothing( UNUR_PAR *parameters, double smoothing );
/* */

int unur_vempk_chg_smoothing( UNUR_GEN *generator, double smoothing );
/* 
   Set and change the smoothing factor.
   The smoothing factor controlles how ``smooth'' the resulting density
   estimation will be. A smoothing factor equal to 0 results in naive
   resampling. A very large smoothing factor (together with the
   variance correction) results in a density which is approximately
   equal to the kernel.
   Default is 1 which results in a smoothing parameter minimising
   the MISE (mean integrated squared error) if the data are not too
   far away from normal.
*/

int unur_vempk_set_varcor( UNUR_PAR *parameters, int varcor );
/* */

int unur_vempk_chg_varcor( UNUR_GEN *generator, int varcor );
/*
   Switch variance correction in generator on/off.
   If @var{varcor} is TRUE then the variance of the used
   density estimation is the same as the sample variance. However this 
   increases the MISE of the estimation a little bit.

   Default is TRUE.
*/

/* =END */

/*---------------------------------------------------------------------------*/
