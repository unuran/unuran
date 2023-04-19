/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Distribution objects of standard distributions without logPDF            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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

#include <unuran.h>
#include <unuran_config.h>
#include <unur_struct.h>
#include "testdistributions.h"

/*---------------------------------------------------------------------------*/
/* Continuous univariate distributions without logPDF                        */
/*---------------------------------------------------------------------------*/
#define DISTR distr->data.cont
/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_beta_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for beta distribution */
  distr = unur_distr_beta( params, n_params );

  /* name of distribution */
  distr->name = "beta_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_beta_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cauchy_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for cauchy distribution */
  distr = unur_distr_cauchy( params, n_params );

  /* name of distribution */
  distr->name = "cauchy_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_cauchy_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_exponential_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for exponential distribution */
  distr = unur_distr_exponential( params, n_params );

  /* name of distribution */
  distr->name = "exponential_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_exponential_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_gamma_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for gamma distribution */
  distr = unur_distr_gamma( params, n_params );

  /* name of distribution */
  distr->name = "gamma_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_gamma_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_laplace_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for laplace distribution */
  distr = unur_distr_laplace( params, n_params );

  /* name of distribution */
  distr->name = "laplace_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_laplace_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_normal_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for normal distribution */
  distr = unur_distr_normal( params, n_params );

  /* name of distribution */
  distr->name = "normal_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_normal_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_powerexponential_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for powerexponential distribution */
  distr = unur_distr_powerexponential( params, n_params );

  /* name of distribution */
  distr->name = "powerexponential_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_powerexponential_wo_logpdf() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_uniform_wo_logpdf( const double *params, int n_params )
{
  register struct unur_distr *distr;

  /* get distribution object for uniform distribution */
  distr = unur_distr_uniform( params, n_params );

  /* name of distribution */
  distr->name = "uniform_wo_logpdf";

  /* disable logPDF */
  DISTR.logpdf  = NULL;
  DISTR.dlogpdf = NULL;

  /* return pointer to object */
  return distr;

} /* end of unur_distr_uniform_wo_logpdf() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
