/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Distribution objects of multivariate distributions having a              *
 *  correlation matrix with constant off-diagonal elements                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Use correlation matrix of the form:                                      *
 *                                                                           *
 *      +-            -+                                                     *
 *      |  1 r r ... r |                                                     *
 *      |  r 1 r ... r |                                                     *
 *      |  r r 1 ... r |                                                     *
 *      |  ...         |                                                     *
 *      |  r r r ... 1 |                                                     *
 *      +-            -+                                                     *
 *                                                                           *
 *     in this case the inverse matrix is given as                           *
 *                                                                           *
 *      +-              -+                                                   *
 *      |  a b b ... b b |                                                   *
 *      |  b a b ... b b |                                                   *
 *      |  b b a ... b b |                                                   *
 *      |        ...     |                                                   *
 *      |  b b b ... a b |                                                   *
 *      |  b b b ... b a |                                                   *
 *      +-              -+                                                   *
 *                                                                           *
 *     with                                                                  *
 *                                                                           *
 *       a = (1+(d-2)*r) / (1+(d-2)*r-(d-1)*r^2)                             *
 *       b = - r / (1+(d-2)*r-(d-1)*r^2)                                     *
 *                                                                           *
 *     and the determinant of the covariance matrix is                       *
 *                                                                           *
 *       det = (1-r)^(d-1) * (1+(d-1)*r)                                     *
 *                                                                           *
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

#include <unuran.h>
#include <float.h>
#include "testdistributions.h"

/*---------------------------------------------------------------------------*/

/* Set correlation matrix with constant off-diagonal elements and its inverse */
int _unur_vc_set_corrmatrix_constantrho( UNUR_DISTR *distr, int dim, double rho ); 

/*---------------------------------------------------------------------------*/

int 
_unur_vc_set_corrmatrix_constantrho( UNUR_DISTR *distr, int dim, double rho )
     /*----------------------------------------------------------------------*/
     /* Create and set correlation matrix with constant off-diagonal         */
     /* elements and its inverse.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dim   ... dimension of random vector                               */
     /*   rho   ... correlation between consecutive elements in process      */
     /*                                                                      */
     /* return:                                                              */
     /*   0 ... success                                                      */
     /*   1 ... failure                                                      */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  
  double *covar, *covar_inv;
  double a, b, denominator;
  int i,j;
  int error = 1;
  
  /* checking parameters */ 
  if (rho == 0.) /* nothing to do */
    return 0;

  /* entries of matrix */
  denominator = 1.+(dim-2)*rho-(dim-1)*rho*rho;    
  a = (1.+(dim-2)*rho)/denominator;
  b = -rho/denominator;

  if ( rho<0. || rho>=1. || dim<1)
    /* invalid arguments */
    return 1;

  /* allocate memory */
  covar = malloc( dim * dim * sizeof(double) );
  covar_inv = malloc( dim * dim * sizeof(double) );
  if (covar == NULL || covar_inv == NULL) {
    if (covar) free (covar);
    if (covar_inv) free (covar_inv);
    return 1;
  }

  /* create correlation matrix */
  for (i=0; i<dim; i++) 
    for (j=0; j<dim; j++)
      covar[idx(i,j)] = (i==j)? 1.: rho;

  /* create inverse of correlation matrix */
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++) 
      covar_inv[idx(i,j)] = (i==j) ? a : b;
    
  
  /* set covariance matrix and its inverse */
  if ( (unur_distr_cvec_set_covar( distr, covar ) == UNUR_SUCCESS) &&
       (unur_distr_cvec_set_covar_inv( distr, covar_inv ) == UNUR_SUCCESS) )
    error = 0;

  /* free working space */
  free (covar); free (covar_inv); 

  /* return errorcode */
  return error;    

#undef idx
} /* end of _unur_vc_set_corrmatrix_constantrho() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal_constantrho(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (corr-matrix : constant off-diagonal elements)  */
     /*---------------------------------------------------------------------------*/
{
  struct unur_distr *distr;

  /* checking parameters */ 
  if ( rho<0. || rho>=1. || dim<1)
    return NULL;

  /* get distribution object for multinormal distribution */
  distr = unur_distr_multinormal( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multinormal_constantrho");
  
  /* set the correlation matrix and its inverse */
  if (_unur_vc_set_corrmatrix_constantrho(distr, dim, rho) != 0) {
    /* error */
    unur_distr_free(distr); distr = NULL; 
  }
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal_constantrho() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multicauchy_constantrho(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (corr-matrix : constant off-diagonal elements)  */
     /*---------------------------------------------------------------------------*/
{
  struct unur_distr *distr;
  
  /* checking parameters */ 
  if ( rho<0. || rho>=1. || dim<1)
    return NULL;

  /* get distribution object for multicauchy distribution */
  distr = unur_distr_multicauchy( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multicauchy_constantrho");
  
  /* set the correlation matrix and its inverse */
  if (_unur_vc_set_corrmatrix_constantrho(distr, dim, rho) != 0) {
    /* error */
    unur_distr_free(distr); distr = NULL; 
  }
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multicauchy_constantrho() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multistudent_constantrho(int dim, double df, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multistudent distribution (corr-matrix : constant off-diagonal elements) */
     /*---------------------------------------------------------------------------*/
{
  struct unur_distr *distr;
  
  /* checking parameters */ 
  if ( rho<0. || rho>=1. || dim<1)
    return NULL;

  /* get distribution object for multistudent distribution */
  distr = unur_distr_multistudent( dim, df, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multistudent_constantrho");
  
  /* set the correlation matrix and its inverse */
  if (_unur_vc_set_corrmatrix_constantrho(distr, dim, rho) != 0) {
    /* error */
    unur_distr_free(distr); distr = NULL; 
  }
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multistudent_constantrho() */

/*---------------------------------------------------------------------------*/
