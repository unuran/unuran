/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Distribution objects of multivariate distributions having special        *
 *  correlation matrices.                                                    *
 *                                                                           *
 *****************************************************************************
 
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

/*  Multinormal distribution (corr-matrix with equal off-diagonal elements)  */
UNUR_DISTR *unur_distr_multinormal_constant_rho(int dim, const double *mean, double rho);
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

void 
_unur_vc_special_correlations_set_ar1(struct unur_distr *distr, int dim, double rho) 
/* 
       Setting the covariance matrix to be
    
        1       r       r^2     ... r^(d-1) 
        r       1       r       ... r^(d-2) 
        r^2     r       1       ... r^(d-3) 
                                ... 
        r^(d-1) r^(d-2) r^(d-3) ... 1 
       
      in this case the inverse matrix is given as 
      
        a b 0 0 ... 0 0
        b c b 0 ... 0 0
        0 b c b ... 0 0
                ... 
        0 0 0 0 ... c b
        0 0 0 0 ... b a

      with 
        
        a = 1 / (1-r^2)
        b = - r / (1-r^2)
        c = (1+r^2) / (1-r^2)
    
      and the determinant of the covariance matrix is
      
        det = (1-r^2)^(d-1) 
*/    
{
#define idx(a,b) ((a)*dim+(b))
  
  double *covar = NULL;
  double *covar_inv = NULL;
  
  int i,j;
  double a,b,c,denominator;
  
  /* checking parameters */ 
  denominator=1-rho*rho;    
  if (fabs(denominator)<DBL_EPSILON || rho<0 || rho>=1 || dim<1) {
    distr = NULL;    
    return;
  }
      
  /* setting the covariance matrix */
  covar = malloc( dim * dim * sizeof(double) );
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = (i==j) ? 1.: pow(rho, abs(i-j));
  }}       
  unur_distr_cvec_set_covar( distr, covar );

    
  /* setting the inverse covariance matrix */
  covar_inv = malloc( dim * dim * sizeof(double) );
    
  a=1./denominator;
  b=-rho/denominator;
  c=(1+rho*rho)/denominator;
  
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = 0.;
    if (i==j &&  (i==0 || i==(dim-1))) covar_inv[idx(i,j)] = a ;
    if (i==j && !(i==0 || i==(dim-1))) covar_inv[idx(i,j)] = c ;
    if (abs(i-j)==1) covar_inv[idx(i,j)] = b ;
  }}   
  unur_distr_cvec_set_covar_inv( distr, covar_inv );
     
  free(covar); free(covar_inv);
    
  return;    

#undef idx
} /* end of _unur_vc_special_correlations_set_ar1() */

/*---------------------------------------------------------------------------*/

void 
_unur_vc_special_correlations_set_constant_rho(struct unur_distr *distr, int dim, double rho) 
/* 
       Setting the covariance matrix to be
    
        1 r r ... r 
        r 1 r ... r 
        r r 1 ... r 
              ... 
        r r r ... 1 
       
      in this case the inverse matrix is given as 
      
        a b b ... b
        b a b ... b
        b b a ... b
              ... 
        b b b ... a

      with 
        
        a = (1+(d-2)*r) / (1+(d-2)*r-(d-2)*r^2)
        b = - r / (1+(d-2)*r-(d-2)*r^2)
      
      and the determinant of the covariance matrix is
      
        det = (1-r)^(d-1) * (1+(d-1)*r) 
*/    
{
#define idx(a,b) ((a)*dim+(b))
  
  double *covar = NULL;
  double *covar_inv = NULL;
  
  int i,j;
  double a,b,denominator;
  
  /* checking parameters */ 
  denominator=1+(dim-2)*rho-(dim-1)*rho*rho;
  if (fabs(denominator)<DBL_EPSILON || rho<0 || rho>=1 || dim<1) {
    distr = NULL;    
    return;
  }
      
  /* setting the covariance matrix */
  covar = malloc( dim * dim * sizeof(double) );
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = (i==j) ? 1.: rho;
  }}       
  unur_distr_cvec_set_covar( distr, covar );

    
  /* setting the inverse covariance matrix */
  covar_inv = malloc( dim * dim * sizeof(double) );
      
  a=(1.+(dim-2)*rho)/denominator;
  b=-rho/denominator;
    
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = (i==j) ? a: b;
  }} 
  unur_distr_cvec_set_covar_inv( distr, covar_inv );
     
  free(covar); free(covar_inv);
    
  return;    

#undef idx
} /* end of _unur_vc_special_correlations_set_constant_rho() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal_ar1(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (corr-matrix from AR1 process)                  */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multinormal distribution */
  distr = unur_distr_multinormal( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multinormal_ar1");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_ar1(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal_ar1() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multicauchy_ar1(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (corr-matrix from AR1 process)                  */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multicauchy distribution */
  distr = unur_distr_multicauchy( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multicauchy_ar1");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_ar1(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multicauchy_ar1() */

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multistudent_ar1(int dim, double df, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multistudent distribution (corr-matrix from AR1 process)                 */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multistudent distribution */
  distr = unur_distr_multistudent( dim, df, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multistudent_ar1");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_ar1(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multistudent_ar1() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal_constant_rho(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (equal off diagonal elements of corr-matrix)    */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multinormal distribution */
  distr = unur_distr_multinormal( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multinormal_constant_rho");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_constant_rho(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multinormal_constant_rho() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multicauchy_constant_rho(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (equal off diagonal elements of corr-matrix)    */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multicauchy distribution */
  distr = unur_distr_multicauchy( dim, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multicauchy_constant_rho");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_constant_rho(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multicauchy_constant_rho() */

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multistudent_constant_rho(int dim, double df, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multistudent distribution (equal off diagonal elements of corr-matrix)   */
     /*---------------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  
  /* get distribution object for multistudent distribution */
  distr = unur_distr_multistudent( dim, df, mean, NULL );

  /* set the name of distribution */
  unur_distr_set_name(distr, "multistudent_constant_rho");
  
  /* set the correlation matrix and it's inverse */
  _unur_vc_special_correlations_set_constant_rho(distr, dim, rho);
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_multistudent_constant_rho() */

/*---------------------------------------------------------------------------*/



#undef DISTR
/*---------------------------------------------------------------------------*/
