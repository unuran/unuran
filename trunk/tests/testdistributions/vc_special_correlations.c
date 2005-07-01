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
#include <unur_source.h>
#include "testdistributions.h"

#define DISTR distr->data.cvec
#define LOGNORMCONSTANT (distr->data.cvec.norm_constant)

/*  Multinormal distribution (corr-matrix with equal off-diagonal elements)  */
UNUR_DISTR *unur_distr_multinormal_constant_rho(int dim, const double *mean, double rho);
/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_multinormal_ar1(int dim, const double *mean, double rho)
     /*---------------------------------------------------------------------------*/
     /*  Multinormal distribution (corr-matrix from AR1 process)                  */
     /*---------------------------------------------------------------------------*/
     /* 
       The covariance matrix is of the following form
    
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
	
      (this determinant occurs in the normalization constant)	
    */    
{
#define idx(a,b) ((a)*dim+(b))
  
  register struct unur_distr *distr;
  
  double *covar;
  double *covar_inv;
  
  int i,j;
  double a,b,c,denominator;
  double det_covar;
  
  /* checking parameters */
  if (dim<1) {
  
  }
  
  if (rho<0 || rho>=1) {

  }

  denominator=1-rho*rho;    
  if (fabs(denominator)<UNUR_EPSILON) {
    _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"cannot compute inverse of covariance");
    return NULL;    
  }
      

  /* get distribution object for multinormal distribution */
  distr = unur_distr_multinormal( dim, mean, NULL );

  /* name of distribution */
  distr->name = "multinormal_ar1";

  
  /* setting the covariance matrix */
  covar = _unur_xmalloc( dim * dim * sizeof(double) );
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = (i==j) ? 1.: pow(rho, abs(i-j));
  }}       
  unur_distr_cvec_set_covar( distr, covar );

    
  /* setting the inverse covariance matrix */
  covar_inv = _unur_xmalloc( dim * dim * sizeof(double) );
    
  a=1./denominator;
  b=-rho/denominator;
  c=(1+rho*rho)/denominator;
  
  for (i=0; i<dim; i++) {
  for (j=0; j<dim; j++) {
    covar_inv[idx(i,j)] = 0.;
    if (i==j &&  (i==0 || i==(dim-1))) DISTR.covar_inv[idx(i,j)] = a ;
    if (i==j && !(i==0 || i==(dim-1))) DISTR.covar_inv[idx(i,j)] = c ;
    if (abs(i-j)==1) DISTR.covar_inv[idx(i,j)] = b ;
  }}   
  unur_distr_cvec_set_covar_inv( distr, covar_inv );
  

  /* update log of normalization constant */
  det_covar = pow(1-rho*rho, dim-1);
  LOGNORMCONSTANT = - ( dim * log(2 * M_PI) + log(det_covar) ) / 2.;  
  
  free(covar); free(covar_inv);
  
  /* return pointer to object */
  return distr;

#undef idx
} /* end of unur_distr_multinormal_ar1() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
