/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cvec.c                                                       *
 *                                                                           *
 *   manipulate multivariate continuous distribution objects                 *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distributions/unur_stddistr.h>
#include "distr_source.h"
#include "distr.h"
#include "cvec.h"
#include <utils/matrix_source.h>

/*---------------------------------------------------------------------------*/

static const char unknown_distr_name[] = "unknown";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec

/* Inverse of covariance matrix cannot be computed it is ill-conditioned.    */
/* We use the threshold |det(A)| / (dim * ||A||) < COVARIANCE_DETMIN.        */
/* (see _unur_matrix_invert() in ./src/utils/matrix.c)                       */
#define COVARIANCE_DETMIN  (1.e-10) 

/*---------------------------------------------------------------------------*/

static void _unur_distr_cvec_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** mulitvariate continuous distributions                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cvec_new( int dim )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: multivariate continuous with given PDF                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... number of components of random vector (dimension)          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  int i;

  /* check dimension for new parameter for distribution */
  if (dim < 2) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"dimension < 2");
    return NULL;
  }

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CVEC);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CVEC;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = dim;   /* mulitvariant */

  /* name of distribution */
  distr->name = unknown_distr_name;
  distr->name_str = NULL;

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cvec_free;

  /* clone */
  distr->clone = _unur_distr_cvec_clone;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;   /* pointer to PDF                                */
  DISTR.dpdf      = NULL;   /* pointer to gradient of PDF                    */
  DISTR.init      = NULL;   /* pointer to special init routine (default: none) */
  DISTR.mean      = NULL;   /* pointer to mean vector (default: not known)   */
  DISTR.covar     = NULL;   /* pointer to covariance matrix (default: not known) */
  DISTR.cholesky  = NULL;   /* pointer to cholesky factor (default: not computed) */
  DISTR.covar_inv = NULL;   /* pointer to inverse covariance matrix 
			       (default: not computed ) */
  DISTR.rankcorr  = NULL;   /* rank correlation of distribution              */
  DISTR.marginals = NULL;   /* array of pointers to marginal distributions   */

  /* initialize parameters of the PDF                                        */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    DISTR.n_params[i] = 0;
    DISTR.params[i] = NULL;
  }

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for PDF
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode       = NULL;         /* location of mode (default: not known)  */
  DISTR.volume     = INFINITY;     /* area below PDF (default: not known)    */


  distr->set = 0u;          /* no parameters set                             */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cvec_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cvec_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.cvec

  struct unur_distr *clone;
  int i,len;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* allocate memory */
  clone = _unur_malloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy data about distribution */
  if (DISTR.mean) {
    CLONE.mean = _unur_malloc( distr->dim * sizeof(double) );
    memcpy( CLONE.mean, DISTR.mean, distr->dim * sizeof(double) );
  }

  if (DISTR.covar) {
    CLONE.covar = _unur_malloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.covar, DISTR.covar, distr->dim * distr->dim * sizeof(double) );
  }

  if (DISTR.cholesky) {
    CLONE.cholesky = _unur_malloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.cholesky, DISTR.cholesky, distr->dim * distr->dim * sizeof(double) );
  }

  if (DISTR.covar_inv) {
    CLONE.covar_inv = _unur_malloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.covar_inv, DISTR.covar_inv, distr->dim * distr->dim * sizeof(double) );
  }

  if (DISTR.rankcorr) {
    CLONE.rankcorr = _unur_malloc( distr->dim * distr->dim * sizeof(double) );
    memcpy( CLONE.rankcorr, DISTR.rankcorr, distr->dim * distr->dim * sizeof(double) );
  }

  if (DISTR.mode) {
    CLONE.mode = _unur_malloc( distr->dim * sizeof(double) );
    memcpy( CLONE.mode, DISTR.mode, distr->dim * sizeof(double) );
  }

  if (DISTR.marginals) {
    CLONE.marginals = _unur_malloc (distr->dim * sizeof(struct unur_distr *));
    /* There are (should be) only two possibilities: 
       either all entries in the array point to the same distribution object;
          (set by unur_distr_cvec_set_marginals() call)
       or each entry has its own copy of some distribution object.
          (set by unur_distr_cvec_set_marginal_array() call)
    */
    if (DISTR.marginals[0] == DISTR.marginals[1]) {
      CLONE.marginals[0] = _unur_distr_clone( DISTR.marginals[0] );
      for (i=1; i<distr->dim; i++)
	CLONE.marginals[i] = CLONE.marginals[0];
    }
    else {
      for (i=0; i<distr->dim; i++) 
	CLONE.marginals[i] = _unur_distr_clone( DISTR.marginals[i] );
    }
  }
  
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    if (DISTR.params[i]) {
      CLONE.params[i] = _unur_malloc( DISTR.n_params[i] * sizeof(double) );
      memcpy( CLONE.params[i], DISTR.params[i], DISTR.n_params[i] * sizeof(double) );
    }
  }

  /* copy user name for distribution */
  if (distr->name_str) {
    len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_malloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_cvec_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cvec_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;

  COOKIE_CHECK(distr,CK_DISTR_CVEC,RETURN_VOID);

  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    if (DISTR.params[i]) free( DISTR.params[i] );

  if (DISTR.mean)      free(DISTR.mean); 
  if (DISTR.covar)     free(DISTR.covar);
  if (DISTR.covar_inv) free(DISTR.covar_inv);
  if (DISTR.cholesky)  free(DISTR.cholesky);
  if (DISTR.rankcorr)  free(DISTR.rankcorr);

  if (DISTR.mode)  free(DISTR.mode);

  if (DISTR.marginals) {
    /* There are (should be) only two possibilities: 
       either all entries in the array point to the same distribution object;
          (set by unur_distr_cvec_set_marginals() call)
       or each entry has its own copy of some distribution object.
          (set by unur_distr_cvec_set_marginal_array() call)
    */
    if (DISTR.marginals[0] == DISTR.marginals[1]) {
      _unur_distr_free(DISTR.marginals[0]);
    }
    else {
      for (i=0; i<distr->dim; i++) 
	_unur_distr_free(DISTR.marginals[i]);
    }
    free (DISTR.marginals);
  }

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_cvec_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CVEC *pdf )
     /*----------------------------------------------------------------------*/
     /* set PDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to PDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdf, UNUR_ERR_NULL);
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  /* we do not allow overwriting a PDF */
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return UNUR_SUCCESS;

} /* end of unur_distr_cvec_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_dpdf( struct unur_distr *distr, UNUR_VFUNCT_CVEC *dpdf )
     /*----------------------------------------------------------------------*/
     /* set gradient of PDF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dpdf  ... pointer to gradient of PDF                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a dPDF */
  if (DISTR.dpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_dpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CVEC *
unur_distr_cvec_get_pdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to PDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to PDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  return DISTR.pdf;
} /* end of unur_distr_cvec_get_pdf() */

/*---------------------------------------------------------------------------*/

UNUR_VFUNCT_CVEC *
unur_distr_cvec_get_dpdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to gradient of PDF of distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to gradient of PDF                                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  return DISTR.dpdf;
} /* end of unur_distr_cvec_get_dpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cvec_eval_pdf( const double *x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PDF of distribution at x                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for pdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );

  if (DISTR.pdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cvec_PDF(x,distr);
} /* end of unur_distr_cvec_eval_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_eval_dpdf( double *result, const double *x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate gradient of PDF of distribution at x                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   result ... to store grad (PDF(x))                                  */
     /*   x      ... argument for dPDF                                       */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  return _unur_cvec_dPDF(result,x,distr);
} /* end of unur_distr_cvec_eval_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_mean( struct unur_distr *distr, const double *mean )
     /*----------------------------------------------------------------------*/
     /* set mean vector of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mean  ... mean vector of distribution                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  /* we have to allocate memory first */
  if (DISTR.mean == NULL)
    DISTR.mean = _unur_malloc( distr->dim * sizeof(double) );

  if (mean)
    /* mean vector given --> copy */
    memcpy( DISTR.mean, mean, distr->dim * sizeof(double) );

  else  /* mean == NULL --> use zero vector instead */
    for (i=0; i<distr->dim; i++)
      DISTR.mean[i] = 0.;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MEAN;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_mean() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_mean( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mean vector of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to mean of distribution                                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* mean vector known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MEAN) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"mean");
    return NULL;
  }

  return DISTR.mean;

} /* end of unur_distr_cvec_get_mean() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_covar( struct unur_distr *distr, const double *covar )
     /*----------------------------------------------------------------------*/
     /* set covariance matrix of distribution.                               */
     /* as a side effect it also computes its cholesky factor.               */ 
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   covar ... covariance matrix of distribution                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int i,j;
  int dim;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  dim = distr->dim;

  /* mark as unknown */
  distr->set &= ~(UNUR_DISTR_SET_COVAR | UNUR_DISTR_SET_CHOLESKY | UNUR_DISTR_SET_COVAR_INV);

  /* we have to allocate memory first */
  if (DISTR.covar == NULL)
    DISTR.covar = _unur_malloc( dim * dim * sizeof(double) );
  if (DISTR.cholesky == NULL)
    DISTR.cholesky = _unur_malloc( dim * dim * sizeof(double) );   

  /* if covar == NULL --> use identity matrix */
  if (covar==NULL) { 
    for (i=0; i<dim; i++) { 
      for (j=0; j<dim; j++) {
         DISTR.covar[idx(i,j)] = (i==j) ? 1. : 0.;
         DISTR.cholesky[idx(i,j)] = (i==j) ? 1. : 0.;
      } 
    } 
  } 

  /* covariance matrix given --> copy data */
  else {
    
    /* check covariance matrix: diagonal entries > 0 */
    for (i=0; i<dim*dim; i+= dim+1)
      if (covar[i] <= 0.) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"variance <= 0");
	return UNUR_ERR_DISTR_DOMAIN;
      }

    /* check for symmetry */
    for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++)
	if (!_unur_FP_same(covar[i*dim+j],covar[j*dim+i])) {
	  _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,
	              "covariance matrix not symmetric");
	  return UNUR_ERR_DISTR_DOMAIN;
	}

    /* copy data */
    memcpy( DISTR.covar, covar, dim * dim * sizeof(double) );

    /* compute Cholesky decomposition and check for positive definitness */
    if (_unur_matrix_cholesky_decomposition(dim, covar, DISTR.cholesky) != UNUR_SUCCESS) {
      _unur_error(distr->name, UNUR_ERR_DISTR_DOMAIN, 
		  "covariance matrix not positive definite");
      return UNUR_ERR_DISTR_DOMAIN;      
    }

  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_COVAR | UNUR_DISTR_SET_CHOLESKY;

  /* o.k. */
  return UNUR_SUCCESS;

#undef idx
} /* end of unur_distr_cvec_set_covar() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_covar( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get covariance matrix of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to covariance matrix of distribution                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* covariance matrix known ? */
  if ( !(distr->set & UNUR_DISTR_SET_COVAR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix");
    return NULL;
  }

  return DISTR.covar;

} /* end of unur_distr_cvec_get_covar() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_cholesky( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get cholesky factor of the covariance matrix of distribution         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to cholesky factor                                         */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* covariance matrix known ? */
  if ( !(distr->set & UNUR_DISTR_SET_CHOLESKY) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix");
    return NULL;
  }

  return DISTR.cholesky;

} /* end of unur_distr_cvec_get_covar_cholesky() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_covar_inv ( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get inverse covariance matrix of distribution                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to inverse of covariance matrix                            */
     /*----------------------------------------------------------------------*/
{
  double det; /* determinant of covariance matrix */
  int dim;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  dim = distr->dim;

  /* covariance matrix known ? */
  if ( !(distr->set & UNUR_DISTR_SET_COVAR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"covariance matrix not known");
    return NULL;
  }

  /* allocate memory */
  if (DISTR.covar_inv == NULL)
    DISTR.covar_inv = _unur_malloc( dim * dim * sizeof(double) );   

  /* calculate inverse covariance matrix */
  if ( !(distr->set & UNUR_DISTR_SET_COVAR_INV) )
    if (_unur_matrix_invert_matrix(dim, DISTR.covar, COVARIANCE_DETMIN, DISTR.covar_inv, &det) != UNUR_SUCCESS) {
      _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"cannot compute inverse of covariance");
      return NULL;
    }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_COVAR_INV;

  return DISTR.covar_inv;

} /* end of unur_distr_cvec_get_covar_inv() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_rankcorr( struct unur_distr *distr, const double *rankcorr )
     /*----------------------------------------------------------------------*/
     /* Set rank-correlation matrix of distribution.                         */
     /* The given matrix is checked for symmetry and positive definitness.   */
     /* The diagonal entries must be equal to 1.                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   rankcorr ... rankcorrelation matrix of distribution                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int i,j;
  int dim;
  double *cholesky;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  dim = distr->dim;

  /* mark as unknown */
  distr->set &= ~(UNUR_DISTR_SET_RANKCORR);

  /* we have to allocate memory first */
  if (DISTR.rankcorr == NULL)
    DISTR.rankcorr = _unur_malloc( dim * dim * sizeof(double) );

  /* if rankcorr == NULL --> use identity matrix */
  if (rankcorr==NULL) { 
    for (i=0; i<dim; i++)
      for (j=0; j<dim; j++)
         DISTR.rankcorr[idx(i,j)] = (i==j) ? 1. : 0.;
  } 

  /* rankcorriance matrix given --> copy data */
  else {
    
    /* check rankcorriance matrix: diagonal entries == 1 */
    for (i=0; i<dim*dim; i+= dim+1) {
      if (!_unur_FP_same(rankcorr[i],1)) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"diagonals != 1");
	return UNUR_ERR_DISTR_DOMAIN;
      }
    }

    /* check for symmetry */
    for (i=0; i<dim; i++)
      for (j=i+1; j<dim; j++)
	if (!_unur_FP_same(rankcorr[i*dim+j],rankcorr[j*dim+i])) {
	  _unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,
	              "rank-correlation matrix not symmetric");
	  return UNUR_ERR_DISTR_DOMAIN;
	}

    /* copy data */
    memcpy( DISTR.rankcorr, rankcorr, dim * dim * sizeof(double) );

    /* compute Cholesky decomposition and check for positive definitness */
    cholesky = _unur_malloc( dim * dim * sizeof(double) );   
    if (_unur_matrix_cholesky_decomposition(dim, rankcorr, cholesky) != UNUR_SUCCESS) {
      _unur_error(distr->name, UNUR_ERR_DISTR_DOMAIN, 
		  "rankcorriance matrix not positive definite");
      free(cholesky);
      return UNUR_ERR_DISTR_DOMAIN;      
    }
    free(cholesky);
    
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_RANKCORR;

  /* o.k. */
  return UNUR_SUCCESS;

#undef idx
} /* end of unur_distr_cvec_set_rankcorr() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_rankcorr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get rank-correlation matrix of distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to rank-correlation matrix of distribution                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* rankcorriance matrix known ? */
  if ( !(distr->set & UNUR_DISTR_SET_RANKCORR) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"rank-correlation matrix");
    return NULL;
  }

  return DISTR.rankcorr;

} /* end of unur_distr_cvec_get_rankcorr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_marginals ( struct unur_distr *distr, struct unur_distr *marginal)
     /*----------------------------------------------------------------------*/
     /* Copy marinal distribution into distribution object.                  */
     /* Only one local copy is made.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   marginal ... pointer to marginal object                            */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *clone;
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, marginal, UNUR_ERR_NULL );
  _unur_check_distr_object( marginal, CONT, UNUR_ERR_DISTR_INVALID );

  /* make copy of marginal distribution object */
  clone = _unur_distr_clone( marginal );

  /* allocate memory for array */
  DISTR.marginals = _unur_malloc (distr->dim * sizeof(struct unur_distr *));

  /* copy pointer */
  for (i=0; i<distr->dim; i++)
    DISTR.marginals[i] = clone;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MARGINAL;

  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_marginals() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_marginal_array ( struct unur_distr *distr, struct unur_distr **marginals)
     /*----------------------------------------------------------------------*/
     /* Copy marinal distributions into distribution object.                 */
     /* For each dimension a new copy is made even if the pointer in         */
     /* the array marginals coincide.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr     ... pointer to distribution object                       */
     /*   marginals ... pointer to array of marginal objects                 */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, marginals, UNUR_ERR_NULL );

  for (i=0; i<distr->dim; i++) {
    _unur_check_NULL( distr->name, *(marginals+i), UNUR_ERR_NULL );
    _unur_check_distr_object( *(marginals+i), CONT, UNUR_ERR_DISTR_INVALID );
  }
    
  /* allocate memory for array */
  DISTR.marginals = _unur_malloc (distr->dim * sizeof(struct unur_distr *));

  /* make copy of marginal distribution object */
  for (i=0; i<distr->dim; i++) 
    DISTR.marginals[i] = _unur_distr_clone( *(marginals+i) );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MARGINAL;

  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_marginal_array() */

/*---------------------------------------------------------------------------*/

const struct unur_distr *
unur_distr_cvec_get_marginal( const struct unur_distr *distr, int n )
     /*----------------------------------------------------------------------*/
     /* Copy marinal distribution into distribution object.                  */
     /* Only one local copy is made.                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   marginal ... pointer to marginal object                            */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  if (n<=0 || n > distr->dim) {
    _unur_error(NULL,UNUR_ERR_DISTR_GET,"n not in 1 .. dim");
    return NULL;
  }

  /* mean vector known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MARGINAL) ) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_GET,"marginals");
    return NULL;
  }

  _unur_check_NULL( distr->name, DISTR.marginals, NULL );

  /* return marginal distribution */
  return (DISTR.marginals[n-1]);
} /* end of unur_distr_cvec_get_marginal() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_pdfparams( struct unur_distr *distr, int par, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set parameters for distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   par      ... which parameter is set                                */
     /*   params   ... parameter array with number `par'                     */
     /*   n_params ... length of parameter array                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( NULL, params, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  /* allocate memory */
  _unur_realloc( DISTR.params[par], n_params * sizeof(double) );

  /* copy parameters */
  memcpy( DISTR.params[par], params, n_params*sizeof(double) );

  /* set length of array */
  DISTR.n_params[par] = n_params;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_get_pdfparams( const struct unur_distr *distr, int par, const double **params )
     /*----------------------------------------------------------------------*/
     /* get number of PDF parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   par      ... which parameter is read                               */
     /*   params   ... pointer to parameter array with number `par'          */
     /*                                                                      */
     /* return:                                                              */
     /*   length of parameter array with number `par'                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  /* check new parameter for distribution */
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    *params = NULL;
    return 0;
  }
  
  *params = DISTR.params[par];

  return (*params) ? DISTR.n_params[par] : 0;
} /* end of unur_distr_cvec_get_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_mode( struct unur_distr *distr, const double *mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of PDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( distr->name, mode, UNUR_ERR_NULL );

  /* mode already set ? */
  if (DISTR.mode == NULL) {
    /* we have to allocate memory first */
    DISTR.mode = _unur_malloc( distr->dim * sizeof(double) );
  }

  /* copy data */
  memcpy( DISTR.mode, mode, distr->dim * sizeof(double) );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cvec_set_mode() */

/*---------------------------------------------------------------------------*/

const double *
unur_distr_cvec_get_mode( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to mode of distribution                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
    return NULL;
  }

  return DISTR.mode;

} /* end of unur_distr_cvec_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_pdfvol( struct unur_distr *distr, double volume )
     /*----------------------------------------------------------------------*/
     /* set volume below PDF                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   volume ... volume below PDF                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CVEC, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (volume <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"PDF volume <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  DISTR.volume = volume;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFVOLUME;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cvec_set_pdfvol() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cvec_get_pdfvol( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get volume below PDF of distribution                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   volume below PDF of distribution                                   */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CVEC, INFINITY );

  /* volume known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PDFVOLUME) ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"volume");
    return INFINITY;
  }

  return DISTR.volume;

} /* end of unur_distr_cvec_get_pdfvol() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cvec_debug( const struct unur_distr *distr, const char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  double *mat;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_CVEC,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous multivariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tdimension = %d\n",genid,distr->dim);
  fprintf(log,"%s:\n",genid);

  /* mean vector */
  mat = ((distr->set & UNUR_DISTR_SET_MEAN) && DISTR.mean) ? DISTR.mean : NULL;
  _unur_matrix_print_vector( distr->dim, mat, "\tmean vector =", log, genid, "\t   ");

  /* covariance matrix */
  mat = ((distr->set & UNUR_DISTR_SET_COVAR) && DISTR.covar) ? DISTR.covar : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\tcovariance matrix =", log, genid, "\t   ");

  /* cholesky matrix */
  mat = ((distr->set & UNUR_DISTR_SET_CHOLESKY) && DISTR.cholesky) ? DISTR.cholesky : NULL;
  _unur_matrix_print_matrix( distr->dim, mat, "\tcholesky factor (of covariance matrix) =", log, genid, "\t   ");

  /* marginal distributions */
  fprintf(log,"%s:\tmarginal distributions: ",genid);
  if (distr->set & UNUR_DISTR_SET_MARGINAL) {
    fprintf(log,"\n");
    if (DISTR.marginals[0] == DISTR.marginals[1]) {
      fprintf(log,"%s: all mariginals [1-%d]:\n",genid,distr->dim);
      _unur_distr_cont_debug( DISTR.marginals[0], genid );
    }
    else {
      int i;
      for (i=0; i<distr->dim; i++) {
	fprintf(log,"%s: mariginal [%d]:\n",genid,i+1);
	_unur_distr_cont_debug( DISTR.marginals[i], genid );
      }
    }
  }
  else {
    fprintf(log," [unknown]\n");
  }
  fprintf(log,"%s:\n",genid);

} /* end of _unur_distr_cvec_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

