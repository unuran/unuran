/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr_cvec.c                                                 *
 *                                                                           *
 *   manipulate multivariate continuous distribution objects                 *
 *                                                                           *
 *   return:                                                                 *
 *     1 ... on success                                                      *
 *     0 ... on error                                                        *
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

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

static const char unknown_distr_name[] = "unknown";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cvec

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
     /* type: multivariate continuous with given p.d.f.                      */
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

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cvec_free;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.dpdf      = NULL;          /* pointer to gradient of p.d.f.          */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  /* initialize parameters of the p.d.f.                                     */
  DISTR.mean  = NULL;              /* default is zero vector                 */
  DISTR.covar = NULL;              /* default is identity matrix             */

  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++) {
    DISTR.n_params[i] = 0;
    DISTR.params[i] = NULL;
  }

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.d.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode       = NULL;         /* location of mode (default: not known)  */
  DISTR.volume     = INFINITY;     /* area below p.d.f. (default: not known) */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cvec_new() */

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

  COOKIE_CHECK(distr,CK_DISTR_CVEC,/*void*/);

  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    if (DISTR.params[i]) free( DISTR.params[i] );

  if (DISTR.mean)  free(DISTR.mean); 
  if (DISTR.covar) free(DISTR.covar);

  if (DISTR.mode && DISTR.mode != DISTR.mean)
    /* only free mode if it does not point to the mean vector */
    free(DISTR.mode);

  free( distr );

} /* end of unur_distr_cvec_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CVEC *pdf )
     /*----------------------------------------------------------------------*/
     /* set p.d.f. of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );
  _unur_check_NULL( distr->name,pdf,0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  /* we do not allow overwriting a pdf */
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of pdf not allowed");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return 1;

} /* end of unur_distr_cvec_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_dpdf( struct unur_distr *distr, UNUR_VFUNCT_CVEC *dpdf )
     /*----------------------------------------------------------------------*/
     /* set gradient of p.d.f. of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dpdf  ... pointer to gradient of p.d.f.                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );
  _unur_check_NULL( distr->name,dpdf,0 );
  _unur_check_distr_object( distr, CVEC, 0 );
  
  /* we do not allow overwriting a dpdf */
  if (DISTR.dpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dpdf not allowed");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return 1;
} /* end of unur_distr_cvec_set_dpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CVEC *
unur_distr_cvec_get_pdf( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to p.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to p.d.f.                                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  return DISTR.pdf;
} /* end of unur_distr_cvec_get_pdf() */

/*---------------------------------------------------------------------------*/

UNUR_VFUNCT_CVEC *
unur_distr_cvec_get_dpdf( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to gradient of p.d.f. of distribution                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to gradient of p.d.f.                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, CVEC, NULL );

  return DISTR.dpdf;
} /* end of unur_distr_cvec_get_dpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cvec_eval_pdf( double *x, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate p.d.f. of distribution at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for pdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pdf(x)                                                             */
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
unur_distr_cvec_eval_dpdf( double *result, double *x, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate gradient of p.d.f. of distribution at x                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   result ... to store grad (pdf(x))                                  */
     /*   x      ... argument for dpdf                                       */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  return _unur_cvec_dPDF(result,x,distr);
} /* end of unur_distr_cvec_eval_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_set_mean( struct unur_distr *distr, double *mean )
     /*----------------------------------------------------------------------*/
     /* set mean vector of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mean  ... mean vector of distribution                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  if (mean) {
    if (DISTR.mean == NULL)
      /* we have to allocate memory first */
      DISTR.mean = _unur_malloc( distr->dim * sizeof(double) );
    /* copy data */
    memcpy( DISTR.mean, mean, distr->dim * sizeof(double) );
  }

  else { /* mean == NULL */
    /* we do not need the allocated memory any more */
    if (DISTR.mean) free (DISTR.mean);
  }

  if ( (distr->set & UNUR_DISTR_SET_MODE) &&
       DISTR.mode == DISTR.mean) 
    /* the mode is equal to the mean vector */
    DISTR.mode = DISTR.mean;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MEAN;

  /* o.k. */
  return 1;
} /* end of unur_distr_cvec_set_mean() */

/*---------------------------------------------------------------------------*/

double *
unur_distr_cvec_get_mean( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mean vector of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to mean of distribution                                    */
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
unur_distr_cvec_set_covar( struct unur_distr *distr, double *covar )
     /*----------------------------------------------------------------------*/
     /* set covariance matrix of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   covar ... covariance matrix of distribution                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  int i,j;
  int dim;

  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  dim = distr->dim;

  /* check covariance matrix: diagonal entries > 0 */
  if (covar) {
    for (i=0; i<dim*dim; i+= dim+1)
      if (covar[i] <= 0.) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"variance <= 0");
	return 0;
      }
  }

  /* check for symmetry */
  for (i=0; i<dim; i++)
    for (j=i+1; j<dim; j++)
      if (!_unur_FP_equal(covar[i*dim+j],covar[j*dim+i])) {
	_unur_error(distr->name ,UNUR_ERR_DISTR_DOMAIN,"covariance matrix not symmetric");
	return 0;
      }
	
  /* there is no check for positive definitness yet.        */
  /* this can be done easily during cholesky decomposition. */

  if (covar) {
    if (DISTR.covar == NULL)
      /* we have to allocate memory first */
      DISTR.covar = _unur_malloc( dim * dim * sizeof(double) );
    /* copy data */
    memcpy( DISTR.covar, covar, dim * dim * sizeof(double) );
  }

  else { /* covar == NULL */
    /* we do not need the allocated memory any more */
    if (DISTR.covar) free (DISTR.covar);
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_COVAR;

  /* o.k. */
  return 1;
} /* end of unur_distr_cvec_set_covar() */

/*---------------------------------------------------------------------------*/

double *
unur_distr_cvec_get_covar( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get covariance matrix of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to covariance matrix of distribution                       */
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

int
unur_distr_cvec_set_pdfparams( struct unur_distr *distr, int par, double *params, int n_params )
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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_NULL( NULL, params, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  /* check new parameter for distribution */
  if (par < 0 || par >= UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
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
  return 1;
} /* end of unur_distr_cvec_set_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cvec_get_pdfparams( struct unur_distr *distr, int par, double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pdf parameters and sets pointer to array params[] of   */
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
unur_distr_cvec_set_mode( struct unur_distr *distr, double *mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of p.d.f.                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  /* mode already set ? */
  if (DISTR.mode == NULL) {
    /* we have to allocate memory first */
    DISTR.mode = _unur_malloc( distr->dim );
  }

  /* copy data */
  memcpy( DISTR.mode, mode, distr->dim * sizeof(double) );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return 1;
} /* end of unur_distr_cvec_set_mode() */

/*---------------------------------------------------------------------------*/

double *
unur_distr_cvec_get_mode( struct unur_distr *distr )
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
     /* set volume below p.d.f.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   volume ... volume below p.d.f.                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  /* check new parameter for distribution */
  if (volume <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pdf volume <= 0");
    return 0;
  }

  DISTR.volume = volume;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFVOLUME;

  /* o.k. */
  return 1;

} /* end of unur_distr_cvec_set_pdfvol() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cvec_get_pdfvol( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get volume below p.d.f. of distribution                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   volume below p.d.f. of distribution                                */
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

void
_unur_distr_cvec_debug( struct unur_distr *distr, char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i,j;

  /* check arguments */
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_CVEC,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous multivariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tdimension = %d\n",genid,distr->dim);
  fprintf(log,"%s:\n",genid);

  fprintf(log,"%s:\tmean vector =\n",genid);
  fprintf(log,"%s:\t   ( %g",genid,(DISTR.mean ? DISTR.mean[0] : 0.));
  for (i=1; i<distr->dim; i++) 
    fprintf(log,", %g",(DISTR.mean ? DISTR.mean[i] : 0.));
  fprintf(log,")\t%s\n",(DISTR.mean ? "" : "[NULL]"));
  fprintf(log,"%s:\n",genid);
  
  fprintf(log,"%s:\tcovariance matrix = %s\n",genid,(DISTR.covar ? "" : "[NULL]"));
  for (j=0; j<distr->dim; j++) {
    fprintf(log,"%s:\t   (%7.4f",genid,(DISTR.covar ? DISTR.covar[distr->dim*j] : (j==0?1.:0.)));
    for (i=1; i<distr->dim; i++) 
      fprintf(log,",%7.4f",(DISTR.covar ? DISTR.covar[distr->dim*j+i] : (i==j?1.:0.)));
    fprintf(log,")\n");
  }
    
  fprintf(log,"%s:\n",genid);

} /* end of _unur_distr_cvec_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

