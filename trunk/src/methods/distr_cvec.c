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

/*****************************************************************************/
/**                                                                         **/
/** mulitvariate continuous distributions                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cvec_new(  int dim )
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

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.dpdf      = NULL;          /* pointer to gradient of p.d.f.          */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pdf        */
  /* initialize parameters of the p.d.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = NULL;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.d.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode       = NULL;         /* location of mode (default: not known)  */
  DISTR.volume     = INFINITY;     /* area below p.d.f. (default: not known) */

  DISTR.upd_mode   = NULL;         /* funct for computing mode               */
  DISTR.upd_volume = NULL;         /* funct for computing area               */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cvec_new() */

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
unur_distr_cvec_set_pdfparams( struct unur_distr *distr, int par, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set parameters for distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   par      ... which parameter is set                                */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
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
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number of pdf parameters                                           */
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

  return (*params) ? 1 : 0;
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

int 
unur_distr_cvec_upd_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute mode of distribution (if possible)                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  if (DISTR.upd_mode == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  /* compute mode */
  if ((DISTR.upd_mode)(distr)) {
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE;
    return 1;
  }
  else {
    /* computing of mode failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"");
    return 0;
  }

} /* end of unur_distr_cvec_upd_mode() */
  
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
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return NULL;
    }
    else {
      /* compute mode */
      unur_distr_cvec_upd_mode( distr );
    }
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

int 
unur_distr_cvec_upd_pdfvol( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute volume below p.d.f. of distribution (if possible)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CVEC, 0 );

  if (DISTR.upd_volume == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  /* compute volume */
  if (!(DISTR.upd_volume)(distr)) {
    /* computing of volume failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"");
    return 0;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFVOLUME;

  return 1;
} /* end of unur_distr_cvec_upd_pdfvol() */
  
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
    /* try to compute volume */
    if (DISTR.upd_volume == NULL) {
      /* no function to compute volume available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"volume");
      return INFINITY;
    }
    else {
      /* compute volume */
      unur_distr_cvec_upd_pdfvol( distr );
    }
  }

  return DISTR.volume;

} /* end of unur_distr_cvec_get_pdfvol() */

/*---------------------------------------------------------------------------*/

#if 0
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
  int i;

  /* check arguments */
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_CVEC,/*void*/);

  log = unur_get_stream();

  /* is this a derived distribution */
  if (distr->base) {
    switch (distr->id) {
    case UNUR_DISTR_CORDER:
      _unur_distr_corder_debug(distr,genid);
      return;
    default:
      _unur_warning(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return;
    }
  }

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tp.d.f with %d argument(s)\n",genid,DISTR.n_params);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(log,"%s:\tmode = %g\n",genid,DISTR.mode);
  else
    fprintf(log,"%s:\tmode unknown\n",genid);

  fprintf(log,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]);
  _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);

  fprintf(log,"\n%s:\tarea below p.d.f. = %g",genid,DISTR.area);
  _unur_print_if_default(distr,UNUR_DISTR_SET_PDFAREA);
  fprintf(log,"\n%s:\n",genid);

} /* end of _unur_distr_cvec_debug() */

#endif 
/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

