/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr.c                                                      *
 *                                                                           *
 *   manipulate distribution obejct                                          *
 *                                                                           *
 *   PARAMETER: struct unur_distr *                                          *
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
#include <source_stddistr.h>

/*---------------------------------------------------------------------------*/

static const char unknown_distr_name[] = "unknown";

/*---------------------------------------------------------------------------*/

/* create empty distribution object for ...                                  */
inline struct unur_distr *_unur_distr_cont_new( void );  /* univ. continuous */
inline struct unur_distr *_unur_distr_discr_new( void ); /* univ. discrete   */
inline struct unur_distr *_unur_distr_demp_new( void );  /* emp. univ. discr.*/

inline int unur_distr_discr_set_prob( struct unur_distr *distr, double *prob, int n_prob );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** routines for all distribution objects                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_new( unsigned int type )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   type ... type of distribution                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  switch (type) {
  case UNUR_DISTR_CONT:
    return _unur_distr_cont_new();
  case UNUR_DISTR_DISCR:
    return _unur_distr_discr_new();
  case UNUR_DISTR_DEMP:
    return _unur_distr_demp_new();
  default:
    _unur_error(NULL,UNUR_ERR_DISTR_UNKNOWN,"");
    return NULL;
  }

} /* end of unur_distr_new() */

/*---------------------------------------------------------------------------*/

void
unur_distr_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,/*void*/);
    break;
  case UNUR_DISTR_DISCR:
    COOKIE_CHECK(distr,CK_DISTR_DISCR,/*void*/);
    break;
  case UNUR_DISTR_DEMP:
    COOKIE_CHECK(distr,CK_DISTR_DEMP,/*void*/);
    if (distr->data.demp.prob) free( distr->data.demp.prob );
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_DISTR_UNKNOWN,"");
  }

  free( distr );

} /* end of unur_distr_free() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_set_name( struct unur_distr *distr, const char *name )
     /*----------------------------------------------------------------------*/
     /* set name of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   name  ... name of distribution                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  distr->name = name;
  return 1;
} /* end of unur_distr_set_name() */

/*---------------------------------------------------------------------------*/

const char *
unur_distr_get_name( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get name of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   name of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return distr->name;
} /* end of unur_distr_get_name() */

/*---------------------------------------------------------------------------*/

unsigned int unur_distr_get_type( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get type of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   type of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return (distr->type);
} /* end of unur_distr_get_type() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_cont( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate continuous.                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate continuous                                     */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_CONT) ? 1 : 0);
} /* end of unur_distr_is_cont() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_is_discr( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* TRUE if distribution is univariate discrete.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if univariate discrete                                       */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );

  return ((distr->type == UNUR_DISTR_DISCR) ? 1 : 0);
} /* end of unur_distr_is_discr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous distributions                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define DISTR distr->data.cont
/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cont_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate continuous with given p.d.f.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
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

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.dpdf      = NULL;          /* pointer to derivative of p.d.f.        */
  DISTR.cdf       = NULL;          /* pointer to c.d.f.                      */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pdf        */
  /* initialize parameters of the p.d.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS+1; i++)
    DISTR.params[i] = 0.;

  DISTR.mode      = INFINITY;      /* location of mode (default: not known)  */
  DISTR.area      = INFINITY;      /* area below p.d.f. (default: not known) */
  DISTR.domain[0] = -INFINITY;     /* left boundary of domain                */
  DISTR.domain[1] = INFINITY;      /* right boundary of domain               */

  DISTR.upd_mode  = NULL;          /* funct for computing mode               */
  DISTR.upd_area  = NULL;          /* funct for computing area               */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of _unur_distr_cont_new() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdf( struct unur_distr *distr, void *pdf )
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
  _unur_check_distr_object( distr, CONT, 0 );

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return 1;

} /* end of unur_distr_cont_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_dpdf( struct unur_distr *distr, void *dpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of p.d.f. of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to derivative of p.d.f.                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );
  _unur_check_NULL( distr->name,dpdf,0 );
  _unur_check_distr_object( distr, CONT, 0 );
  
  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return 1;
} /* end of unur_distr_cont_set_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdf( struct unur_distr *distr, void *cdf )
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
  _unur_check_NULL( distr->name,cdf,0 );
  _unur_check_distr_object( distr, CONT, 0 );
  
  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return 1;
} /* end of unur_distr_cont_set_cdf() */

/*---------------------------------------------------------------------------*/

void *unur_distr_cont_get_pdf( struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.pdf;
} /* end of unur_distr_cont_get_pdf() */

/*---------------------------------------------------------------------------*/

void *unur_distr_cont_get_dpdf( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to derivative of p.d.f. of distribution                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to derivative of p.d.f.                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.dpdf;
} /* end of unur_distr_cont_get_dpdf() */

/*---------------------------------------------------------------------------*/

void *unur_distr_cont_get_cdf( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to c.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to c.d.f.                                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.cdf;
} /* end of unur_distr_cont_get_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_pdf( double x, struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.pdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_PDF(x,distr);
} /* end of unur_distr_cont_eval_pdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_dpdf( double x, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate derivative of p.d.f. of distribution at x                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for dpdf                                        */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   (pdf(x))'                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_dPDF(x,distr);
} /* end of unur_distr_cont_eval_dpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_cdf( double x, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate c.d.f. of distribution at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   cdf(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.cdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_CDF(x,distr);
} /* end of unur_distr_cont_eval_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfparams( struct unur_distr *distr, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set array of parameters for distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
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
  _unur_check_distr_object( distr, CONT, 0 );
  if (n_params>0) _unur_check_NULL(distr->name,params,0);

  /* check new parameter for generator */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* copy parameters */
  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );

  /* o.k. */
  return 1;
} /* end of unur_distr_cont_set_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_pdfparams( struct unur_distr *distr, double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pdf parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number if pdf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return -1                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, -1 );
  _unur_check_distr_object( distr, CONT, -1 );

  *params = (DISTR.n_params) ? DISTR.params : NULL;
  return DISTR.n_params;

} /* end of unur_distr_cont_get_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_domain( struct unur_distr *distr, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  DISTR.domain[0] = left;
  DISTR.domain[1] = right;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_DOMAIN;

  /* if distr is an object for a standard distribution, */
  /* we might have truncated the distribution!          */
  distr->set &= ~UNUR_DISTR_SET_STDDOMAIN;

  /* derived parameters like mode, area, etc. might be wrong now! */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;

  /* o.k. */
  return 1;

} /* end of unur_distr_cont_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_domain( struct unur_distr *distr, double *left, double *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   if no boundaries have been set +/- INFINITY is returned.           */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = -INFINITY;
  *right = INFINITY;

  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return 1;
} /* end of unur_distr_cont_get_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_mode( struct unur_distr *distr, double mode )
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
  _unur_check_distr_object( distr, CONT, 0 );

  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return 1;
} /* end of unur_distr_cont_set_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_mode( struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, 0 );

  if (DISTR.upd_mode == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  /* compute mode */
  DISTR.mode = (DISTR.upd_mode)(distr);

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  return 1;
} /* end of unur_distr_cont_upd_mode() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   mode of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INFINITY;
    }
    else {
      /* compute mode */
      DISTR.mode = (DISTR.upd_mode)(distr);
      /* changelog */
      distr->set |= UNUR_DISTR_SET_MODE;
    }
  }

  return DISTR.mode;

} /* end of unur_distr_cont_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfarea( struct unur_distr *distr, double area )
     /*----------------------------------------------------------------------*/
     /* set area below p.d.f.                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   area  ... area below p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pdf area <= 0");
    return 0;
  }

  DISTR.area = area;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  /* o.k. */
  return 1;

} /* end of unur_distr_cont_set_pdfarea() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute area below p.d.f. of distribution (if possible)        */
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
  _unur_check_distr_object( distr, CONT, 0 );

  if (DISTR.upd_area == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  /* compute mode */
  DISTR.area = (DISTR.upd_area)(distr);

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  return 1;
} /* end of unur_distr_cont_upd_pdfarea() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get area below p.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   area below p.d.f. of distribution                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PDFAREA) ) {
    /* try to compute area */
    if (DISTR.upd_area == NULL) {
      /* no function to compute area available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"area");
      return INFINITY;
    }
    else {
      /* compute mode */
      DISTR.area = (DISTR.upd_area)(distr);
      /* changelog */
      distr->set |= UNUR_DISTR_SET_PDFAREA;
    }
  }

  return DISTR.area;

} /* end of unur_distr_cont_get_pdfarea() */

/*---------------------------------------------------------------------------*/


/*****************************************************************************/

void
_unur_distr_cont_debug( struct unur_distr *distr, char *genid )
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
  COOKIE_CHECK(distr,CK_DISTR_CONT,/*void*/);

  log = unur_get_stream();

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

} /* end of _unur_distr_cont_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** empirical univariate discrete distributions                             **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define DISTR distr->data.demp
/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_demp_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate discete                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_DEMP);

  /* set type of distribution */
  distr->type = UNUR_DISTR_DEMP;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* set defaults                                                            */

  /* finite probability vector */
  DISTR.prob      = NULL;          /* probability vector                     */
  DISTR.n_prob    = 0;             /* length of probability vector           */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of _unur_distr_demp_new() */

/*---------------------------------------------------------------------------*/

int
unur_distr_demp_set_prob( struct unur_distr *distr, double *prob, int n_prob )
     /*----------------------------------------------------------------------*/
     /* set probability vector for distribution                              */
     /* set array of parameters for distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*   prob    ... pointer to probability vector                          */
     /*   n_prob  ... length of probability vector                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DEMP, 0 );

  /* check new parameter for generator */
  if (n_prob < 0) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"length of p.v.");
    return 0;
  }
  /* we do not check non-negativity of p.v.
     (it is cheaper to do it when unur_init() is called */

  /* allocate memory for probability vector */
  DISTR.prob = _unur_malloc( n_prob * sizeof(double) );
  if (!DISTR.prob) return 0;

  /* copy probability vector */
  memcpy( DISTR.prob, prob, n_prob * sizeof(double) );
  DISTR.n_prob = n_prob;

  /* set name for distribution */
  distr->name = "(empirical)";

  /* o.k. */
  return 1;
} /* end of unur_distr_demp_set_prob() */

/*****************************************************************************/

void
_unur_distr_demp_debug( struct unur_distr *distr, char *genid, int printvector )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*   printvector ... print probability vector if not 0                  */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_DEMP,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  if (DISTR.n_prob>0) {
    /* probability vector given */
    fprintf(log,"%s:\tprobability vector of length %d",genid,DISTR.n_prob);
    if (printvector) {
      for (i=0; i<DISTR.n_prob; i++) {
	if (i%10 == 0)
	  fprintf(log,"\n%s:\t",genid);
	fprintf(log,"  %.5f",DISTR.prob[i]);
      }
    }
    fprintf(log,"\n%s:\n",genid);
  }
} /* end of _unur_distr_demp_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate discrete distributions                                       **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define DISTR distr->data.discr
/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_discr_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate discete                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  register int i;

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_DISCR);

  /* set type of distribution */
  distr->type = UNUR_DISTR_DISCR;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* set defaults                                                            */

  /* finite probability vector */
  DISTR.prob      = NULL;          /* probability vector                     */
  DISTR.n_prob    = 0;             /* length of probability vector           */

  /* probability mass function */
  DISTR.pmf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.cdf       = NULL;          /* pointer to c.d.f.                      */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pmf        */
  /* initialize parameters of the p.m.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS+1; i++)
    DISTR.params[i] = 0.;

  /* DISTR.mode      = 0.;            location of mode                       */
  DISTR.area      = 1.;            /* area below p.m.f.                      */
  DISTR.domain[0] = 1.;            /* left boundary of domain                */
  DISTR.domain[1] = INFINITY;      /* right boundary of domain               */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of _unur_distr_discr_new() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

void
_unur_distr_discr_debug( struct unur_distr *distr, char *genid, int printvector )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*   printvector ... print probability vector if not 0                  */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_DISCR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  if (DISTR.n_prob>0) {
    /* probability vector given */
    fprintf(log,"%s:\tprobability vector of length %d",genid,DISTR.n_prob);
    if (printvector) {
      for (i=0; i<DISTR.n_prob; i++) {
	if (i%10 == 0)
	  fprintf(log,"\n%s:\t",genid);
	fprintf(log,"  %.5f",DISTR.prob[i]);
      }
    }
    fprintf(log,"\n%s:\n",genid);
  }
  else if (DISTR.pmf) {
    /* probability mass function given */
    fprintf(log,"%s:\tp.m.f with %d argument(s)\n",genid,DISTR.n_params);
    for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);

    /*      if (distr->set & UNUR_DISTR_SET_MODE) */
    /*        fprintf(log,"%s:\tmode = %g\n",genid,DISTR.mode); */
    /*      else */
    /*        fprintf(log,"%s:\tmode unknown\n",genid); */

    /*    fprintf(log,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]); */
    /*    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN); */

    /*      fprintf(log,"\n%s:\tarea below p.d.f. = %g",genid,DISTR.area); */
    /*      _unur_print_if_default(distr,UNUR_DISTR_SET_PDFAREA); */

    fprintf(log,"%s:\n",genid);
  }

} /* end of _unur_distr_discr_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/

