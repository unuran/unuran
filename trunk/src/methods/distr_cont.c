/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr_cont.c                                                 *
 *                                                                           *
 *   manipulate univariate continuous distribution objects                   *
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

#define DISTR distr->data.cont

/* for derived distributions (e.g. order statistics):
   data of underlying distributions */
#define BASE  distr->base->data.cont

/*---------------------------------------------------------------------------*/

static void _unur_distr_cont_free( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* destroy distribution object.                                              */
/*---------------------------------------------------------------------------*/

static double fmaxbr( UNUR_FUNCT_CONT *f, UNUR_DISTR *distr, double a, double b, double tol );
/*---------------------------------------------------------------------------*/
/* find maximum of function with Brent's algorithm.                          */
/*---------------------------------------------------------------------------*/

static int find_bracket( UNUR_FUNCT_CONT *f, UNUR_DISTR *distr, double *a, double *b );
/*---------------------------------------------------------------------------*/
/* find bracket for Brent's algorithm.                                       */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous distributions                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cont_new( void )
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

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cont_free;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.dpdf      = NULL;          /* pointer to derivative of p.d.f.        */
  DISTR.cdf       = NULL;          /* pointer to c.d.f.                      */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pdf        */
  /* initialize parameters of the p.d.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.d.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode      = INFINITY;      /* location of mode (default: not known)  */
  DISTR.area      = INFINITY;      /* area below p.d.f. (default: not known) */

  DISTR.trunc[0] = DISTR.domain[0] = -INFINITY; /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INFINITY;  /* right boundary of domain  */

  DISTR.upd_mode  = NULL;          /* funct for computing mode               */
  DISTR.upd_area  = NULL;          /* funct for computing area               */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cont_new() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cont_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  if (distr) free( distr );
} /* end of unur_distr_cont_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CONT *pdf )
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

  /* we do not allow overwriting a pdf */
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of pdf not allowed");
    return 0;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return 0;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return 1;

} /* end of unur_distr_cont_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_dpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of p.d.f. of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dpdf   ... pointer to derivative of p.d.f.                         */
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
  
  /* we do not allow overwriting a dpdf */
  if (DISTR.dpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dpdf not allowed");
    return 0;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return 0;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return 1;
} /* end of unur_distr_cont_set_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdf( struct unur_distr *distr, UNUR_FUNCT_CONT *cdf )
     /*----------------------------------------------------------------------*/
     /* set p.d.f. of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   cdf   ... pointer to c.d.f.                                        */
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
  
  /* we do not allow overwriting a cdf */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of cdf not allowed");
    return 0;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return 0;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return 1;
} /* end of unur_distr_cont_set_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_pdf( struct unur_distr *distr )
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

UNUR_FUNCT_CONT *
unur_distr_cont_get_dpdf( struct unur_distr *distr )
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

UNUR_FUNCT_CONT *
unur_distr_cont_get_cdf( struct unur_distr *distr )
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

  /* check new parameter for distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* copy parameters */

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       the parameters for the underlying distributions are set */
    BASE.n_params = n_params;
    if (n_params) memcpy( BASE.params, params, n_params*sizeof(double) );
  }
  else {
    DISTR.n_params = n_params;
    if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  }

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
     /*   number of pdf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       the parameters for the underlying distributions are returned */
    *params = (BASE.n_params) ? BASE.params : NULL;
    return BASE.n_params;
  }
  else {
    *params = (DISTR.n_params) ? DISTR.params : NULL;
    return DISTR.n_params;
  }

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

  /* check new parameter for distribution */
  if (left >= right) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }
  if (left < DISTR.domain[0] || right > DISTR.domain[1]) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain exceeds old domain, not allowed");
    return 0;
  }

  /* store data */
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = right;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_DOMAIN;

  /* if distr is an object for a standard distribution, this   */
  /* not the original domain of it. (not a "standard domain")  */
  /* However, since we have changed the domain, we assume      */
  /* that this is not a truncated distribution.                */
  /* At last we have to mark all derived parameters as unknown */
  distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_TRUNCATED | 
		  UNUR_DISTR_SET_MASK_DERIVED );

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       we also set the domain for the underlying distribution */
    BASE.trunc[0] = BASE.domain[0] = left;
    BASE.trunc[1] = BASE.domain[1] = right;
  }

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
unur_distr_cont_get_truncated( struct unur_distr *distr, double *left, double *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the truncated domain of the        */
     /* distribution                                                         */
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
  *left  = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[0] : DISTR.domain[0];
  *right = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[1] : DISTR.domain[1];

  return 1;
} /* end of unur_distr_cont_get_truncated() */

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
      unur_distr_cont_upd_mode( distr );
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

  /* check new parameter for distribution */
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

  /* compute area */
  if (!((DISTR.upd_area)(distr)) || DISTR.area <= 0.) {
    /* computing of area failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"upd area <= 0");
    DISTR.area = 1.;   /* avoid possible floating point exceptions */
    distr->set &= ~UNUR_DISTR_SET_PDFAREA;
    return 0;
  }

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

  /* area known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PDFAREA) ) {
    /* try to compute area */
    if (DISTR.upd_area == NULL) {
      /* no function to compute area available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"area");
      return INFINITY;
    }
    else {
      /* compute area */
      unur_distr_cont_upd_pdfarea( distr );
    }
  }

  return DISTR.area;

} /* end of unur_distr_cont_get_pdfarea() */

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

} /* end of _unur_distr_cont_debug() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_distr_cont_find_mode(struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* find mode of univariate continuous distribution numerically          */
     /*                                                                      */
     /* This is achieved by the following steps:                             */
     /* -- Determine a interval containing the mode and a third              */
     /*    point within ( x0< x1 < x2 )                                      */
     /*    ++ Determine a region where to search the mode; this will be the  */
     /*       interval [mode-100, mode+100] if this is no contrdiction to    */
     /*       the given domain;  `mode' is the best known approx to the mode */
     /*    ++ Find a point in the interval with a positve pdf                */
     /*       This is done by two geometric sequences of MAX_SRCH elements   */
     /*       and find two other points within the given domain              */
     /*    ++ Unbounded domains: refine x0, x1, x2 until:                    */
     /*       pdf(x0) < pdf(x1) > pdf(x2)                                    */
     /* -- invoke a maximization-routine to determine the mode               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{

  #define MAX_SRCH (50)
 
  int i;

  double x[3];   /* mode (and x[2]) should be between x[0] and  x[2] ...*/   
  double fx[3];  /* ... and the respective funtion values */
  double mode;   /* (approximative) mode of the distribution  */
  double mode_l; /* lower bound for mode search */
  double mode_u; /* upper bound for mode search */
  double step;

  int unbound_left; 
  int unbound_right; 


  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;


  mode = unur_distr_cont_get_mode(distr); /* yet best guess for mode */


  /* determine where to look for the mode */
  
  /* unbounded domain */
  if ( _unur_FP_is_minus_infinity(DISTR.domain[0]) &&
       _unur_FP_is_infinity(DISTR.domain[1]) ){

    unbound_left = 1;
    unbound_right = 1;

    x[1]  = mode;
    fx[1] = DISTR.pdf(x[1], distr);    
    mode_l = mode - 100.0;
    mode_u = mode + 100.0;

  }
  /* domain unbounded on the right */
  else if ( ! _unur_FP_is_minus_infinity(DISTR.domain[0]) &&
              _unur_FP_is_infinity(DISTR.domain[1]) ){

    unbound_left = 0;
    unbound_right = 1;

    if ( mode >= DISTR.domain[0] ){
      x[1]  = mode;
      fx[1] = DISTR.pdf(x[1], distr);
      mode_l = DISTR.domain[0];
      mode_u = 2 * mode - DISTR.domain[0];
    }
    else{
      x[1] = DISTR.domain[0] + 100.0;
      fx[1] = DISTR.pdf(x[1], distr);
      mode_l = DISTR.domain[0];
      mode_u = x[1] + 100.0;
    }

  }
  /* domain unbounded on the left */
  else if ( _unur_FP_is_minus_infinity(DISTR.domain[0]) &&
          ! _unur_FP_is_infinity(DISTR.domain[1]) ){

    unbound_left = 1;
    unbound_right = 0;

    if ( mode <= DISTR.domain[1] ){
      x[1]  = mode;
      fx[1] = DISTR.pdf(x[1], distr);
      mode_l = DISTR.domain[1] - 2 * mode;
      mode_u = DISTR.domain[1];
    }
    else{
      x[1] = DISTR.domain[1] - 100.0;
      fx[1] = DISTR.pdf(x[1], distr);
      mode_l = x[1] - 100.0;
      mode_u = DISTR.domain[1];      
    }

  }
  /* domain is bounded */
  else {  

    unbound_left = 0;
    unbound_right = 0;

    if ( mode >= DISTR.domain[0] && mode <= DISTR.domain[1] ){
      x[1]  = mode;
      fx[1] = DISTR.pdf(x[1], distr);
    }
    else{
      x[1] = DISTR.domain[0]/2.0 + DISTR.domain[1]/2.0;
      fx[1] = DISTR.pdf(x[1], distr);
    }
    mode_l = DISTR.domain[0];
    mode_u = DISTR.domain[1];

  }

  mode = x[1];  /* not exact mode -- best guess */


  /* find point with pdf > 0.0 -- max MAX_SRCH trials */

  /* search on the left side */
  step = pow(x[1]-mode_l, 1.0/MAX_SRCH);
  i = 0;  
  while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
    x[1]  = mode - pow(step, i);
    fx[1] = DISTR.pdf(x[1], distr);
    i++;
  }
	 
  /* search on the right side */
  if( _unur_FP_same(0.0, fx[1]) ){
    step = pow(mode_u-x[1], 1.0/MAX_SRCH);
    i = 0;
    while (i <= MAX_SRCH && _unur_FP_same(0.0, fx[1]) ){
      x[1]  = mode + pow(step, i);
      fx[1] = DISTR.pdf(x[1], distr);
      i++;
    }
  }
  
  /* no success -- exit routine  */   
  if( _unur_FP_same(fx[1], 0.0) )
     return 0;  /* can't find mode in flat region  */

  /* x[1] has pdf > 0 or routines already terminated */ 



  /* determine 3 points in the given domain --
     at least one with pdf > 0                  */
  if ( unbound_left ){

    x[2] = x[1];       fx[2] = fx[1];
    x[1] = x[2] - 1.0; fx[1] = DISTR.pdf(x[1], distr);
    x[0] = x[2] - 2.0; fx[0] = DISTR.pdf(x[0], distr);

  }
  else if ( unbound_right ){

    x[0] = x[1];       fx[0] = fx[1];
    x[1] = x[0] + 1.0; fx[1] = DISTR.pdf(x[1], distr);
    x[2] = x[0] + 2.0; fx[2] = DISTR.pdf(x[2], distr);

  }
  else{      /* bounded */

    x[0] = DISTR.domain[0];  fx[0] = DISTR.pdf(x[0], distr);
    x[2] = DISTR.domain[1];  fx[2] = DISTR.pdf(x[2], distr);

    if ( _unur_FP_same(x[1], DISTR.domain[0])  ||
         _unur_FP_same(x[1], DISTR.domain[1]) ){
      x[1]  = DISTR.domain[0]/2.0 + DISTR.domain[1]/2.0; 
      fx[1] = DISTR.pdf(x[1], distr);
    }
    
  }
  /* points x[i] with their function values determined */


  /* find interval containing the mode of the distribution */

/*    printf("-- x0: %f, fx0: %e\n", x[0], fx[0]); */
/*    printf("-- x1: %f, fx1: %e\n", x[1], fx[1]); */
/*    printf("-- x2: %f, fx2: %e\n", x[2], fx[2]); */


  step = 1.0;
  if ( unbound_right ){
    while(fx[0] <= fx[1] && fx[1] <= fx[2]){ /* on the left side of the mode */

      step *= 2.0;
      x[0]  = x[1]; fx[0] = fx[1];
      x[1]  = x[2]; fx[1] = fx[2];
      x[2] += step; fx[2] = DISTR.pdf(x[2], distr);   
    }
  }

  step = 1.0;  /* reset step size */
  if ( unbound_left ){
    while(fx[0] >= fx[1] && fx[1] >= fx[2]){ /* on the right side of the mode */

      step *= 2.0;
      x[2]  = x[1]; fx[2] = fx[1];
      x[1]  = x[0]; fx[1] = fx[0];
      x[0] -= step; fx[0] = DISTR.pdf(x[0], distr);

    }
  }

  /* now: the mode is between x[0] and x[2]   */



/*    printf("x0: %f, fx0: %e\n", x[0], fx[0]); */
/*    printf("x1: %f, fx1: %e\n", x[1], fx[1]); */
/*    printf("x2: %f, fx2: %e\n", x[2], fx[2]); */


  // TODO: invoke optimization


  /* o.k. */
  return 1;
} /* end of _unur_distr_cont_find_mode() */

/*---------------------------------------------------------------------------*/

/*
 *****************************************************************************
 *	    		    C math library                                   *
 * function FMINBR - one-dimensional search for a function minimum           *
 *			  over the given range                               *
 *                                                                           *
 * Author: Oleg Keselyov.                                                    *
 *                                                                           *
 * modified by Josef Leydold (documentation unchanged)                       *
 *                                                                           *
 * Input                                                                     *
 *	double fminbr(a,b,f,tol)                                             *
 *	double a; 			Minimum will be seeked for over      *
 *	double b;  			a range [a,b], a being < b.          *
 *	double (*f)(double x);		Name of the function whose minimum   *
 *					will be seeked for                   *
 *	double tol;			Acceptable tolerance for the minimum *
 *					location. It have to be positive     *
 *					(e.g. may be specified as EPSILON)   *
 *                                                                           *
 * Output                                                                    *
 *	Fminbr returns an estimate for the minimum location with accuracy    *
 *	3*SQRT_EPSILON*abs(x) + tol.                                         *
 *	The function always obtains a local minimum which coincides with     *
 *	the global one only if a function under investigation being          *
 *	unimodular.                                                          *
 *	If a function being examined possesses no local minimum within       *
 *	the given range, Fminbr returns 'a' (if f(a) < f(b)), otherwise      *
 *	it returns the right range boundary value b.                         *
 *                                                                           *
 * Algorithm                                                                 *
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical    *
 *	computations. M., Mir, 1980, p.202 of the Russian edition            *
 *                                                                           *
 *	The function makes use of the "gold section" procedure combined with *
 *	the parabolic interpolation.                                         *
 *	At every step program operates three abscissae - x,v, and w.         *
 *	x - the last and the best approximation to the minimum location,     *
 *	    i.e. f(x) <= f(a) or/and f(x) <= f(b)                            *
 *	    (if the function f has a local minimum in (a,b), then the both   *
 *	    conditions are fulfiled after one or two steps).                 *
 *	v,w are previous approximations to the minimum location. They may    *
 *	coincide with a, b, or x (although the algorithm tries to make all   *
 *	u, v, and w distinct). Points x, v, and w are used to construct      *
 *	interpolating parabola whose minimum will be treated as a new        *
 *	approximation to the minimum location if the former falls within     *
 *	[a,b] and reduces the range enveloping minimum more efficient than   *
 *	the gold section procedure.                                          *
 *	When f(x) has a second derivative positive at the minimum location   *
 *	(not coinciding with a or b) the procedure converges superlinearly   *
 *	at a rate order about 1.324                                          *
 *                                                                           *
 *****************************************************************************/

#define SQRT_EPSILON FLT_MAX

double
fmaxbr(f_invest,distr,a,b,tol)          /* An estimate to the min location   */
     UNUR_FUNCT_CONT *f_invest;         /* Function under investigation      */
     UNUR_DISTR *distr;                 /* distribution object with params   */
     double a;                          /* Left border | of the range	     */
     double b;                          /* Right border| the min is seeked   */
     double tol;                        /* Acceptable tolerance              */
{
#define f(x) (-(*f_invest)((x),distr))  /* minimize negative function        */

  double x,v,w;                         /* Abscissae, descr. see above       */
  double fx;                            /* f(x)                              */
  double fv;                            /* f(v)                              */
  double fw;                            /* f(w)                              */
  const double r = (3.-sqrt(5.0))/2;    /* Gold section ratio                */

  /* check arguments */
  CHECK_NULL(distr,0.);
  CHECK_NULL(f_invest,0.);
  COOKIE_CHECK(distr,CK_DISTR_CONT,0.);
  if (tol <= 0. || b <= a) {
    _unur_error(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0.;
  }

  /** TODO: abbruch nach zu vielen schritten!! **/


  v = a + r*(b-a);  fv = f(v);          /* First step - always gold section  */
  x = v;  w = v;
  fx=fv;  fw=fv;

  for(;;) {             /* Main iteration loop	*/
    double range = b-a;                 /* Range over which the minimum is   */
                                        /* seeked for		             */
    double middle_range = (a+b)/2;
    double tol_act =                    /* Actual tolerance                  */
		SQRT_EPSILON*fabs(x) + tol/3;
    double new_step;                    /* Step at this iteration            */
       

    if( fabs(x-middle_range) + range/2 <= 2*tol_act )
      return x;                         /* Acceptable approx. is found       */

                                        /* Obtain the gold section step      */
    new_step = r * ( x<middle_range ? b-x : a-x );


                                 /* Decide if the interpolation can be tried */
    if( fabs(x-w) >= tol_act ) {        /* If x and w are distinct           */
                                        /* interpolatiom may be tried	     */
      register double p;                /* Interpolation step is calculated  */
      register double q;                /* as p/q; division operation        */
                                        /* is delayed until last moment      */
      register double t;

      t = (x-w) * (fx-fv);
      q = (x-v) * (fx-fw);
      p = (x-v)*q - (x-w)*t;
      q = 2*(q-t);

      if( q>(double)0 )                 /* q was calculated with the         */
	p = -p;                         /* opposite sign; make q positive    */
      else                              /* and assign possible minus to	p    */
	q = -q;

      if( fabs(p) < fabs(new_step*q) && /* If x + p/q falls in [a,b] not too */
	  p > q*(a-x+2*tol_act)      && /* close to a and b, and isn't too   */
	  p < q*(b-x-2*tol_act)  )      /* large, it is accepted.            */
	new_step = p/q;                 /* If p/q is too large then the	     */
                                        /* gold section procedure can reduce */
                                        /* [a,b] range to more extent.       */
    }

    if( fabs(new_step) < tol_act ) {    /* Adjust the step to be not less    */
      if( new_step > (double)0 )        /* than tolerance.                   */
	new_step = tol_act;
      else
	new_step = -tol_act;
    }
                                 /* Obtain the next approximation to min     */
    {                            /* and reduce the enveloping range.         */
      register double t = x + new_step;	/* Tentative point for the min       */
      register double ft = f(t);

      if( ft <= fx ) {                  /* t is a better approximation	     */
	if( t < x )   
	  b = x;                        /* Reduce the range so that          */
	else                            /* t would fall within it            */
	  a = x;
      
	v = w;  w = x;  x = t;          /* Assign the best approx to x	     */
	fv=fw;  fw=fx;  fx=ft;
      }
      else {                            /* x remains the better approx       */
	if( t < x )
	  a = t;                        /* Reduce the range enclosing x      */
	else
	  b = t;
      
        if( ft <= fw || w==x ) {
	  v = w;  w = t;
	  fv=fw;  fw=ft;
        }
        else if( ft<=fv || v==x || v==w ) {
	  v = t;
	  fv=ft;
        }
      }
    }                   /* ----- end-of-block ----- */
  }                /* ===== End of loop ===== */

#undef f
} /* end of fmaxbr() */

#undef SQRT_EPSILON

/*---------------------------------------------------------------------------*/

int
find_bracket( UNUR_FUNCT_CONT *f_invest, UNUR_DISTR *distr, double *a, double *b )
     /*----------------------------------------------------------------------*/
     /* find mode of univariate continuous distribution numerically          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   f_invest ... function under investigation                          */
     /*   distr    ... pointer to distribution object  with parameters       */
     /*   a        ... left border   | of the range the max                  */
     /*   b        ... right border  | is seeked in fmaxbr                   */
     /*   tol      ... acceptable tolerance                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
#define f(x) (-(*f_invest)((x),distr))  /* minimize negative function        */

  double l, m, r;            /* left, middle, and right point of bracket     */
  double fl, fm, fr;         /* function value at l, m, and r, resp.         */
  int bd_left, bd_right;     /* left|right boundary of bracket found         */

  /* setup */
  
  /* we use the domain boundaries and the stored mode as starting points */
  l = DISTR.domain[0];
  r = DISTR.domain[1];

  /* mark boundaries as (in)finite */
  bd_left  = (_unur_FP_is_minus_infinity(DISTR.domain[0])) ? FALSE : TRUE;
  bd_right = (_unur_FP_is_infinity(DISTR.domain[1]))       ? FALSE : TRUE;

  /* for a bounded domain there is nothing to do */
  if (bd_left && bd_right) {
    *a = l;  *b = r;  return 1;
  }

  /* domain not bounded --> search */
 
  /* neither l nor r must be infinite */
  if (!bd_left)  l = (r < -1.) ? r-1. : -1;  /* start with arbitrary point instead */
  if (!bd_right) r = (l >  1.) ? l+1. :  1;

  fl = f(l);  fr = f(r);

  /** TODO: continue ... **/

  /* o.k. */
  return 1;

#undef f
} /* end of find_bracket() */


/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/








