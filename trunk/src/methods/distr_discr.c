/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr_discr.c                                                *
 *                                                                           *
 *   manipulate univariate discrete distribution objects                     *
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

#define DISTR distr->data.discr

/*---------------------------------------------------------------------------*/

static void _unur_distr_discr_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate discrete distributions                                       **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_discr_new( void )
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

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_discr_free;

  /* set defaults                                                            */

  /* probability mass function */
  DISTR.pmf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.cdf       = NULL;          /* pointer to c.d.f.                      */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pmf        */
  /* initialize parameters of the p.m.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.m.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  /* DISTR.domain[0] = ?;             left boundary of domain                */
  /* DISTR.domain[1] = ?;             right boundary of domain               */

  /* DISTR.mode      = 0.;            location of mode                       */

  DISTR.sum     = 1.;              /* sum over p.m.f.                        */
  DISTR.upd_sum = NULL;            /* funct for computing sum                */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_discr_new() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_discr_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  if (distr) free( distr );
} /* end of unur_distr_discr_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmf( struct unur_distr *distr, UNUR_FUNCT_DISCR *pmf )
     /*----------------------------------------------------------------------*/
     /* set p.m.f. of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pmf   ... pointer to p.m.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );
  _unur_check_NULL( distr->name,pmf,0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* we do not allow overwriting a pmf */
  if (DISTR.pmf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of pmf not allowed");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.pmf = pmf;
  return 1;

} /* end of unur_distr_discr_set_pmf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_cdf( struct unur_distr *distr, UNUR_FUNCT_DISCR *cdf )
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
  _unur_check_distr_object( distr, DISCR, 0 );
  
  /* we do not allow overwriting a cdf */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of cdf not allowed");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return 1;
} /* end of unur_distr_discr_set_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_pmf( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to p.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to p.m.f.                                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.pmf;
} /* end of unur_distr_discr_get_pmf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_cdf( struct unur_distr *distr )
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
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.cdf;
} /* end of unur_distr_discr_get_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pmf( int k, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate p.m.f. of distribution at k                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for pmf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pmf(k)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.pmf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_discr_PMF(k,distr);
} /* end of unur_distr_discr_eval_pmf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_cdf( int k, struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate c.d.f. of distribution at k                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   cdf(k)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.cdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_discr_CDF(k,distr);
} /* end of unur_distr_discr_eval_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfparams( struct unur_distr *distr, double *params, int n_params )
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
  _unur_check_distr_object( distr, DISCR, 0 );
  if (n_params>0) _unur_check_NULL(distr->name,params,0);

  /* check new parameter for distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  /* copy parameters */
  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );

  /* o.k. */
  return 1;
} /* end of unur_distr_discr_set_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_pmfparams( struct unur_distr *distr, double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pmf parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number of pmf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  *params = (DISTR.n_params) ? DISTR.params : NULL;
  return DISTR.n_params;

} /* end of unur_distr_discr_get_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfsum( struct unur_distr *distr, double sum )
     /*----------------------------------------------------------------------*/
     /* set sum over p.m.f.                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   sum   ... sum over p.d.f.                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* check new parameter for distribution */
  if (sum <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pmf sum <= 0");
    return 0;
  }

  DISTR.sum = sum;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PMFSUM;

  /* o.k. */
  return 1;

} /* end of unur_distr_discr_set_pmfsum() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_upd_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute sum over p.m.f. of distribution (if possible)          */
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
  _unur_check_distr_object( distr, DISCR, 0 );

  if (DISTR.upd_sum == NULL) {
    /* no function to compute sum available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return 0;
  }

  /* compute sum */
  DISTR.sum = (DISTR.upd_sum)(distr);

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PMFSUM;

  return 1;
} /* end of unur_distr_discr_upd_pmfsum() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_discr_get_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get sum over p.m.f. of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   sum over p.m.f. of distribution                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PMFSUM) ) {
    /* try to compute sum */
    if (DISTR.upd_sum == NULL) {
      /* no function to compute sum available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"sum");
      return INFINITY;
    }
    else {
      /* compute sum */
      DISTR.sum = (DISTR.upd_sum)(distr);
      /* changelog */
      distr->set |= UNUR_DISTR_SET_PMFSUM;
    }
  }

  return DISTR.sum;

} /* end of unur_distr_discr_get_pmfsum() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

void
_unur_distr_discr_debug( struct unur_distr *distr, char *genid )
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
  COOKIE_CHECK(distr,CK_DISTR_DISCR,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tp.m.f with %d argument(s)\n",genid,DISTR.n_params);
  for( i=0; i<DISTR.n_params; i++ )
    fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);

  /*      if (distr->set & UNUR_DISTR_SET_MODE) */
  /*        fprintf(log,"%s:\tmode = %g\n",genid,DISTR.mode); */
  /*      else */
  /*        fprintf(log,"%s:\tmode unknown\n",genid); */
  
  /*    fprintf(log,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]); */
  /*    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN); */
  
  /*      fprintf(log,"\n%s:\tsum over p.m.f. = %g",genid,DISTR.sum); */
  /*      _unur_print_if_default(distr,UNUR_DISTR_SET_PMFSUM); */
  
  fprintf(log,"%s:\n",genid);

} /* end of _unur_distr_discr_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
