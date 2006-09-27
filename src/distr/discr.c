/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: discr.c                                                           *
 *                                                                           *
 *   manipulate univariate discrete distribution objects                     *
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
#include <parser/functparser_source.h>
#include "distr_source.h"
#include "distr.h"
#include "discr.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr

/*---------------------------------------------------------------------------*/

static double _unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for PMF.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for CDF.                                           */
/*---------------------------------------------------------------------------*/

static void _unur_distr_discr_free( struct unur_distr *distr );

/*---------------------------------------------------------------------------*/

static int _unur_distr_discr_find_mode( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* find mode of unimodal probability vector numerically by bisection.        */
/*---------------------------------------------------------------------------*/

inline double unur_distr_discr_eval_pv(int k, const UNUR_DISTR *distribution );
/*---------------------------------------------------------------------------*/
/* declare function inline.                                                  */
/*---------------------------------------------------------------------------*/


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

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_DISCR);

  /* set type of distribution */
  distr->type = UNUR_DISTR_DISCR;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* destructor */
  distr->destroy = _unur_distr_discr_free;

  /* clone */
  distr->clone = _unur_distr_discr_clone;

  /* set defaults                                                            */

  /* finite probability vector */
  DISTR.pv        = NULL;          /* probability vector (PV)                */
  DISTR.n_pv      = 0;             /* length of PV                           */

  /* probability mass function */
  DISTR.pmf       = NULL;          /* pointer to PMF                         */
  DISTR.cdf       = NULL;          /* pointer to CDF                         */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.set_params= NULL;          /* funct for setting parameters and domain*/

  DISTR.n_params  = 0;             /* number of parameters of the pmf        */
  /* initialize parameters of the PMF                                        */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for PMF
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.trunc[0] = DISTR.domain[0] = 0;         /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INT_MAX;   /* right boundary of domain  */

  DISTR.mode     = 0;              /* location of mode                       */
  DISTR.upd_mode = _unur_distr_discr_find_mode;  /* funct for computing mode */

  DISTR.sum     = 1.;              /* sum over PMF                           */
  DISTR.upd_sum = NULL;            /* funct for computing sum                */

  DISTR.pmftree    = NULL;         /* pointer to function tree for PMF       */
  DISTR.cdftree    = NULL;         /* pointer to function tree for CDF       */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_discr_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_discr_clone( const struct unur_distr *distr )
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
#define CLONE clone->data.discr

  struct unur_distr *clone;
  int len;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy function trees into generator object (when there is one) */
  CLONE.pmftree  = (DISTR.pmftree) ? _unur_fstr_dup_tree(DISTR.pmftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree) ? _unur_fstr_dup_tree(DISTR.cdftree) : NULL;

  /* copy probability vector into generator object (when there is one) */
  if (DISTR.pv) {
    CLONE.pv = _unur_xmalloc( DISTR.n_pv * sizeof(double) );
    memcpy( CLONE.pv, DISTR.pv, DISTR.n_pv * sizeof(double) );
  }

  /* copy user name for distribution */
  if (distr->name_str) {
    len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_discr_clone() */

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
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;
  _unur_check_distr_object( distr, DISCR, RETURN_VOID );

  if (DISTR.pmftree)  _unur_fstr_free(DISTR.pmftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);

  if (DISTR.pv) free( DISTR.pv );

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_discr_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pv( struct unur_distr *distr, const double *pv, int n_pv )
     /*----------------------------------------------------------------------*/
     /* set probability vector for distribution                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*   pv      ... pointer to PV                                          */
     /*   n_pv    ... length of PV                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* it is not possible to set a PV when a PMF is given. */
  if (DISTR.pmf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PMF given, cannot set PV");
    return UNUR_ERR_DISTR_SET;
  }

  /* check new parameter for distribution */
  if (n_pv < 0) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV");
    return UNUR_ERR_DISTR_SET;
  }

  /* n_pv must not be too large */
  if ( (DISTR.domain[0] > 0) && ((unsigned)DISTR.domain[0] + (unsigned)n_pv > INT_MAX) ) {
    /* n_pv too large, causes overflow */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV too large, overflow");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;

  /* we do not check non-negativity of p.v.
     (it is cheaper to do it when unur_init() is called */

  /* allocate memory for probability vector */
  DISTR.pv = _unur_xmalloc( n_pv * sizeof(double) );
  if (!DISTR.pv) return UNUR_ERR_MALLOC;

  /* copy probability vector */
  memcpy( DISTR.pv, pv, n_pv * sizeof(double) );
  DISTR.n_pv = n_pv;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pv() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_make_pv( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* compute probability vector                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*                                                                      */
     /* return:                                                              */
     /*   length of probability vector                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  double *pv;          /* pointer to probability vector */
  int n_pv;            /* length of PV */
  double cdf;          /* cumulated sum of PV */
  double thresh_cdf;   /* threshold for truncating PV */
  int valid;           /* whether cumputed PV is valid */
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* PMF or CDF required */
  if ( DISTR.pmf == NULL && DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PMF or CDF");
    return 0;
  }

  /* if there exists a PV, it has to be removed */
  if (DISTR.pv != NULL)
    free(DISTR.pv);

  /* compute PV */

  if ((unsigned)DISTR.domain[1] - (unsigned)DISTR.domain[0] < UNUR_MAX_AUTO_PV ) {

    /* first case: bounded domain */
    n_pv = DISTR.domain[1] - DISTR.domain[0] + 1;
    pv = _unur_xmalloc( n_pv * sizeof(double) );
    if (DISTR.pmf) {
      for (i=0; i<n_pv; i++)
	pv[i] = _unur_discr_PMF(DISTR.domain[0]+i,distr);
    }
    else if (DISTR.pmf) {
      int cdf_old = 0.;
      int cdf;
      for (i=0; i<n_pv; i++) {
	cdf = _unur_discr_CDF(DISTR.domain[0]+i,distr);
	pv[i] = cdf - cdf_old;
	cdf_old = cdf;
      }
    }
    valid = TRUE;
  }

  else {
    /* second case: domain too big but sum over PMF given       */
    /* we chop off the trailing part of the distribution        */
#define MALLOC_SIZE 1000 /* allocate 1000 doubles at once       */

    int n_alloc;         /* number of doubles allocated         */
    int max_alloc;       /* maximal number of allocated doubles */
    int size_alloc;      /* size of allocated blocks            */

    /* get maximal size of PV */
    if ( (DISTR.domain[0] <= 0) || (INT_MAX - DISTR.domain[0] >= UNUR_MAX_AUTO_PV - 1) ) {
      /* we can have a PV of length UNUR_MAX_AUTO_PV */
      size_alloc = MALLOC_SIZE;
      max_alloc = UNUR_MAX_AUTO_PV;
    }
    else { /* length of PV must be shorter than UNUR_MAX_AUTO_PV */
      size_alloc = max_alloc = INT_MAX - DISTR.domain[0];
    }

    /* init counter */
    n_pv = 0;
    pv = NULL;
    valid = FALSE;  /* created PV is empty yet and not valid */
    cdf = 0.;       /* cumulated sum of PV */
    /* threshold for truncating PV */
    thresh_cdf = (distr->set & UNUR_DISTR_SET_PMFSUM) ? (1.-1.e-8)*DISTR.sum : INFINITY;

    /* compute PV */
    for (n_alloc = size_alloc; n_alloc <= max_alloc; n_alloc += size_alloc) {
      pv = _unur_xrealloc( pv, n_alloc * sizeof(double) );

      if (DISTR.pmf) {
	for (i=0; i<size_alloc; i++) {
	  cdf += pv[n_pv] = _unur_discr_PMF(DISTR.domain[0]+n_pv,distr);
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }
      else if (DISTR.cdf) {
	int cdf_old = 0.;
	for (i=0; i<size_alloc; i++) {
	  cdf = _unur_discr_CDF(DISTR.domain[0]+i,distr);
	  pv[n_pv] = cdf - cdf_old;
	  cdf_old = cdf;
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }	  
      if (cdf > thresh_cdf) break;
    }

    if (distr->set & UNUR_DISTR_SET_PMFSUM) {
      /* make a warning if computed PV might not be valid */
      if (valid != TRUE)
	/* not successful */
	_unur_warning(distr->name,UNUR_ERR_DISTR_GET,"PV truncated");
    }
    else { /* PMFSUM not known */
      /* assume we have the important part of distribution */
      valid = TRUE;
      DISTR.sum = cdf;
      distr->set |= UNUR_DISTR_SET_PMFSUM;
    }
    
#undef MALLOC_SIZE
  }

  /* store vector */
  DISTR.pv = pv;
  DISTR.n_pv = n_pv;
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;

  /* o.k. */
  return (valid) ? n_pv : -n_pv;
} /* end of unur_distr_discr_make_pv() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_get_pv( const struct unur_distr *distr, const double **pv )
     /*----------------------------------------------------------------------*/
     /* get length of probability vector and set pointer to probability      */
     /* vector                                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   pv       ... pointer to probability vector                         */
     /*                                                                      */
     /* return:                                                              */
     /*   length of probability vector                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  *pv = (DISTR.pv) ? DISTR.pv : NULL;
  return DISTR.n_pv;

} /* end of unur_distr_discr_get_pv() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pv( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* returns the value of the probability vector at k or, if there is no  */
     /* probability vector defined, evaluates  the pmf                       */
     /*                                                                      */
     /* parampeters:                                                         */
     /*  k     ... argument for probability vector of pmf                    */
     /*  distr ... pointer to distribution object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   pv[k] or pmf(k)                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.pv != NULL)
    /* use probability vector */
    return (DISTR.pv[k-DISTR.domain[0]]);

  if (DISTR.pmf != NULL)
    /* use PMF */
    return _unur_discr_PMF(k,distr);

  /* else: data missing */
  _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
  return INFINITY;

} /* end of unur_distr_discr_eval_pv() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmf( struct unur_distr *distr, UNUR_FUNCT_DISCR *pmf )
     /*----------------------------------------------------------------------*/
     /* set PMF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pmf   ... pointer to PMF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr,  UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pmf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* it is not possible to set a PMF when a PV is given. */
  if (DISTR.pv != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PV given, cannot set PMF");
    return UNUR_ERR_DISTR_SET;
  }

  /* we do not allow overwriting a PMF */
  if (DISTR.pmf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.pmf = pmf;
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_pmf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_cdf( struct unur_distr *distr, UNUR_FUNCT_DISCR *cdf )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   cdf   ... pointer to CDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, cdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_pmf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to PMF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to PMF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.pmf;
} /* end of unur_distr_discr_get_pmf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_cdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to CDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to CDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.cdf;
} /* end of unur_distr_discr_get_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pmf( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PMF of distribution at k                                    */
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
unur_distr_discr_eval_cdf( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of distribution at k                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF(k)                                                             */
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
unur_distr_discr_set_pmfstr( struct unur_distr *distr, const char *pmfstr )
     /*----------------------------------------------------------------------*/
     /* set PMF of distribution via a string interface                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   pmfstr ... string that describes function term of PMF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pmfstr, UNUR_ERR_NULL );

  /* it is not possible to set a PMF when a PV is given. */
  if (DISTR.pv != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PV given, cannot set PMF");
    return UNUR_ERR_DISTR_SET;
  }

  /* we do not allow overwriting a PMF */
  if (DISTR.pmf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_DATA;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse PMF string */
  if ( (DISTR.pmftree = _unur_fstr2tree(pmfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pmf  = _unur_distr_discr_eval_pmf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pmfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution via a string interface                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   cdfstr ... string that describes function term of CDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_DATA;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL )
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;

  /* set evaluation function */
  DISTR.cdf  = _unur_distr_discr_eval_cdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_cdfstr() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for PMF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for PMF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PMF at k                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  return ((DISTR.pmftree) ? _unur_fstr_eval_tree(DISTR.pmftree,(double)k) : 0.);
} /* end of _unur_distr_discr_eval_pmf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for CDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at k                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,(double)k) : 0.);
} /* end of _unur_distr_discr_eval_cdf_tree() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_discr_get_pmfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get PMF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.pmftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.pmftree,"x","PMF",TRUE);
} /* end of unur_distr_discr_get_pmfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_discr_get_cdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get CDF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} /* end of unur_distr_discr_get_cdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfparams( struct unur_distr *distr, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set array of parameters for distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name, params, UNUR_ERR_NULL);

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* even if the set routine fails, the derived parameters are
     marked as unknown. but this is o.k. since in this case something
     has been wrong. */

  /* use special routine for setting parameters
     (if there is one) */

  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));

  /* otherwise simply copy parameters */

  /* but first check number of new parameter for the distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_pmfparams( const struct unur_distr *distr, const double **params )
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
unur_distr_discr_set_domain( struct unur_distr *distr, int left, int right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (left >= right) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* store data */
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = (DISTR.pv == NULL) ? right : left+DISTR.n_pv-1;

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

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_domain( const struct unur_distr *distr, int *left, int *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity            */
     /*   if no boundaries have been set [INT_MIN, INT_MAX] is returned.     */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = INT_MIN;
  *right = INT_MAX;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_get_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_mode( struct unur_distr *distr, int mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of PMF                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_upd_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute mode of distribution (if possible)                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  if (DISTR.upd_mode == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute mode */
  if ((DISTR.upd_mode)(distr)==UNUR_SUCCESS) {
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE;
    return UNUR_SUCCESS;
  }
  else {
    /* computing of mode failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

} /* end of unur_distr_discr_upd_mode() */
  
/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_mode( struct unur_distr *distr )
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
  _unur_check_NULL( NULL, distr, INT_MAX );
  _unur_check_distr_object( distr, DISCR, INT_MAX );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INT_MAX;
    }
    else {
      /* compute mode */
      if (unur_distr_discr_upd_mode(distr)!=UNUR_SUCCESS) {
	/* finding mode not successfully */
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INT_MAX;
      }
    }
  }

  return DISTR.mode;

} /* end of unur_distr_discr_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfsum( struct unur_distr *distr, double sum )
     /*----------------------------------------------------------------------*/
     /* set sum over PMF                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   sum   ... sum over PMF                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (sum <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"pmf sum <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  DISTR.sum = sum;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PMFSUM;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_pmfsum() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_upd_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute sum over PMF of distribution (if possible)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_SET );

  if (DISTR.upd_sum == NULL) {
    /* no function to compute sum available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute sum */
  if ((DISTR.upd_sum)(distr)==UNUR_SUCCESS) {
    /* changelog */
    distr->set |= UNUR_DISTR_SET_PMFSUM;
    return UNUR_SUCCESS;
  }
  else
    return UNUR_ERR_DISTR_DATA;

} /* end of unur_distr_discr_upd_pmfsum() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_discr_get_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get sum over PMF of distribution                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   sum over PMF of distribution                                       */
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
      if ((DISTR.upd_sum)(distr)==UNUR_SUCCESS)
	/* changelog */
	distr->set |= UNUR_DISTR_SET_PMFSUM;
      else
	return INFINITY;
    }
  }

  return DISTR.sum;

} /* end of unur_distr_discr_get_pmfsum() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_discr_debug( const struct unur_distr *distr, const char *genid, int printvector )
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
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_DISCR,RETURN_VOID);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  if ( DISTR.pmf ) {
    /* have probability mass function */
    fprintf(log,"%s:\tPMF with %d argument(s)\n",genid,DISTR.n_params);
    for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);
  }

  if (DISTR.n_pv>0) {
    /* have probability vector */
    fprintf(log,"%s:\tprobability vector of length %d",genid,DISTR.n_pv);
    if (printvector) {
      for (i=0; i<DISTR.n_pv; i++) {
	if (i%10 == 0)
	  fprintf(log,"\n%s:\t",genid);
	fprintf(log,"  %.5f",DISTR.pv[i]);
      }
    }
    fprintf(log,"\n%s:\n",genid);
  }

  /* domain */
  if ( DISTR.pmf ) {
    /* have probability mass function */
    fprintf(log,"%s:\tdomain for pmf = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[1]);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(log,"\n%s:\n",genid);
  }

  if (DISTR.n_pv>0) {
    /* have probability vector */
    fprintf(log,"%s:\tdomain for pv = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[0]-1+DISTR.n_pv);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(log,"\n%s:\n",genid);
  }

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(log,"%s:\tmode = %d\n",genid,DISTR.mode);
  else
    fprintf(log,"%s:\tmode unknown\n",genid);
  
  if (distr->set & UNUR_DISTR_SET_PMFSUM)
    fprintf(log,"%s:\tsum over PMF = %g\n",genid,DISTR.sum);
  else
    fprintf(log,"%s:\tsum over PMF unknown\n",genid);
  fprintf(log,"%s:\n",genid);
  
} /* end of _unur_distr_discr_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_distr_discr_find_mode(struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /*  find mode of a probability vector by bisection                      */
     /*                                                                      */
     /*  Any two of the three points x[i] must always differ at least by one */
     /*  If no further point xnew can be fitted between, the mode is found   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define sgn(a)  ( (a) >= (0) ? ( (a==0)?(0):(1) ) : (-1) )
#define max_pos3(a,b,c) ( (a) >= (b) ? ( ((a) >= (c)) ? (0) : (2) ) :\
                                     ( ((b) >= (c)) ? (1) : (2) ) )

#define INT1         (1) 
#define INT2         (2) 
#define INT3         (3)
#define UNDEFINED    (0)               
#define X2_BORDER    (1)              
#define XNEW_BORDER  (2)             

  int bisect;                     /* for choosing an interval               */
  int interval;                   /* interval containing xnew               */
  int mode;                       /* mode                                   */
  int x[3], xnew;                 /* mode between x[0] and x[1]             */
  double fx[3], fxnew;            /* ... and the respective function values */
  int xtmp = INT_MAX;
  double fxtmp = FLT_MAX;    /** TODO: FLT_MAX must be much smaller than DBL_MAX **/

  const double r = (3.-sqrt(5.))/2.;       /* sectio aurea                  */


  /* check arguments */
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
 
  /* derive three distinct points */
  x[0] = DISTR.domain[0];
  x[1] =DISTR.domain[1];
  fx[0] = unur_distr_discr_eval_pv(x[0], distr);
  fx[1] = unur_distr_discr_eval_pv(x[1], distr);


  if ( x[0] == x[1] ){            /* domain contains only one point         */
    mode = x[0];
  }
  else if ( x[1] == x[0] + 1 ){   /* domain contains only two points        */
    mode = (fx[0] >= fx[1]) ? x[0] : x[1];
  }
  else{                           /* domain contains at least three points  */

    x[2]  = (int) (r*x[0] + (1-r)*x[1]);
    if ( x[2] == x[0] )
      x[2]++;
    if ( x[2] == x[1] )
      x[2]--;
    fx[2] = unur_distr_discr_eval_pv(x[2], distr);


    /* at least one of the x[i] should have a positive function value  */
    if (fx[0] == 0.0 && fx[1] == 0.0 ){
      int i=1;
      while (fx[2] == 0. && i < 100){
	x[2]  = (x[1]/100)*i + (x[0]/100)*(100-i); /* integers !!! */
        fx[2] = unur_distr_discr_eval_pv(x[2], distr);
        i++;
      } 
    }
    if (fx[2] == 0.){  /* no success */
      _unur_error(distr->name,UNUR_ERR_DISTR_DATA,
         "In find_mode(): no positive entry in PV found during 100 trials");
      return UNUR_ERR_DISTR_DATA;  
    }

    /* x[i] are now initialized -- at least one entry is > 0  
       and no two of the x[i] are identical                  */ 


    while (1) {

      /* terminating the program legally */
      if ( (x[2]-x[0]) == sgn(x[2]-x[0]) &&
	   (x[1]-x[2]) == sgn(x[1]-x[2])    ){
	mode = x[ max_pos3(fx[0], fx[1], fx[2]) ];
	break;   /* mode found */
      }


      /* find xnew not identical with any of the x[i] */ 
      xnew  = (int) (r*x[0] + (1-r)*x[2]);

      if ( xnew == x[0] ){
	  xnew += sgn(x[2]-x[0]);
	  if (xnew == x[2]){
	    xnew = x[2];
	    x[2] += sgn(x[2]-x[0]);  /* cant be = x[1] */
	    fx[2] = unur_distr_discr_eval_pv(x[2], distr);
	  }
      }
      if ( xnew == x[2] ){
	  xnew -= sgn(x[2]-x[0]);
	  if (xnew == x[0]){
	    xnew = x[2];
	    x[2] += sgn(x[2]-x[0]);  /* cant be = x[1] */
	    fx[2] = unur_distr_discr_eval_pv(x[2], distr);
	  }
      }

      fxnew = unur_distr_discr_eval_pv(xnew, distr);

      /* Information of point xnew isn't enough to
         refine interval containig the mode -- determine new xnew    */
      if ( _unur_FP_same(fx[2], fxnew ) &&
           (! x[0] > x[2]) && (! x[1] > x[2]) ){

	interval = -1;  /* should be impossible when entering switch */
	if ( abs(x[1]-x[2]) > 1 ){
	  xtmp = x[1]/2 + x[2]/2;
	  fxtmp = unur_distr_discr_eval_pv(xtmp, distr);
	  interval = UNDEFINED;
          if ( ! _unur_FP_same(fxtmp, fx[2]) )
	    interval = INT3;
	}
	if ( abs(xnew-x[0]) > 1 ){
	  xtmp = xnew/2 + x[0]/2;
	  fxtmp = unur_distr_discr_eval_pv(xtmp, distr);
	  interval = UNDEFINED;
          if ( ! _unur_FP_same(fxtmp, fx[2]) )
	    interval = INT1;
        }
	if ( abs(x[2]-xnew) > 1 ){
	  xtmp = x[2]/2 + xnew/2;
	  fxtmp = unur_distr_discr_eval_pv(xtmp, distr);
	  interval = UNDEFINED;
          if ( ! _unur_FP_same(fxtmp, fx[2]) )
	    interval = INT2;
	}

	switch ( interval ){
	case INT1:
	  xnew = xtmp; fxnew = fxtmp;
	  break;
	case INT2:
	  xnew = xtmp; fxnew = fxtmp;
	  break;
	case INT3:
	  xnew = x[2]; fxnew = fx[2];
	  x[2] = xtmp; fx[2] = fxtmp;
	  break;
	case UNDEFINED:
	  return UNUR_ERR_DISTR_DATA;  /* mode not found -- exit */
	default:
	  _unur_error(distr->name, UNUR_ERR_SHOULD_NOT_HAPPEN,"");
	  return UNUR_ERR_SHOULD_NOT_HAPPEN;
	} /* end of switch (interval) */ 

      }   /* flat region left */


      /* regular bisection */
      bisect = -1; /* should be impossibe lwhen entering switch */
      if (      fxnew > fx[0] && fxnew > fx[2] ){
	bisect = X2_BORDER;
      }
      else if ( fxnew > fx[1] && fxnew > fx[2] ){
	bisect = X2_BORDER;
      }
      else if ( fx[2] > fxnew && fx[2] > fx[1] ){
	bisect = XNEW_BORDER;
      }
      else if ( fx[2] > fxnew && fx[2] > fx[0] ){
	bisect = XNEW_BORDER;
      }
      else if ( fx[0] > fxnew ){
	bisect = X2_BORDER;
      }
      else if ( fx[1] > fx[2] ){
	bisect = XNEW_BORDER;
      }
      else if ( _unur_FP_same(fx[0], fxnew) && fxnew < fx[2] ){
	bisect = XNEW_BORDER;
      }
      else if ( _unur_FP_same(fx[0], fxnew) && fxnew > fx[2] ){
	bisect = X2_BORDER;
      }
      else if ( _unur_FP_same(fx[2], fx[1]) && fxnew < fx[2] ){
	bisect = XNEW_BORDER;
      }
      else if ( _unur_FP_same(fx[2], fx[1]) && fxnew < fx[2] ){
	bisect = X2_BORDER;
      }
      else {
       _unur_error(distr->name, UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      }
      
      switch ( bisect ) {
      case XNEW_BORDER:
	x[0] = x[1];  fx[0] = fx[1];
	x[1] = xnew;  fx[1] = fxnew;
	break;
      case X2_BORDER:
	x[1] = x[2];  fx[1] = fx[2];
	x[2] = xnew;  fx[2] = fxnew;
	break;
      default:
	_unur_error(distr->name, UNUR_ERR_SHOULD_NOT_HAPPEN,"");
	return UNUR_ERR_SHOULD_NOT_HAPPEN;
      } /* end of switch (bisect) */


    } /* while (1) end */

  }  /* else (at least 3 points) end */
  
  /* mode successfully computed */
  DISTR.mode = mode;
  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE; 
  
  /* o.k. */
  return UNUR_SUCCESS;
  
#undef sgn
#undef max_pos3
#undef INT1
#undef INT2
#undef INT3
#undef UNDEFINED
#undef X2_BORDER
#undef XNEW_BORDER
} /* end of _unur_distr_discr_find_mode() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
