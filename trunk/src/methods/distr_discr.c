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

// static int _unur_distr_discr_find_mode( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* find mode of unimodal probability vector numerically by bisection         */
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

  /* finite probability vector */
  DISTR.prob      = NULL;          /* probability vector                     */
  DISTR.n_prob    = 0;             /* length of probability vector           */

  /* probability mass function */
  DISTR.pmf       = NULL;          /* pointer to PMF                         */
  DISTR.cdf       = NULL;          /* pointer to CDF                         */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pmf        */
  /* initialize parameters of the PMF                                        */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for PMF
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.trunc[0] = DISTR.domain[0] = INT_MIN;   /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INT_MAX;   /* right boundary of domain  */

  DISTR.mode     = 0;              /* location of mode                       */
  DISTR.upd_mode = NULL;           /* funct for computing mode               */

  DISTR.sum     = 1.;              /* sum over PMF                           */
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
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;

  COOKIE_CHECK(distr,CK_DISTR_DISCR,/*void*/);

  if (DISTR.prob) free( DISTR.prob );

  free( distr );

} /* end of unur_distr_discr_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_prob( struct unur_distr *distr, double *prob, int n_prob )
     /*----------------------------------------------------------------------*/
     /* set probability vector for distribution                              */
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
  _unur_check_distr_object( distr, DISCR, 0 );

  /* it is not possible to set a PV when a PMF is given. */
  if (DISTR.pmf != NULL) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"PMF given, cannot set PV");
    return 0;
  }

  /* check new parameter for distribution */
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

  /* o.k. */
  return 1;
} /* end of unur_distr_discr_set_prob() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_get_prob( struct unur_distr *distr, double **prob )
     /*----------------------------------------------------------------------*/
     /* get length of probability vector and set pointer to probability      */
     /* vector.                                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   prob     ... pointer to probability vector                         */
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

  *prob = (DISTR.prob) ? DISTR.prob : NULL;
  return DISTR.n_prob;

} /* end of unur_distr_discr_get_prob() */

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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL,distr,0 );
  _unur_check_NULL( distr->name,pmf,0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* it is not possible to set a PMF when a PV is given. */
  if (DISTR.prob != NULL) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"PV given, cannot set PMF");
    return 0;
  }

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
     /* set CDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   cdf   ... pointer to CDF                                           */
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
  
  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
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
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.pmf;
} /* end of unur_distr_discr_get_pmf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_cdf( struct unur_distr *distr )
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
  _unur_check_NULL( NULL,distr,NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.cdf;
} /* end of unur_distr_discr_get_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pmf( int k, struct unur_distr *distr )
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
unur_distr_discr_eval_cdf( int k, struct unur_distr *distr )
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
unur_distr_discr_set_domain( struct unur_distr *distr, int left, int right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

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

  /* o.k. */
  return 1;

} /* end of unur_distr_discr_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_domain( struct unur_distr *distr, int *left, int *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
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
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return 1;
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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return 1;
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
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

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
      unur_distr_discr_upd_mode( distr );
    }
  }

  return DISTR.mode;

} /* end of unur_distr_discr_get_mode() */

/*---------------------------------------------------------------------------*/

int 
_unur_distr_discr_find_mode(struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /*  find mode of a probability vector by bisection                      */
     /*                                                                      */
     /*  Any two of the three points x[i] must always differ at least by one */
     /*  If no further point xnew can be fitted between, the mode is found   */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{

  #define sgn(a)  ( (a) >= (0) ? ( (a==0)?(0):(1) ) : (-1) )
  #define max_pos3(a,b,c) ( (a) >= (b) ? ( ((a) >= (c)) ? (1) : (3) ) :\
                                     ( ((b) >= (c)) ? (2) : (3) ) )

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
  double fxtmp = FLT_MAX;

  const double r = (3.-sqrt(5.))/2.;       /* sectio aurea                  */


  /* check arguments */
  CHECK_NULL( distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );
  if (DISTR.prob == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,
           "probability vektor required for finding mode numerically");
    return 0;
  }

  mode = INT_MAX;

  /* derive three distinct points */

  x[0] = DISTR.domain[0];
  x[1] = DISTR.domain[0] + DISTR.n_prob - 1;
  if (x[1] < x[0])
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,
		"Overflow: x[1] > INT_MAX");
  fx[0] = DISTR.prob[x[0]];
  fx[1] = DISTR.prob[x[1]];
  
  if ( x[0] == x[1] ){            /* domain contains only one point         */
    mode = x[0];
  }
  else if ( x[1] - x[0] == 1 ){   /* domain contains only two points        */
    mode = (fx[0] >= fx[1]) ? x[0] : x[1];
  }
  else{                           /* domain contains at least three points  */

    x[2]  = (int) (r*x[0] + (1-r)*x[1]);
    if ( x[2] == x[0] )
      x[2]++;
    if ( x[2] == x[1] )
      x[2]--;
    fx[2] = DISTR.prob[x[2]];

    /* at least one of the x[i] should have a positive function value  */
    if (x[0] == 0.0 && x[1] == 0.0 ){
      int i=1;
      while (fx[2] == 0 && i < 100){
	x[2]  = (x[1]/100)*i + (x[0]/100)*(100-i); /* integers !!! */
        fx[2] = DISTR.prob[x[2]];
        i++;
      } 
    }
    if (fx[2] == 0){  /* no success */
      _unur_error(distr->name,UNUR_ERR_DISTR_DATA,
         "In find_mode(): no positive entry in probability vector found
          during 100 trials");
      unur_distr_discr_set_mode(distr, INT_MAX);
      return 0;  
    }

    /* x[i] are now initialized -- at least one entry is > 0  
       and no two of the x[i] are identical                  */ 

    while (1){

      /* terminating the program legally */
      if ( (x[2]-x[0]) == sgn(x[2]-x[0]) &&
	   (x[1]-x[2]) == sgn(x[1]-x[2])    ){
	mode = x[ max_pos3(fx[1], fx[2], fx[3]) ];
	break;   /* mode found */
      }


      /* find xnew not identical with any of the x[i] */ 
      xnew  = (int) (r*x[0] + (1-r)*x[2]);

      if ( xnew == x[0] ){
	  xnew += sgn(x[2]-x[0]);
	  if (xnew == x[2]){
	    xnew = x[2];
	    x[2] += sgn(x[2]-x[0]);  /* cant be = x[1] */
	    fx[2] = DISTR.prob[x[2]];
	  }
      }
      if ( xnew == x[2] ){
	  xnew -= sgn(x[2]-x[0]);
	  if (xnew == x[0]){
	    xnew = x[2];
	    x[2] += sgn(x[2]-x[0]);  /* cant be = x[1] */
	    fx[2] = DISTR.prob[x[2]];
	  }
      }

      fxnew = DISTR.prob[xnew];


      /* Information of point xnew isn't enough to
         refine interval containig the mode -- determine new xnew    */
      if ( _unur_FP_same(fx[2], fxnew ) &&
           (! x[0] > x[2]) && (! x[1] > x[2]) ){

	interval = -1;  /* should be impossible when entering switch */
	if ( abs(x[1]-x[2]) > 1 ){
	  xtmp = x[1]/2 + x[2]/2;
	  fxtmp = DISTR.prob[xtmp];
	  interval = UNDEFINED;
          if ( ! _unur_FP_same(fxtmp, fx[2]) )
	    interval = INT3;
	}
	if ( abs(xnew-x[0]) > 1 ){
	  xtmp = xnew/2 + x[0]/2;
	  fxtmp = DISTR.prob[xtmp];
	  interval = UNDEFINED;
          if ( ! _unur_FP_same(fxtmp, fx[2]) )
	    interval = INT1;
        }
	if ( abs(x[2]-xnew) > 1 ){
	  xtmp = x[2]/2 + xnew/2;
	  fxtmp = DISTR.prob[xtmp];
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
	  unur_distr_discr_set_mode(distr, INT_MAX);
	  return 0;  /* mode not found -- exit */
	  break;
	default:
	  _unur_error(distr->name, UNUR_ERR_SHOULD_NOT_HAPPEN,"");
	  unur_distr_discr_set_mode(distr, INT_MAX);
	  return 0;
	  break;
	} /* end of switch (interval) */ 

      }   /* flat region left */

      /* regular bisection */
      bisect = -1; /* should be impossibe lwhen entering switch */
      if ( fxnew > fx[0] && fxnew > fx[2] ){
	bisect = X2_BORDER;
      }
      else if ( fxnew > fx[1] && fxnew > fx[2] ){
	bisect = XNEW_BORDER;
      }
      else if ( fx[2] > fxnew && fx[2] > fx[1])
	bisect = XNEW_BORDER;
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
      else{
       _unur_error(distr->name, UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      }
      
      switch ( bisect ){
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
	unur_distr_discr_set_mode(distr, INT_MAX);
	return 0;
      } /* end of switch (bisect) */


    } /* while (1) end */

  }  /* else (at least 3 points) end */
     
    unur_distr_discr_set_mode(distr, mode);
    /* o.k. */
    return 1;

} /* end of _unur_distr_discr_find_mode() */


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
     /* (re-) compute sum over PMF of distribution (if possible)             */
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
      DISTR.sum = (DISTR.upd_sum)(distr);
      /* changelog */
      distr->set |= UNUR_DISTR_SET_PMFSUM;
    }
  }

  return DISTR.sum;

} /* end of unur_distr_discr_get_pmfsum() */

/*---------------------------------------------------------------------------*/




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

  if ( DISTR.pmf ) {
    /* have probability mass function */
    fprintf(log,"%s:\tPMF with %d argument(s)\n",genid,DISTR.n_params);
    for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);
  }

  if (DISTR.n_prob>0) {
    /* have probability vector */
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

  /* domain */
  if ( DISTR.pmf )
    /* have probability mass function */
    fprintf(log,"%s:\tdomain for pmf = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[1]);

  if (DISTR.n_prob>0)
    /* have probability vector */
    fprintf(log,"%s:\tdomain for pv = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[0]-1+DISTR.n_prob);

  _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
  fprintf(log,"\n%s:\n",genid);

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(log,"%s:\tmode = %d\n",genid,DISTR.mode);
  else
    fprintf(log,"%s:\tmode unknown\n",genid);
  
  fprintf(log,"\n%s:\tsum over PMF = %g",genid,DISTR.sum);
  _unur_print_if_default(distr,UNUR_DISTR_SET_PMFSUM);
  
} /* end of _unur_distr_discr_debug() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/









