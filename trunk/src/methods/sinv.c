/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      sinv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    Spline approximation for INVerse of CDF                      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the CDF                                                   *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      pointer to PDF and dPDF                                              *
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
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Numerical inversion is a method for generating random variables          *
 *  using the CDF (and in case of newton's method the PDF).                  *
 *                                                                           *
 *  THEOREM:                                                                 *
 *     Let X be a random variable with CDF F(x).                             *
 *     Then the F(X) are  uniformly distributed.                             *
 *                                                                           *
 *  COROLLARY:                                                               *
 *     Starting with uniformly distributed random variables U,               *
 *     the F^(-1)(U) have F(x) as CDF.                                       *
 *                                                                           *
 *  Starting with an U, the task is to find a X fulfilling:                  *
 *    F(X) - U = 0.                                                          *
 *                                                                           *
 *  Numerical algorithms to find zeros that are used in SINV are variants of * 
 *  newton's method (damped newton to guarantee improvement) and             *
 *  the regula falsi ( stabilized regula falsi preserving sign change; at    *
 *  first an interval with sign change is determined).                       *
 *                                                                           *
 *  In both cases it is possible to specify the maximal number of            *
 *  iterations, a desired accuracy in X and starting values for the          *
 *  algorithms.                                                              *
 *  Instead of starting values it is also possible to use a table            *
 *  containing suitable starting values.                                     *
 *  If neither the table nor explicit starting values are used,              *
 *  SINV chooses as starting values:                                         *
 *     newton's method:  x:     CDF(x) = 0.5                                 *
 *     regula falsi:     x1,x2: CDF(x1) = 1 - CDF(x2) = 0.05                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

#define SINV_MAX_U_LENGTH  (0.05)   /* maximal value for |u_i - u_{i-1}|     */
#define SINV_TAILCUTOFF    (1.e-10) /* Tail is cut off as max(TAILCUTOFF,e_u/0.1)*/
/** TODO: maschinenabhaengige konstante !!!! **/

#define SINV_XDEVIATION    (0.05)
/** TODO: beschreibung: how many percent may the approximation be away from the center of the 
interval that it is still accepted as splitting point (for the sake
of saves of CDF-evaluations.

konstante 0.005 optimiert fuer normalverteilung eu 1.e-10. Jedenfalls ist die
  Wahl von x1=arithmetic mean of the intervalborders "robuster".
Groesser XDEVIATION liefert eher mehr intervalle bei eu>1.e-10 aber weniger
CDF-Auswertungen im setup. 
bei eu=1.e-12 ist 0.3  am besten
insgesamt scheint am Anfang die Halbierung in x-Richtung, spaeter dann die 
Approximative in u-richtung besser zu sein*/

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define SINV_DEBUG_TABLE        0x00000010u   /* print table                 */
#define SINV_DEBUG_CHG          0x00001000u   /* print changed parameters    */
#define SINV_DEBUG_SAMPLE       0x01000000u   /* trace sampling              */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define SINV_SET_ORDER          0x001u  /* order of polynomial               */
#define SINV_SET_U_RESOLUTION   0x002u  /* maximal error in u                */
#define SINV_SET_GUIDEFACTOR    0x004u  /* relative size of guide table      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "SINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_sinv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_sinv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_sinv_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create the table with splines                                             */
/*---------------------------------------------------------------------------*/

static struct unur_sinv_interval *_unur_sinv_interval_new( struct unur_gen *gen, double p, double u );
/*---------------------------------------------------------------------------*/
/* make a new interval with node (u=F(p),p).                                 */
/*---------------------------------------------------------------------------*/

static struct unur_sinv_interval *_unur_sinv_interval_adapt( struct unur_gen *gen, 
							     struct unur_sinv_interval *iv );
/*---------------------------------------------------------------------------*/
/* check parameters in interval and split or truncate where necessary.       */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_interval_is_monotone( struct unur_gen *gen, struct unur_sinv_interval *iv );
/*---------------------------------------------------------------------------*/
/* check whether the given interval is monotone.                             */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_interval_parameter( struct unur_gen *gen, struct unur_sinv_interval *iv );
/*---------------------------------------------------------------------------*/
/* compute all parameter for interval (spline coefficients).                 */
/*---------------------------------------------------------------------------*/

static double _unur_sinv_eval_polynomial( double x, double *coeff, int order );
/*---------------------------------------------------------------------------*/
/* evaluate polynomial.                                                      */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_list_to_array( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy list of intervals into double array.                                 */
/*---------------------------------------------------------------------------*/

static int _unur_sinv_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_intervals( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into logfile.               */
/*---------------------------------------------------------------------------*/

static void _unur_sinv_debug_chg_truncated( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* trace changes of the truncated domain.                                    */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.sinv        /* data for parameter object         */
#define GEN       gen->data.sinv        /* data for generator object         */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define CDF(x)    _unur_cont_CDF((x),(gen->distr))  /* call to CDF           */
#define PDF(x)    _unur_cont_PDF((x),(gen->distr))  /* call to PDF           */
#define dPDF(x)   _unur_cont_dPDF((x),(gen->distr)) /* call to derivative of PDF */

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_sinv_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.cdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"CDF"); return NULL; }

  /* if default variant is Newton's method, then we also need the PDF ! */

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_SINV_PAR);

  /* copy input */
  par->distr   = distr;         /* pointer to distribution object            */

  /* set default values */
  PAR.order = (DISTR_IN.pdf) ? 3 : 1; /* order of polynomial                 */
  PAR.u_resolution = 1.0e-8;    /* maximal error allowed in u-direction      */
  PAR.guide_factor = 1.;        /* size of guide table / number of intervals */
  PAR.bleft = -1.e20;           /* left border of the computational domain   */
  PAR.bright = 1.e20;           /* right border of the computational domain  */

  par->method   = UNUR_METH_SINV; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_sinv_init;

  return par;

} /* end of unur_sinv_new() */

/*****************************************************************************/

int
unur_sinv_set_order( struct unur_par *par, int order)
     /*----------------------------------------------------------------------*/
     /* Set order of Hermite interpolation.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   order  ... order of interpolation polynome                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SINV );

  /* check new parameter for generator */
  if (order!=1 && order!=3 && order!=5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"order");
    return 0;
  }

  if (order > 1 && par->distr->data.cont.pdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    return 0;
  }

  if (order > 3 && par->distr->data.cont.dpdf == NULL) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dPDF");
    return 0;
  }

  /* store date */
  PAR.order = order;

  /* changelog */
  par->set |= SINV_SET_ORDER;

  return 1;

} /* end of unur_sinv_set_order() */

/*---------------------------------------------------------------------------*/

int
unur_sinv_set_u_resolution( struct unur_par *par, double u_resolution )
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   u_resolution ... maximal error in u                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SINV );

  /* check new parameter for generator */
  if (u_resolution < DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"u-resolution");
    return 0;
  }

  /* store date */
  PAR.u_resolution = u_resolution;

  /* changelog */
  par->set |= SINV_SET_U_RESOLUTION;

  return 1;

} /* end of unur_sinv_set_u_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_sinv_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,SINV );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return 0;
  }

  /* store date */
  PAR.guide_factor = factor;

  /* changelog */
  par->set |= SINV_SET_GUIDEFACTOR;

  return 1;

} /* end of unur_sinv_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int 
unur_sinv_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /* the new domain should not exceed the original domain given by        */
     /* unur_distr_cont_set_domain(). Otherwise it is truncated.             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  double Umin, Umax;

  /* check arguments */
  CHECK_NULL(gen, 0);
  _unur_check_gen_object(gen, SINV);

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* set bounds of U -- in respect to given bounds */
  Umin = (left > -INFINITY) ? CDF(left)  : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;

  /* check result */
  if (Umin > Umax) {
    /* this is a serios error that should not happen */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

  if (_unur_FP_equal(Umin,Umax)) {
    /* CDF values very close */
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (Umin == 0. || _unur_FP_same(Umax,1.)) {
      /* this is very bad */
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return 0;
    }
  }

  /* copy new boundaries into generator object */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN.Umin = Umin;
  GEN.Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & SINV_DEBUG_CHG) 
    _unur_sinv_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_sinv_chg_truncated() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

struct unur_gen *
_unur_sinv_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_SINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_SINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_sinv_create(par);
  if (!gen) { free(par); return NULL; }

  /* free parameters */
  free(par);
  
  /* domain not truncated at init */
  DISTR.trunc[0] = DISTR.domain[0];
  DISTR.trunc[1] = DISTR.domain[1];

  /* set bounds of U -- in respect to given bounds                          */
  GEN.Umin = GEN.CDFmin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
  GEN.Umax = GEN.CDFmax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;

  if (GEN.CDFmin > GEN.CDFmax) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF not increasing");
    _unur_sinv_free(gen); return NULL;
  }

  /* compute splines */
  if (!_unur_sinv_create_table(gen)) {
    _unur_sinv_free(gen);
    return NULL;
  }

  /* copy linked list into array */
  _unur_sinv_list_to_array( gen );

  /* make initial guide table */
  _unur_sinv_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_sinv_debug_init(gen);
#endif

  return gen;

} /* end of _unur_sinv_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_sinv_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_SINV_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_SINV_GEN);

  /* copy distribution object into generator object */
  gen->distr = _unur_distr_cont_clone( par->distr );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_sinv_sample;

  gen->destroy = _unur_sinv_free;
  gen->clone = _unur_sinv_clone;

  /* copy parameters into generator object */
  GEN.order = PAR.order;            /* order of polynomial                   */
  GEN.u_resolution = PAR.u_resolution; /* maximal error in u-direction       */
  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables      */
  GEN.bleft  = max(PAR.bleft,DISTR.domain[0]);
  GEN.bright = min(PAR.bright,DISTR.domain[1]);

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */

  /* initialize variables */
  GEN.N = 0;
  GEN.guide_size = 0; 
  GEN.iv = NULL;
  GEN.intervals = NULL;
  GEN.guide = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_sinv_create() */

/*---------------------------------------------------------------------------*/

int
_unur_sinv_create_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a table of splines                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer generator object                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  struct unur_sinv_interval *iv;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_SINV_GEN,0);
  
  /* first and terminating interval */
  GEN.iv = _unur_sinv_interval_new(gen,GEN.bleft,CDF(GEN.bleft));
  if (GEN.iv == NULL) return 0;
  GEN.iv->next = _unur_sinv_interval_new(gen,GEN.bright,CDF(GEN.bright));
  if (GEN.iv->next == NULL) return 0;

  for (iv=GEN.iv; iv->next!=NULL; ) {
    COOKIE_CHECK(iv,CK_SINV_IV,0);
    iv = _unur_sinv_interval_adapt(gen,iv);
  }

  /* last interval is only used to store right boundary */
  iv->spline[0] = iv->p;

  /* o.k. */
  return 1;
}  /* end of _unur_sinv_create_table() */

/*---------------------------------------------------------------------------*/

static struct unur_sinv_interval *
_unur_sinv_interval_adapt( struct unur_gen *gen, struct unur_sinv_interval *iv )
     /*----------------------------------------------------------------------*/
     /* check parameters in interval and split or truncate where necessary.  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   iv       ... if splitted                                           */
     /*   iv->next ... if interval was o.k.                                  */
     /*----------------------------------------------------------------------*/
{
  double p_new;   /* new design point */
  struct unur_sinv_interval *iv_new, *iv_tmp;
  double x, Fx;

  /* 1st check: right most interval (of at least 2)
     with CDF greater than 1.- TAILCUTOFF */

  iv_tmp = iv->next->next;
  if(iv_tmp && iv->next->u > 1. - min(SINV_TAILCUTOFF, 0.1*GEN.u_resolution)) {
    /* chop off right hand tail */
    free (iv_tmp);
    iv->next->next = NULL;
    GEN.N--;
    /* update right boundary */
    GEN.bright = iv->next->p;
    return iv;
  }

  /* 2nd check: is the left most interval (of at least 2) 
     with CDF less than TAILCUTOFF */

  if (iv==GEN.iv && iv->next->next && iv->next->u < min(SINV_TAILCUTOFF, 0.1*GEN.u_resolution)) {
    /* chop off left hand tail */
    iv_tmp = GEN.iv;
    GEN.iv = iv->next;
    free (iv_tmp);
    GEN.N--;
    /* update left boundary */
    GEN.bleft = GEN.iv->p;
    return GEN.iv;
  }

  /* 3rd check: |u_i - u_{i-1}| must not exceed threshold value */
  /* 4th check: monotonicity                                    */

  if ( (iv->next->u - iv->u > SINV_MAX_U_LENGTH) ||
       (! _unur_sinv_interval_is_monotone(gen,iv)) ) {
    /* mean of x-interval as splitting point */
    p_new = 0.5 * (iv->next->p + iv->p);
    /* insert new interval into linked list */
    iv_new = _unur_sinv_interval_new(gen,p_new,CDF(p_new));
    iv_new->next = iv->next;
    iv->next = iv_new;
    return iv;
  }

  /* compute coefficients for spline (only necessary if monotone) */
  _unur_sinv_interval_parameter(gen,iv);

  /* 5th check: error in u-direction */

  /* compute approximate value for inverse CDF in center of interval */
  x = _unur_sinv_eval_polynomial( 0.5, iv->spline, GEN.order );
  Fx = CDF(x);

  if (fabs(Fx - 0.5*(iv->next->u + iv->u)) > GEN.u_resolution) {
    /* error in u-direction too large */
    p_new = 0.5 * (iv->next->p + iv->p);
    /* if possible we use the point x instead of p_new */
    if(fabs(p_new-x)< SINV_XDEVIATION * (iv->next->p - iv->p))
      iv_new = _unur_sinv_interval_new(gen,x,Fx);
    else
      iv_new = _unur_sinv_interval_new(gen,p_new,CDF(p_new));
    iv_new->next = iv->next;
    iv->next = iv_new;
    return iv;
  }

  /* interval o.k. */
  return iv->next;

} /* end of _unur_sinv_interval_adapt() */

/*---------------------------------------------------------------------------*/

static struct unur_sinv_interval *
_unur_sinv_interval_new( struct unur_gen *gen, double p, double u )
     /*----------------------------------------------------------------------*/
     /* make a new interval with node (u=F(p),p).                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   p   ... left design point of new interval                          */
     /*   u   ... value of CDF at p, u=CDF(p)                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_sinv_interval *iv;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SINV_GEN,NULL);

  /* first check u */
  if (u<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"CDF(x) < 0.");
    return NULL;
  }

  /* we need new interval */
  iv = _unur_malloc( sizeof(struct unur_sinv_interval) );
  COOKIE_SET(iv,CK_SINV_IV);

  /* compute and store data */
  switch (GEN.order) {
  case 5:
    iv->df = dPDF(p);
  case 3:
    iv->f = PDF(p);
  case 1:
    iv->p = p;
    iv->u = u;
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    free(iv);
    return NULL;
  }

  iv->next = NULL;  /* add eol marker */
  ++(GEN.N);   /* increment counter for intervals */

  /* o.k. */
  return iv;

} /* end of _unur_sinv_interval_new() */

/*---------------------------------------------------------------------------*/

int 
_unur_sinv_interval_is_monotone( struct unur_gen *gen, struct unur_sinv_interval *iv )
     /*----------------------------------------------------------------------*/
     /* check whether the given interval is monotone.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if monotone                                                  */
     /*   0 ... otherwise                                                    */
     /*----------------------------------------------------------------------*/
{

  double bound;

  switch (GEN.order) {
  case 5:
    /* the monotone check is in the moment only implemented for order 3
       as approximation we use the same check for order 5*/
  case 3:
    /* difference quotient */
    bound = 3.*(iv->next->p - iv->p)/(iv->next->u - iv->u);
    return (1./iv->next->f > bound || 1./iv->f > bound) ? 0 : 1;
  case 1:
    /* linear interpolation is always monotone */
  default:  /* we assume that we have checked GEN.order very often till now */
    return 1;
  }

} /* end of _unur_sinv_interval_is_monotone() */

/*---------------------------------------------------------------------------*/

int
_unur_sinv_interval_parameter( struct unur_gen *gen, struct unur_sinv_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute all parameter for interval (spline coefficients).            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  double delta_u, delta_p;
  double f1, fs0, fs1, fss0, fss1;

  delta_u = iv->next->u - iv->u;
  delta_p = iv->next->p - iv->p;

  switch (GEN.order) {

  case 1:    /* linear interpolation */
    iv->spline[0] = iv->p;
    iv->spline[1] = delta_p;
    return 1;

  case 3:    /* cubic Hermite interpolation */
    iv->spline[0] = iv->p;
    iv->spline[1] = delta_u / iv->f;
    iv->spline[2] = 3.* delta_p - delta_u * (2./iv->f + 1./iv->next->f);
    iv->spline[3] = -2.* delta_p + delta_u * (1./iv->f + 1./iv->next->f);
    return 1;

  case 5:    /* quintic Hermite interpolation */
    f1   = delta_p;
    fs0  = delta_u / iv->f;      
    fs1  = delta_u / iv->next->f;
    fss0 = -delta_u * delta_u * iv->df / (iv->f * iv->f * iv->f);
    fss1 = -delta_u * delta_u * iv->next->df / (iv->next->f * iv->next->f * iv->next->f);

    iv->spline[0] = iv->p;
    iv->spline[1] = fs0;
    iv->spline[2] = 0.5*fss0;
    iv->spline[3] = 10.*f1 - 6.*fs0 - 4.*fs1 - 1.5*fss0 + 0.5*fss1;
    iv->spline[4] = -15.*f1 + 8.*fs0 + 7.*fs1 + 1.5*fss0 - fss1;
    iv->spline[5] = 6.*f1 - 3.*fs0 - 3.*fs1 - 0.5*fss0 + 0.5*fss1;
    return 1;

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

} /* end of _unur_sinv_interval_parameter() */

/*---------------------------------------------------------------------------*/

double
_unur_sinv_eval_polynomial( double x, double *coeff, int order )
     /*----------------------------------------------------------------------*/
     /* evaluate polynomial using Horner scheme.                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument                                                 */
     /*   coeff ... coefficients of polynomial (increasing order)            */
     /*   order ... order of polynomial                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   value of spline at u                                               */
     /*----------------------------------------------------------------------*/
{
  int i;
  double poly;

  poly = coeff[order];
  for (i=order-1; i>=0; i--)
    poly = x*poly + coeff[i];

  return poly;
} /* end of _unur_sinv_eval_polynomial() */

/*---------------------------------------------------------------------------*/

int
_unur_sinv_list_to_array( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy list of intervals into double array.                            */
     /* the linked list is freed.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  int i; 
  struct unur_sinv_interval *iv, *next;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_SINV_GEN,0);

  /* allocate memory */
  GEN.intervals = malloc( GEN.N*(GEN.order+2) * sizeof(double) );

  i = 0;
  for (iv=GEN.iv; iv!=NULL; iv=next) {
    /* copy */
    GEN.intervals[i] = iv->u;
    memcpy( GEN.intervals+(i+1), &(iv->spline), (GEN.order+1)*sizeof(double) );
    i += GEN.order+2;
    /* and free linked list */
    next = iv->next;
    free(iv);
  }

  /* linked list is now empty */
  GEN.iv = NULL;

  /* o.k. */
  return 1;
} /* end of _unur_sinv_list_to_array() */

/*---------------------------------------------------------------------------*/
static int
_unur_sinv_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1  ... if successful                                               */
     /*   0  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  int i,j, imax;
  double delta_U;

  /* check arguments */
  CHECK_NULL(gen,0);  COOKIE_CHECK(gen,CK_SINV_GEN,0);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN.guide) {
    GEN.guide_size = GEN.N * GEN.guide_factor;
    if (GEN.guide_size <= 0) GEN.guide_size = 1; 
    GEN.guide = _unur_malloc( GEN.guide_size * sizeof(int) );
  }

  delta_U = GEN.CDFmax - GEN.CDFmin;
  imax = (GEN.N-1) * (GEN.order+2);

  /* u value at end of interval */
#define u(i)  (((GEN.intervals[(i)+GEN.order+2])-GEN.CDFmin)/delta_U)

  i = 0;
  GEN.guide[0] = 0;
  for( j=1; j<GEN.guide_size ;j++ ) {
    while( u(i) < (j/(double)GEN.guide_size) )
      i += GEN.order+2;
    if (i >= imax) {
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
      break;
    }
    GEN.guide[j]=i;
  }

#undef u

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide[j] = i;

  /* o.k. */
  return 1;
} /* end of _unur_sinv_make_guide_table() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_sinv_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE clone->data.sinv

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_SINV_GEN,NULL);

  /* allocate memory for generator object */
  clone = _unur_malloc( sizeof(struct unur_gen) );

  /* copy main part */
  memcpy( clone, gen, sizeof(struct unur_gen) );

  /* set generator identifier */
  clone->genid = _unur_set_genid(GENTYPE);

  /* copy distribution object into generator object */
  clone->distr = _unur_distr_cont_clone( gen->distr );

  /* auxiliary generator */
  if (gen->gen_aux) clone->gen_aux = unur_gen_clone( gen->gen_aux );

  /* copy tables for generator object */
  CLONE.intervals = _unur_malloc( GEN.N*(GEN.order+2) * sizeof(double) );
  memcpy( CLONE.intervals, GEN.intervals, GEN.N*(GEN.order+2) * sizeof(double) );
  CLONE.guide = _unur_malloc( GEN.guide_size * sizeof(int) );
  memcpy( CLONE.guide, GEN.guide, GEN.guide_size * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_sinv_clone() */

/*****************************************************************************/

void
_unur_sinv_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_SINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free linked list of intervals */
  if (GEN.iv) {
    struct unur_sinv_interval *iv,*next;
    for (iv = GEN.iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free tables */
  if (GEN.intervals) free (GEN.intervals);
  if (GEN.guide)     free (GEN.guide);

  /* free memory */
  _unur_distr_free(gen->distr);
  _unur_free_genid(gen);

  COOKIE_CLEAR(gen);
  free(gen);

} /* end of _unur_sinv_free() */

/*****************************************************************************/

double
_unur_sinv_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  double U;
  int i;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_SINV_GEN,0.);

  /* sample from U( Umin, Umax ) */
  U = GEN.Umin + _unur_call_urng(gen->urng) * (GEN.Umax - GEN.Umin);

  /* look up in guide table and search for interval */
  i =  GEN.guide[(int) (GEN.guide_size*(U-GEN.CDFmin)/(GEN.CDFmax-GEN.CDFmin))];
  while (U > GEN.intervals[i+GEN.order+2])
    i += GEN.order+2;

  /* rescale uniform random number */
  U = (U-GEN.intervals[i])/(GEN.intervals[i+GEN.order+2] - GEN.intervals[i]);

  /* evaluate polynome */
  return _unur_sinv_eval_polynomial( U, GEN.intervals+i+1, GEN.order );

} /* end of _unur_sinv_sample() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_sinv_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = SINV (Spline approximation of INVerse CDF)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(log,"%s: sampling routine = _unur_sinv_sample\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: order of polynomial = %d",gen->genid,GEN.order);
  _unur_print_if_default(gen,SINV_SET_ORDER);
  fprintf(log,"\n%s: u-resolution = %g",gen->genid,GEN.u_resolution);
  _unur_print_if_default(gen,SINV_SET_U_RESOLUTION);
  fprintf(log,"\n%s: domain of computation = [%g,%g]\n",gen->genid,GEN.bleft,GEN.bright);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.Umin,GEN.Umax);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN.guide_factor);
  _unur_print_if_default(gen,SINV_SET_GUIDEFACTOR);
  fprintf(log,"\n%s:\n",gen->genid);

  _unur_sinv_debug_intervals(gen);

  fflush(stdout);

} /* end of _unur_sinv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_sinv_debug_intervals( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print table of intervals into logfile                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int i,n;
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: Intervals: %d\n",gen->genid,GEN.N-1);

  if (gen->debug & SINV_DEBUG_TABLE) {
    fprintf(log,"%s:   Nr.      u=CDF(p)     p=spline[0]   spline[1]    ...\n",gen->genid);
    for (n=0; n<GEN.N-1; n++) {
      i = n*(GEN.order+2);
      fprintf(log,"%s:[%4d]: %#12.6g  %#12.6g  %#12.6g", gen->genid, n,
	      GEN.intervals[i], GEN.intervals[i+1], GEN.intervals[i+2]);
      if (GEN.order>1)
	fprintf(log,"  %#12.6g  %#12.6g", GEN.intervals[i+3], GEN.intervals[i+4]);
      if (GEN.order>3)
	fprintf(log,"  %#12.6g  %#12.6g", GEN.intervals[i+5], GEN.intervals[i+6]);
      fprintf(log,"\n");
    }
    i = n*(GEN.order+2);
    fprintf(log,"%s:[%4d]: %#12.6g  %#12.6g  (right boundary)\n", gen->genid, n,
	    GEN.intervals[i], GEN.intervals[i+1]);
  }

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_sinv_debug_intervals() */

/*---------------------------------------------------------------------------*/

void 
_unur_sinv_debug_chg_truncated( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_SINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.Umin,GEN.Umax);

} /* end of _unur_sinv_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
