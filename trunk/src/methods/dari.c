/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dari.c                                                       *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    discrete automatic rejection inversion                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PMF of a T(-1/2)-concave distribution;                         *
 *      produce a value x consistent with its PMF.                           *
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
 *   [1] Hoermann W. (199?):                                                 *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DARI_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DARI_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DARI_SET_CFACTOR        0x001u
#define DARI_SET_SIZE           0x002u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DARI"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dari_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute hat and squeezes.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dari_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dari_sample( struct unur_gen *generator );
static int _unur_dari_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_dari_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dari_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       par->data.dari        /* data for parameter object         */
#define GEN       gen->data.dari        /* data for generator object         */
#define DISTR     gen->distr.data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */     

#define PMF(x)    _unur_cont_PMF((x),&(gen->distr))   /* call to PMF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_dari_new( struct unur_distr *distr )
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
     /*                                                                      */
     /* comment:                                                             */
     /*   if the sum over the PMF is not close to 1 it is necessary to       */
     /*   set pmf_sum to an approximate value of its sum (+/- 30 % is ok).  */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }

  /** TODO: brauchst du?? **/
  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    /** TODO: das mit dem mode funktioniert noch nicht so wie geplant.
	die schnittstelle fuer DISCR ist fertig,
	fuer die standard verteilungen in UNURAN habe ich es noch nicht
	gemacht **/
    if (!unur_distr_discr_upd_mode(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return NULL; 
    }
  }

  /** TODO: brauchst du?? **/
  if (!(distr->set & UNUR_DISTR_SET_PMFSUM))
    if (!unur_distr_discr_upd_pmfsum(distr)) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return NULL; 
    }

  /** brauchst du sonst noch was ? **/

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_DARI_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */

  /** TODO: default werte **/

  PAR.c_factor  = 0.664;   /** TODO: ?? **/
          /* optimal value for the normal distribution, which is good for 
	     all bell-shaped densities. The minimax approach for that 
	     transformation has c_factor=2. */

  PAR.squeeze   = 0;          /* no (??) squeezes by default */
  PAR.size      = 10;          /* size of table for speeding up generation      */

  par->method   = UNUR_METH_DARI;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dari_init;

  return par;

} /* end of unur_dari_new() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_dari_set_cfactor( struct unur_par *par, double cfactor )
     /*----------------------------------------------------------------------*/
     /* set factor for position of left and right construction point         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   cfactor ... factor                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* check new parameter for generator */
  /** TODO: welche werte fuer c sind zulaessig / sinnvoll ? 
  zulaessig ist jedes c>0, man koennte damit falsche Flaechenangaben kompensieren.
  Wenn sum genau bekannt ist, ist ein c > 2 (2 ist der minimax approach) so weit
  ich weiss nie sinnvoll. Ich denke aber, das sollte man besser nicht prinzipiell
  verbieten, hoechstens eine warnung.**/
  if (cfactor <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c-factor <= 0");
    return 0;
  }

  if (cfactor > 2.1)
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c-factor > 2 not recommended. skip");

  /* store date */
  PAR.c_factor = cfactor;

  /* changelog */
  par->set |= DARI_SET_CFACTOR;

  return 1;

} /* end of unur_dari_set_cfactor() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_squeeze( struct unur_par *par, char squeeze )
     /*----------------------------------------------------------------------*/
     /* turn on/off using squeezes                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   squeeze ... 0 = no squeeze,  !0 = use squeeze                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );
  
  /* check input */
  _unur_check_par_object( par,DARI );

  /* store data */
  PAR.squeeze = squeeze;

  /* o.k. */
  return 1;

} /* end of unur_dari_set_squeeze() */

/*---------------------------------------------------------------------------*/

int
unur_dari_set_size( struct unur_par *par, int size )
     /*----------------------------------------------------------------------*/
     /* set size of table                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   size   ... table size                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* check parameter */
  /** TODO **/
  if (size <= 0 || size > 1000) {  
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"invalid table size");
    return 0;
  }
  
  /* store data */
  PAR.size = size;

  /* changelog */
  par->set |= DARI_SET_SIZE;

  /* o.k. */
  return 1;
} /* end of unur_dari_set_size() */
  
/*---------------------------------------------------------------------------*/

int
unur_dari_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,DARI );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | DARI_VARFLAG_VERIFY) : (par->variant & (~DARI_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_dari_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL( gen,0 );

  /* check input */
  _unur_check_gen_object( gen,DARI );

  if (verify) {
    /* turn verify mode on */
    gen->variant |= DARI_VARFLAG_VERIFY;
    SAMPLE = _unur_dari_sample_check;
  }
  else {
    /* turn verify mode off */
    gen->variant &= ~DARI_VARFLAG_VERIFY;
    SAMPLE = _unur_dari_sample;
  }

  /* o.k. */
  return 1;

} /* end of unur_dari_chg_verify() */

/*****************************************************************************/

int
unur_dari_chg_pmfparams( struct unur_gen *gen, double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* change array of parameters for distribution                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* IMPORTANT: The given parameters are not checked against domain       */
     /*            errors (in opposition to the unur_<distr>_new() call).    */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  register int i;

  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );
  if (n_params>0) _unur_check_NULL(gen->genid,params,0);
  
  /* check new parameter for generator */
  if (n_params > UNUR_DISTR_MAXPARAMS || n_params < 0 ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* copy parameters */
  DISTR.n_params = n_params;
  for (i=0; i < n_params; i++)
    DISTR.params[i] = params[i];

  /* changelog */
  /* mode and sum might be wrong now! 
     but the user is responsible to change it.
     so we dont say:
     gen->distr.set &= ~(UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_PMFSUM );
     gen->set &= ~DARI_SET_CDFMODE;
  */

  /* o.k. */
  return 1;
} /* end of unur_dari_chg_pmfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_dari_chg_domain( struct unur_gen *gen, int left, int right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity.           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_dari_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_mode( struct unur_gen *gen, int mode )
     /*----------------------------------------------------------------------*/
     /* change mode of distribution                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   mode  ... mode                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );
  
  /* copy parameters */
  DISTR.mode = mode;

  /** no changelog required ? **/

  /* o.k. */
  return 1;
} /* end of unur_dari_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_mode( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute mode of distribution                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );

  return unur_distr_discr_upd_mode( &(gen->distr) );
} /* end of unur_dari_upd_mode() */

/*---------------------------------------------------------------------------*/

int
unur_dari_chg_pmfsum( struct unur_gen *gen, double sum )
     /*----------------------------------------------------------------------*/
     /* change sum over PMF of distribution                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   sum   ... sum                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );
  
  /* check new parameter for generator */
  if (sum <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"sum <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.sum = sum;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_dari_chg_pmfsum() */

/*---------------------------------------------------------------------------*/

int
unur_dari_upd_pmfsum( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* recompute sum over PMF of distribution                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );

  return unur_distr_discr_upd_pmfsum( &(gen->distr) );
} /* end of unur_dari_upd_pmfsum() */

/*****************************************************************************/

struct unur_gen *
_unur_dari_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
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
  if ( par->method != UNUR_METH_DARI ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DARI_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dari_create(par);
  if (!gen) return NULL;

  /* create hat and squeeze (setup procedure) */
  if ( !_unur_dari_hat(gen) ) {
    /* error */
    _unur_dari_free(gen);
    return NULL;
  }

  /* hat successfully created */
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_dari_debug_init(gen);
#endif

  /* free parameters */
  free(par);

  /* o.k. */
  return gen;

} /* end of _unur_dari_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dari_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute hat and squeeze                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{ 
  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_DARI_GEN, 0 );

  /* o.k. */
  return 1;

} /* end of _unur_dari_hat() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_dari_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DARI_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DARI_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->variant & DARI_VARFLAG_VERIFY) ? _unur_dari_sample_check : _unur_dari_sample;
  gen->destroy = _unur_dari_free;

  /* copy some parameters into generator object */

  /** TODO ? **/


  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* initialize parameters */
  /** TODO ? **/

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_dari_create() */

/*---------------------------------------------------------------------------*/

int
unur_dari_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,DARI );

  /* compute universal bounding rectangle */
  return _unur_dari_hat( gen );
} /* end of unur_dari_reinit() */

/*****************************************************************************/

int
_unur_dari_sample( struct unur_gen *gen )
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
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DARI_GEN,0.);

  
  return 1.;
} /* end of _unur_dari_sample() */

/*****************************************************************************/

int
_unur_dari_sample_check( struct unur_gen *gen )
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
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_DARI_GEN,0.);
  
  return 1;
} /* end of _unur_dari_sample_check() */

/*****************************************************************************/

void
_unur_dari_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DARI ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DARI_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of _unur_dari_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dari_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_DARI_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = dari (discrete automatic rejection inversion)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_discr_debug( &(gen->distr), gen->genid, FALSE );

  fprintf(log,"%s: sampling routine = _unur_dari_sample",gen->genid);
  if (gen->variant & DARI_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: Data for hat and squeeze:\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_dari_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
