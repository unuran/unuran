/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    numerical inversion of cumulative distribution function      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the c.d.f.                                                *
 *      newton's method: additional pointer to the p.d.f.                    *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      c.d.f. at mode                                                       *
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
 *   [1] Neumaier A. (to be publishes): Introduction to numerical analysis,  *
 *       Cambridge University Press                                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Numerical inversion is a method for generating random variables          *
 *  using the cdf (and in case of newton's method the pdf).                  *
 *                                                                           *
 *  THEOREM:                                                                 *
 *     Let X be a random variable with cdf F(x).                             *
 *     Then the F(X) are  uniformly distributed.                             *
 *                                                                           *
 *  COROLLARY:                                                               *
 *     Starting with uniformly distributed random variables U,               *
 *     the F^(-1)(U) have F(x) as cdf.                                       *
 *                                                                           *
 *  Starting with an U, the task is to find a X fulfilling:                  *
 *    F(X) - U = 0.                                                          *
 *                                                                           *
 *  Numerical algorithms to find zeros that are used in NINV are variants of * 
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
 *  NINV chooses as starting values:                                         *
 *     newton's method:  x:     cdf(x) = 0.5                                 *
 *     regula falsi:     x1,x2: cdf(x1) = 1 - cdf(x2) = 0.05                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <source_unuran.h>

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* Starting interval includes this percentage of all univariate rand numbers */
/* must be > 0. and < 1.                                                     */
#define INTERVAL_COVERS  (.9)

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define NINV_VARFLAG_NEWTON   0x1u   /* use Newton's method                  */
#define NINV_VARFLAG_REGULA   0x2u   /* use regula falsi (default)           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NINV_DEBUG_SAMPLE       0x01000000u
#define NINV_DEBUG_CHG          0x00001000u   /* print changed parameters    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NINV_SET_MAX_ITER     0x001u   /* number of maximal interations      */
#define NINV_SET_X_RESOLUTION 0x002u   /* maximal relative error in x        */
#define NINV_SET_START        0x004u   /* intervals at start (left/right)    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_sample_regula( struct unur_gen *gen );
static double _unur_ninv_sample_newton( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_regula( struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* algorithm: regula falsi                                                   */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_newton( struct unur_gen *gen, double u);
/*---------------------------------------------------------------------------*/
/* algorithm: newton method                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_create_table(UNUR_GEN *gen);
/*---------------------------------------------------------------------------*/
/* create the table with starting points                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample_regula( struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling (regula falsi).                                            */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample_newton( struct unur_gen *gen, 
					    double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling (newton's method).                                         */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_chg_domain(UNUR_GEN *gen);
/*---------------------------------------------------------------------------*/
/* trace changes of the domain.                                              */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.ninv        /* data for parameter object         */
#define GEN       gen->data.ninv        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]      /* left boundary of domain of distribution  */
#define BD_RIGHT  domain[1]      /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)    _unur_cont_PDF((x),&(gen->distr))   /* call to p.d.f.      */
#define CDF(x)    _unur_cont_CDF((x),&(gen->distr))   /* call to c.d.f.      */

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_ninv_new( struct unur_distr *distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"c.d.f."); return NULL; }

  /* if default variant is Newton's method, then we also need the p.d.f. ! */

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_NINV_PAR);

  /* copy input */
  par->distr   = distr;           /* pointer to distribution object          */

  /* set default values */
  PAR.max_iter  = 40;             /* maximal number of iterations            */
  PAR.rel_x_resolution = 1.0e-8;  /* maximal relative error allowed in x     */

  /* starting points for numerical inversion */
  PAR.s[0]      = 0.0;     /* regula falsi: left boundary of starting interval
			      newton: starting point                         */
  PAR.s[1]      = 0.0;     /* regula falsi: right boundary of starting interval
			      newton: not used                               */
  /* If s1 and s2 are equal a defaults are used, see below */

  PAR.table_on  = FALSE;   /* Do not use a table for starting points
			      by default.                                    */
 
  par->method   = UNUR_METH_NINV;          /* method and default variant     */
  par->variant  = NINV_VARFLAG_REGULA;     /* Use regula falsi as default 
					      method                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_ninv_init;

  return par;

} /* end of unur_ninv_new() */

/*****************************************************************************/

int unur_ninv_set_usenewton( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use Newton's method                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );

  /* check new parameter for generator */
  if (! par->DISTR_IN.pdf) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.d.f.");
    par->variant = NINV_VARFLAG_REGULA;   /* use regula falsi instead  */
    return 0;
 }

  /* store date */
  par->variant = NINV_VARFLAG_NEWTON;

  return 1;

} /* end of unur_ninv_set_usenewton() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_useregula( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* use regula falsi                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );

  /* store date */
  par->variant = NINV_VARFLAG_REGULA;

  return 1;

} /* end of unur_ninv_set_useregula() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_max_iter(UNUR_PAR *par, int max_iter )
     /*----------------------------------------------------------------------*/
     /* set number of maximal iterations                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximal iterations");
    return 0;
  }

  /* store date */
  PAR.max_iter = max_iter;

  /* changelog */
  par->set |= NINV_SET_MAX_ITER;

  return 1;

} /* end of unur_ninv_set_max_iter() */

/*---------------------------------------------------------------------------*/

int unur_ninv_chg_max_iter(UNUR_GEN *gen, int max_iter )
     /*----------------------------------------------------------------------*/
     /* change number of maximal iterations                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   max_iter ...  number of maximal iterations                         */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, 0);

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(gen->genid, UNUR_ERR_PAR_SET, "maximal iterations");
    return 0;
  }

  /* store date */
  GEN.max_iter = max_iter;

  /* changelog */
  gen->set |= NINV_SET_MAX_ITER;

  /* o.k.  */
  return 1;

} /* end of unur_ninv_chg_max_iter() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_x_resolution(UNUR_PAR *par, double x_resolution)
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter for building generator object*/
     /*   x_resolution ... maximal relative error in x                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );

  /* check new parameter for generator */
  if (x_resolution < DBL_EPSILON) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"x resolution");
    return 0;
  }

  /* store date */
  PAR.rel_x_resolution = x_resolution;

  /* changelog */
  par->set |= NINV_SET_X_RESOLUTION;

  return 1;

} /* end of unur_ninv_set_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int unur_ninv_chg_x_resolution(UNUR_GEN *gen, double x_resolution)
     /*----------------------------------------------------------------------*/
     /* set maximal relative error in x                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   x_resolution ... maximal relative error in x                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, 0);

  /* check new parameter for generator */
  if (x_resolution < DBL_EPSILON) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"x resolution");
    return 0;
  }

  /* store date */
  GEN.rel_x_resolution = x_resolution;

  /* changelog */
  gen->set |= NINV_SET_X_RESOLUTION;

  return 1;

} /* end of unur_ninv_chg_x_resolutuion() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_set_start(UNUR_PAR *par, double s1, double s2)
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /* arguments that are used by method are ignored.                       */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*                                                                      */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );

  /* store date */
  if ( s1 <= s2 ){
     PAR.s[0] = s1;
     PAR.s[1] = s2;
  }
  else{
     PAR.s[0] = s2;
     PAR.s[1] = s1;
  }

  /* changelog */
  par->set |= NINV_SET_START;

  return 1;

} /* end of unur_ninv_set_start() */

/*---------------------------------------------------------------------------*/

int unur_ninv_chg_start(UNUR_GEN *gen, double s1, double s2)
     /*----------------------------------------------------------------------*/
     /* set starting points.                                                 */
     /*   Newton:        s1           starting point                         */
     /*   regular falsi: s1, s2       boundary of starting interval          */
     /* arguments that are used by method are ignored.                       */
     /*                                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /*                                                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen, 0);

  /* store date */
  if ( s1 <= s2 ){
     GEN.s[0] = s1;
     GEN.s[1] = s2;
  }
  else{
     GEN.s[0] = s2;
     GEN.s[1] = s1;
  }

 if ( _FP_same(GEN.s[0], GEN.s[1]) && !GEN.table_on) {
      /* length of interval == 0 -> choose bounderies with                   */
      /*  INTERVAL_COVERS *100% chance for sign change in interval           */
      GEN.s[0] = -10.;      /* arbitrary starting value                      */
      GEN.s[1] =  10.;      /* arbitrary starting value                      */
      GEN.s[0] = _unur_ninv_regula(gen, (1.-INTERVAL_COVERS)/2. );
      GEN.s[1] = GEN.s[0] + 10.;   /* arbitrary interval length              */
      GEN.s[1] = _unur_ninv_regula(gen, (1.+INTERVAL_COVERS)/2. );
    }

  /* changelog */
  gen->set |= NINV_SET_START;

  return 1;

} /* end of unur_ninv_chg_start() */

/*---------------------------------------------------------------------------*/

int unur_ninv_set_table(UNUR_PAR *par, int tbl_pnts)
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* the function unur_ninv_set_start() is overruled                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,NINV );
  PAR.table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  PAR.table_on = TRUE;

  return 1;

} /* end of unur_ninv_set_table() */

/*---------------------------------------------------------------------------*/

int unur_ninv_chg_table(UNUR_GEN *gen, int tbl_pnts)
     /*----------------------------------------------------------------------*/
     /* if used, a table is generated to find better startig points          */
     /* set somewhere else will be ignored                                   */
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
  CHECK_NULL(gen, 0);

  free(GEN.table); /* new table will be generated */  
  GEN.table_size = (tbl_pnts >= 10) ? tbl_pnts : 10;
  
  _unur_ninv_create_table(gen);  

  /* ok */
  return 1;

} /* end of unur_ninv_chg_table() */

/*---------------------------------------------------------------------------*/

int
_unur_ninv_create_table(UNUR_GEN *gen)
{
  int i;

  /* check arguments */
  CHECK_NULL(gen, 0);
  _unur_check_gen_object(gen, NINV);

  GEN.table = (double *) malloc(GEN.table_size*sizeof(double));

  GEN.s[0] = - 10.;  /* arbitrary starting values                        */
  GEN.s[1] =   10.;
  GEN.table_on = FALSE;   /* table can't be used to calculate itself     */

  GEN.CDFmin = GEN.Umin;
  GEN.CDFmax = GEN.Umax;

  GEN.table[0] = DISTR.BD_LEFT;
  GEN.table[GEN.table_size-1] = DISTR.BD_RIGHT;
       
  for (i=1; i<GEN.table_size/2; i++){
    GEN.table[i] =
      _unur_ninv_regula(gen, i/(GEN.table_size-1.) );
    GEN.table[GEN.table_size-1-i] = 
      _unur_ninv_regula(gen, ((GEN.table_size-i-1))/(GEN.table_size-1.) );

    GEN.s[0] = (GEN.table[i] <= -INFINITY)? GEN.table[i+1]: GEN.table[i]; 
    GEN.s[1] = (GEN.table[i] >= INFINITY)?
      GEN.table[GEN.table_size-i-2]: GEN.table[GEN.table_size-i-1];
 
  }  /* end of for()                                                     */

  if (GEN.table_size/2 != GEN.table_size/2.){ /* GEN.table_size is odd ? */
    GEN.table[GEN.table_size/2] =
      _unur_ninv_regula(gen, ((GEN.table_size/2))/(GEN.table_size-1.) );
  }

  GEN.table_on = TRUE;

  /* o.k. */
  return 1;

}  /* end of _unur_ninv_create_table() */

/*---------------------------------------------------------------------------*/

int 
unur_ninv_chg_domain(UNUR_GEN *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the distribution  */
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

  /* check arguments */
  CHECK_NULL(gen, 0);
  _unur_check_gen_object(gen, NINV);

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }


  /* copy new boundaries into generator object */
  DISTR.BD_LEFT  = left;
  DISTR.BD_RIGHT = right;
  if (GEN.table_on && left < GEN.table[0]){
    _unur_warning(NULL, UNUR_ERR_DISTR_SET,
        "left boarder of domain exceeds range of table -> truncated");
    DISTR.BD_LEFT = GEN.table[0];
  }
  if (GEN.table_on && right > GEN.table[GEN.table_size-1]){
    _unur_warning(NULL, UNUR_ERR_DISTR_SET,
        "right boarder of domain exceeds range of table -> truncated");
    DISTR.BD_LEFT = GEN.table[GEN.table_size-1];
  }

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_DOMAIN;

  /* set bounds of U -- in respect to given bounds */
  GEN.Umin = (DISTR.BD_LEFT  <= -INFINITY) ? 0.0 : CDF(DISTR.BD_LEFT); 
  GEN.Umax = (DISTR.BD_RIGHT >=  INFINITY) ? 1.0 : CDF(DISTR.BD_RIGHT); 

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug & NINV_DEBUG_CHG) 
    _unur_ninv_debug_chg_domain( gen );
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_ninv_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_ninv_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object(gen, NINV);
  if (n_params>0) CHECK_NULL(params, 0);
  
  /* check new parameter for generator */
  if (n_params > UNUR_DISTR_MAXPARAMS || n_params < 0 ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* copy parameters */
  DISTR.n_params = n_params;
  memcpy(DISTR.params, params, n_params * sizeof(double));

  /* changelog */
  /* mode and area might be wrong now! 
     but the user is responsible to change it.
     so we dont say:
     gen->distr.set &= ~(UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_PDFAREA );
     gen->set &= ~NINV_SET_CDFMODE;
  */

  /* o.k. */
  return 1;
} /* end of unur_ninv_chg_pdfparams() */

/*****************************************************************************/

struct unur_gen *
_unur_ninv_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_NINV ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_ninv_create(par);
  if (!gen) { free(par); return NULL; }

  /* set bounds of U -- in respect to given bounds                          */
  GEN.Umin = (DISTR.BD_LEFT  <= -INFINITY) ? 0.0 : CDF(DISTR.BD_LEFT); 
  GEN.Umax = (DISTR.BD_RIGHT >=  INFINITY) ? 1.0 : CDF(DISTR.BD_RIGHT); 

  /* check arguments */
  switch (par->variant) {

  case NINV_VARFLAG_REGULA:

    if ( _FP_same(GEN.s[0], GEN.s[1]) && !GEN.table_on) {
      /* length of interval == 0 -> choose bounderies with                   */
      /*  INTERVAL_COVERS *100 % chance for sign change in interval          */
      GEN.s[0] = -10.;      /* arbitrary starting value                      */
      GEN.s[1] =  10.;      /* arbitrary starting value                      */
      GEN.s[0] = _unur_ninv_regula(gen, (1.-INTERVAL_COVERS)/2. );
      GEN.s[1] = GEN.s[0] + 10.;   /* arbitrary interval length              */
      GEN.s[1] = _unur_ninv_regula(gen, (1.+INTERVAL_COVERS)/2. );
    }

    break;    /* case REGULA end */

  case NINV_VARFLAG_NEWTON:

    if (_FP_same(GEN.s[0], GEN.s[1]) && !GEN.table_on) {
    /* s0 == s1  -> starting value set to value                              */
    /* such that CDF(value) = .5                                             */
      GEN.s[0] = -9.987655;                /* arbitrary starting values      */
      GEN.s[1] =  9.987655;
      GEN.s[0] = _unur_ninv_regula(gen, 0.5);
    }

    break;    /* case NEWTON end */

  }  /* end of switch  */

  
  GEN.CDFmin = -INFINITY;  /* reset in _unur_ninv_create_table()             */
  GEN.CDFmax =  INFINITY;  /* reset in _unur_ninv_create_table()             */
 
  /* generating the table with potential starting values                     */
  if (GEN.table_on){
    /* Umin and Umax are already set */
    _unur_ninv_create_table(gen);
  }
  else{
       GEN.table = NULL;
  }   /* end of if(GEN.table...)                                             */


#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_ninv_debug_init(par,gen);
#endif

  /* free parameters */
  free(par);
  
  return gen;

} /* end of _unur_ninv_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_ninv_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NINV_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NINV_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  switch (par->variant) {
  case NINV_VARFLAG_NEWTON:
    SAMPLE = _unur_ninv_sample_newton;
    break;
  case NINV_VARFLAG_REGULA: default:
    SAMPLE = _unur_ninv_sample_regula;
    break;
  }

  gen->destroy = _unur_ninv_free;
  gen->reinit = NULL;               /* no re-initialization routine required */

  /* copy parameters into generator object */
  GEN.max_iter = PAR.max_iter;      /* maximal number of iterations          */
  GEN.rel_x_resolution = PAR.rel_x_resolution; /* maximal relative error in x*/
  GEN.table_on = PAR.table_on;      /* useage of table for starting points   */
  GEN.table_size = PAR.table_size;  /* number of points for table            */
  GEN.s[0] = PAR.s[0];              /* staring points                        */
  GEN.s[1] = PAR.s[1];

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* init pointer */
  GEN.table = NULL;

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_ninv_create() */

/*****************************************************************************/

double
_unur_ninv_sample_regula( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use regula falsi)                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_regula(gen, _unur_call_urng(gen->urng) ) ;
}

/*---------------------------------------------------------------------------*/

double 
_unur_ninv_regula( struct unur_gen *gen, double u )
     /*---------------------------------------------------------------------*/
     /*   algorithm: regula falsi                                           */
     /*                                                                     */
     /*   parameters:                                                       */
     /*      gen ... pointer to generator object                            */
     /*      u   ... random number (uniform distribution)                   */
     /*   return:                                                           */
     /*     double (sample from random variate)                             */
     /*                                                                     */
     /*   error:                                                            */
     /*     return 0.                                                       */
     /*---------------------------------------------------------------------*/
{ 
    
  double x1, x2, a, xtmp;/* points for RF                                   */
  double x2abs;          /* absolute value of x2                            */
  double f1, f2, ftmp;   /* function values at x1, x2, xtmp                 */
  double length;         /* oriented length of the interval with sign change*/
  double lengthabs;      /* absolute length of interval                     */
  int  lengthsgn;        /* orientation of the Intervalls                   */
  double step;           /* enlarges interval til sign change found         */
  double dx;             /* RF-stepsize                                     */
  int count = 0;         /* counter for  "no sign change"                   */
  int i;                 /* loop variable, index                            */
    

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_NINV_GEN,0.);
  
  if ( _FP_same( GEN.Umin, GEN.Umax) ){
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"CDF constant");
    return INFINITY;   
  }

  /* initialize starting interval */
  if (GEN.table_on){

    /* 0 <= i < table_size-1  */
    i = (int) ( u * 
      ( (GEN.Umax-GEN.Umin)/(GEN.CDFmax-GEN.CDFmin) +
	GEN.Umin - GEN.CDFmin ) *
      (GEN.table_size-1.0) );

    if ( ! _FP_is_minus_infinity(GEN.table[i]) ){
      x1 = GEN.table[i];
      f1 = GEN.CDFmin + i*((GEN.CDFmax-GEN.CDFmin)/(GEN.table_size-1.0));
    }
    else{
      x1 = 2. * GEN.table[i+1] - GEN.table[i+2];
      f1 = CDF(x1);
    }
    if( ! _FP_is_infinity(GEN.table[i+1]) ){
      x2 = GEN.table[i+1];
      f2 = GEN.CDFmin + (i+1)*((GEN.CDFmax-GEN.CDFmin)/(GEN.table_size-1.0));
    }
    else{
      x2 = 2. * GEN.table[i] - GEN.table[i-1];
      f2 = CDF(x2);
    }

  }
  else { /* no table    */

   x1 =  GEN.s[0];      /* left boudary of interval */
   f1 =  CDF(x1);
   x2 =  GEN.s[1];      /* right boudary of interval*/   
   f2 =  CDF(x2);

  }   /* end of if()    */

  if (x1-x2 >= 0) { 
    xtmp = x1; ftmp = f1;
    x1   = x2; f1   = f2;
    x2 = xtmp + DBL_EPSILON;
    f2 = CDF(x2); 
  }

  /* rescale u in respect to given bounds */
  u = GEN.Umin + u * ( GEN.Umax - GEN.Umin );

  /* compute function value at interval boundaries */
  f1 -= u;
  f2 -= u;


  /* search for interval with changing signs */
  step = 1.;     /* interval too small -> make it bigger ( + 2^n * gap ) */ 
  while ( f1*f2 > 0. ) {
    if ( f1 > 0. ) { /* lower boundary too big */    
      x2  = x1;  
      f2  = f1;
      x1 -= step;   
      f1 = CDF(x1) - u;
    }
    else {         /* upper boundary too small */
      x1  = x2;
      f1  = f2;
      x2 += step;
      f2 = CDF(x2) - u;
    }

    /* increase step width */
    step *= 2.;
  }  /* while end -- interval found */ 


  a     = x1;       /* a und x2 soll immer ZW enthalten */

  /* Sekantensuche, ZW wird beibehalten */
  for (i=0; i < GEN.max_iter; i++) {
    count++;
     
    /* f2 always less (better), otherwise change */
    if ( fabs(f1) < fabs(f2) ) {   /* change */
      xtmp = x1; ftmp = f1;
      x1 = x2;   f1 = f2;
      x2 = xtmp; f2 = ftmp;
    }

    x2abs = fabs(x2);   /* absolute value of x2 */

    /* exact hit   || flat region  */    
    if ( f2 == 0. || _FP_same(f1, f2) )
      break; /* -> finished */

    if ( f1*f2 <= 0) {  /* ZeichenWechsel vorhanden             */
      count = 0;   /* zaehler fuer bisektion wird rueckgestellt */
      a    = x1;   /* [a, x2] enthaelt ZW                       */
    }
    
    length = x2 - a;  /* gerichtete laenge */
    lengthabs = fabs(length);
    lengthsgn = (length < 0.) ? -1. : 1.;
    
    if ( lengthabs <= GEN.rel_x_resolution * x2abs  )
      /* relative x-genauigkeit erreicht -> finished */
      break; /* -> finished */

  
    /* Sekanten-Schritt  oder Bisektion */
    dx = ( f1-f2==0. ) ? length/2. : f2*(x2-x1)/(f2-f1) ;  
    
    /* minimaler schritt */
    if ( fabs(dx) < GEN.rel_x_resolution * x2abs ){
      dx = lengthsgn * 0.99 * GEN.rel_x_resolution * x2abs;
      while (x2 == x2 - dx){ /* dx zu klein */
	if ( dx != 2.*dx)    /* am genauigkeits-limit des rechners */
	  dx = 2.*dx;
        else
	  dx = length/2.;    /* Bisektion */
      }
    }

       
    /* Bisektionsschritt, wenn:                             */  
    /* kein  ZW  || Schritt fuhrt aus Intervall             */
    if ( count > 1 || 
        (lengthabs-GEN.rel_x_resolution*x2abs)/(dx*lengthsgn) <= 1. )
      dx = length/2.; /* Bisektionsschritt */
  

    /* Update der Punkte */    
    x1 = x2;       f1 = f2;
    x2 = x2-dx;    f2 = CDF(x2) - u; 
    
  }  /* for-schleife ende */

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (gen->debug & NINV_DEBUG_SAMPLE)
      _unur_ninv_debug_sample_regula( gen,u,x2,f2,i );
#endif
   return x2;

} /* end of _unur_ninv_sample_regula()  */

/*****************************************************************************/

double _unur_ninv_sample_newton(struct unur_gen *gen)
     /*----------------------------------------------------------------------*/
     /* sample from generator (use newtons method)                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*----------------------------------------------------------------------*/
{
  return _unur_ninv_newton(gen, _unur_call_urng(gen->urng) ) ;
}

/*---------------------------------------------------------------------------*/

double
_unur_ninv_newton( struct unur_gen *gen, double U )
     /*----------------------------------------------------------------------*/
     /* sample from generator (use Newton's method)                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*     U   ... random number (uniform distribution)                     */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{ 
  double x;           /* point for netwon-iteration                   */
  double fx;          /* cdf at x                                     */
  double dfx;         /* pdf at x                                     */
  double fxabs;       /* absolute valuo of fx                         */
  double xtmp, fxtmp; /* temprary variables for x and fx              */
  double xold, fxold; /* remember last values for stopping criterion  */
  double fxtmpabs;    /* fabs of fxtmp                                */
  double damp;        /* damping factor                               */
  double step;        /* helps to escape from flat regions of the cdf */
  int i;              /* counter for for-loop, index                  */
    
  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_NINV_GEN,0.);

  if ( _FP_same( GEN.Umin, GEN.Umax) ){
    _unur_warning(gen->genid, UNUR_ERR_GEN_CONDITION,"CDF constant");
    return INFINITY;   
  }

  /* initialize starting point */
  if (GEN.table_on){
    /* i is between 0 and table_size-1 */
    i  = (int) U *
      ( (GEN.Umax-GEN.Umin)/(GEN.CDFmax-GEN.CDFmin) +
         GEN.Umin - GEN.CDFmin ) *
      ( GEN.table_size - 2 );
    x  = GEN.table[i+1];
    fx = GEN.CDFmin + (i+1)*((GEN.CDFmax-GEN.CDFmin)/(GEN.table_size-1.0));
  }
  else{
    x     = GEN.s[0];
    fx    = CDF(x);
  }


  /* rescale u in respect to given bounds */
  U = U*GEN.Umax + (1.0-U)*GEN.Umin;

  fx   -= U;
  dfx   = PDF(x);
  fxabs = fabs(fx);
  xold  = x;    /* there is no old value yet */
  fxold = fx;   /* there is no old value yet */ 

  damp = 2.;        /* to be halved at least once */  
  step = 1.;
  
  /* begin for-loop:  newton-iteration  */
  for (i=0; i < GEN.max_iter; i++) {

    while (dfx == 0.) {   /* function flat at x */
      
      if (fx == 0.)  /* exact hit -> leave while-loopt */
	break; 

      if (fx > 0.)         /* try another x */
        xtmp  = x - step;   
      else 
        xtmp  = x + step;
         
      fxtmp    = CDF(xtmp) - U;
      fxtmpabs = fabs(fxtmp);

      if ( fxtmpabs < fxabs ){        /* improvement, update x            */
        step = 1.;     /* set back stepsize */
        x     = xtmp;
        fx    = fxtmp;
      }
      else if ( fxtmpabs*fxabs < 0. ){/*step was too large, dont update x */
        step /= 2.;                      
      } 
      else{                           /* step was too short, update x     */
        step *= 2.;    
        x     = xtmp;
        fx    = fxtmp;
      }  

      dfx   = PDF(x);
      fxabs = fabs(fx);     
    }   /* end of while-loop, (flat region left) */

   step = 1.;   /* set back stepsize */

   if (fx == 0.)  /* exact hit -> finished */
     break;


    do{    /* newton-step  (damped if nececcary) */
        damp /= 2.;
        xtmp = x - damp * fx/dfx;
        fxtmp = CDF(xtmp) - U;
    } while ( fabs(fxtmp)-fxabs >= fxabs * GEN.rel_x_resolution ); /* no improvement */

    
    /* updation variables according to newton-step      */
    damp  = 2.;       /* set back factor for damping    */
    xold  = x;        /* remember last value of x       */
    fxold = fx;       /* remember last value of fx      */
    x     = xtmp;     /* update x                       */
    fx    = fxtmp;    /* update function value at x     */
    dfx   = PDF(x);   /* update derivative sof fx at x  */
    fxabs = fabs(fx); /* update absolute value of fx    */
 

    /* stopping criterion */
    if ( fabs(x-xold) <= fabs(x) * GEN.rel_x_resolution )        
      break;   /* no improvement with newton-step -> finished */

  }  /* end of for-loop  (MAXITER reached -> finished) */

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file (in case error) */
    if (gen->debug & NINV_DEBUG_SAMPLE)
      _unur_ninv_debug_sample_newton( gen,U,x,fx,i );
#endif


  return x;

} /* end of _unur_ninv_sample_newton() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

void
_unur_ninv_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NINV ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(GEN.table);
  free(gen);

} /* end of _unur_ninv_free() */

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
_unur_ninv_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(par,/*void*/);  COOKIE_CHECK(par,CK_NINV_PAR,/*void*/);
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = ninv (numerical inversion of c.d.f.)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

  fprintf(log,"%s: sampling routine = _unur_ninv_sample",gen->genid);
  switch (par->variant) {
  case NINV_VARFLAG_NEWTON:
    fprintf(log,"_newton\n");
    break;
  case NINV_VARFLAG_REGULA: default:
    fprintf(log,"_regula\n");
    break;
  }

  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_ninv_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_sample_regula( struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (regula falsi)                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d interations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN.max_iter);

} /* end of _unur_ninv_debug_sample_regula() */

/*---------------------------------------------------------------------------*/

void
_unur_ninv_debug_sample_newton( struct unur_gen *gen, double u, double x, double fx, int iter )
     /*----------------------------------------------------------------------*/
     /* trace sampling (newton's method)                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: u = %8.6f,\t x = %8.6g\t(cdf(x)-u = %8.2g)\t -- %2d interations [%d]\n",
	  gen->genid,u,x,fx,iter,GEN.max_iter);

} /* end of _unur_ninv_debug_sample_newton() */

/*---------------------------------------------------------------------------*/

void 
_unur_ninv_debug_chg_domain( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_NINV_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: domain of (truncated) distribution changed:\n",gen->genid);
  fprintf(log,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.BD_LEFT, DISTR.BD_RIGHT);
  fprintf(log,"%s:\tU in (%g,%g)\n",gen->genid,GEN.Umin,GEN.Umax);

} /* end of _unur_ninv_debug_chg_domain() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
