/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      utdr.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection with three points of contact   *
 *              T(x) = -1/sqrt(x)     (T_c with c= -1/2)                     *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and mode of a T-concave distribution                     *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:  pointer to the density, mode of the density                  *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *      double  il, ir       ... left and right boundary of domain           *
 *                               (default: +/- INFINITY)                     *
 *      double  pdf_area     ... area below p.d.f (need not be 1)            *
 *                               (default: 1.)                               *
 *      double *pdf_param    ... parameters of p.d.f.                        *
 *                               (default: NULL)                             *
 *      int     n_pdf_param  ... number of parameters of p.d.f.              *
 *                               (default: 0)                                *
 *      double  c_factor     ... constant for choosing the constr. points    *
 *                               (default: 0.664)                            *
 *      double  delta_factor ... constant for approx. first derivative       *
 *                               (default: 0.00001)                          *
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
 *   [1] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *       (see Algorithm UTDR)                                                *
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

#define UTDR_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define UTDR_SET_CFACTOR        0x001u
#define UTDR_SET_DELTA          0x002u
#define UTDR_SET_PDFMODE        0x004u   /* pdf at mode is set               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "UTDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_utdr_hat( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute hat and squeezes.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_utdr_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Re-initialize (existing) generator.                                       */
/*---------------------------------------------------------------------------*/

static double _unur_utdr_sample( UNUR_GEN *generator );
/** TODO: static double _unur_utdr_sample_check( UNUR_GEN *generator ); **/
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_utdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_utdr_debug_init( struct unur_gen *gen,
				   double try, double trys, double cfac, int setupok);
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       par->data.utdr        /* data for parameter object         */
#define GEN       gen->data.utdr        /* data for generator object         */
#define DISTR     gen->distr.data.cont  /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x) ((*(DISTR.pdf))((x),DISTR.params,DISTR.n_params))    /* call to p.d.f. */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_utdr_new( struct unur_distr *distr )
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
     /*   if the area below the p.d.f. is not close to 1 it is necessary to  */
     /*   set pdf_area to an approximate value of its area (+/- 30 % is ok). */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"p.d.f."); return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); return NULL; }

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_UTDR_PAR);

  /* copy input */
  par->distr       = distr;   /* pointer to distribution object              */

  /* set default values */
  PAR.c_factor     = 0.664; 
          /* optimal value for the normal distribution, which is good for 
	     all bell-shaped densities. The minimax approach for that 
	     transformation has c_factor=2. */
  PAR.delta_factor = 0.00001;
          /* constant for choosing delta to replace the tangent.
	     default should not be changed if used doubles have at least
	     10 decimal digits precision. */

  PAR.fm        = -1.;                /* p.d.f. at mode (unknown)            */
  PAR.hm        = -1.;                /* square of p.d.f. at mode (unknown)  */

  par->method   = UNUR_METH_UTDR;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_utdr_init;

  return par;

} /* end of unur_utdr_new() */

/*****************************************************************************/

int 
unur_utdr_set_pdfatmode( UNUR_PAR *par, double fmode )
     /*----------------------------------------------------------------------*/
     /* Set pdf at mode. if set the p.d.f. at the mode is never changed.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   fmode ... pdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,UTDR );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"pdf(mode)");
    return 0;
  }

  /* store date */
  PAR.fm = fmode;             /* pdf at mode */
  PAR.hm = -1./sqrt(fmode);   /* transformed pdf at mode */

  /* changelog */
  par->set |= UTDR_SET_PDFMODE;

  return 1;

} /* end of unur_utdr_set_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_cfactor( struct unur_par *par, double cfactor )
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
  _unur_check_par_object( par,UTDR );

  /* check new parameter for generator */
  /** TODO: welche werte fuer c sind zulaessig / sinnvoll ? **/
  if (cfactor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c-factor < 0");
    return 0;
  }

  /* store date */
  PAR.c_factor = cfactor;

  /* changelog */
  par->set |= UTDR_SET_CFACTOR;

  return 1;

} /* end of unur_utdr_set_cfactor() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_delta( struct unur_par *par, double delta )
     /*----------------------------------------------------------------------*/
     /* set factor for replacing tangents by secants                         */
     /* (which then are move above the transformed density)                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   delta ... delta-factor                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,par,0 );

  /* check input */
  _unur_check_par_object( par,UTDR );

  /* check new parameter for generator */
  /** TODO: welche werte fuer delta sind zulaessig / sinnvoll ? **/
  if (delta < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"delta < 0");
    return 0;
  }

  /* store date */
  PAR.delta_factor = delta;

  /* changelog */
  par->set |= UTDR_SET_DELTA;

  return 1;

} /* end of unur_utdr_set_delta() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par,UTDR );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | UTDR_VARFLAG_VERIFY) : (par->variant & (~UTDR_VARFLAG_VERIFY));

  /* o.k. */
  return 1;

} /* end of unur_utdr_set_verify() */

/*****************************************************************************/

int
unur_utdr_chg_pdfparams( struct unur_gen *gen, double *params, int n_params )
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
  _unur_check_gen_object( gen,UTDR );
  if (n_params>0) CHECK_NULL(params,0);
  
  /* check new parameter for generator */
  if (n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return 0;
  }

  /* copy parameters */
  DISTR.n_params = n_params;
  for (i=0; i < n_params; i++)
    DISTR.params[i] = params[i];

  /* changelog */
  /* mode and area might be wrong now! 
     but the user is responsible to change it.
     so we dont say:
     gen->distr.set &= ~(UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_PDFAREA );
     gen->set &= ~UTDR_SET_CDFMODE;
  */

  /* o.k. */
  return 1;
} /* end of unur_utdr_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

int 
unur_utdr_chg_domain( struct unur_gen *gen, double left, double right )
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
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,UTDR );

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return 0;
  }

  /* copy new boundaries into generator object */
  DISTR.BD_LEFT = left;
  DISTR.BD_RIGHT = right;
  GEN.il = left;
  GEN.ir = right;

  /* changelog */
  gen->distr.set |= UNUR_DISTR_SET_DOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
#endif
  
  /* o.k. */
  return 1;
  
} /* end of unur_utdr_chg_domain() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_mode( struct unur_gen *gen, double mode )
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
  _unur_check_gen_object( gen,UTDR );
  
  /* copy parameters */
  DISTR.mode = mode;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_utdr_chg_mode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_pdfatmode( struct unur_gen *gen, double fmode )
     /*----------------------------------------------------------------------*/
     /* change value of pdf at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   fmode ... pdf at mode                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,UTDR );

  /* check new parameter for generator */
  if (fmode <= 0.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"pdf(mode)");
    return 0;
  }

  /* store date */
  GEN.fm = fmode;             /* pdf at mode */
  GEN.hm = -1./sqrt(fmode);   /* transformed pdf at mode */

  /* changelog */
  gen->set |= UTDR_SET_PDFMODE;

  /* o.k. */
  return 1;
} /* end of unur_utdr_chg_pdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_utdr_chg_pdfarea( struct unur_gen *gen, double area )
     /*----------------------------------------------------------------------*/
     /* change area below p.d.f. of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   area  ... area                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,0);
  _unur_check_gen_object( gen,UTDR );
  
  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"area <= 0");
    return 0;
  }

  /* copy parameters */
  DISTR.area = area;

  /* no changelog required */

  /* o.k. */
  return 1;
} /* end of unur_utdr_chg_pdfarea() */

/*****************************************************************************/

struct unur_gen *
_unur_utdr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_UTDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_utdr_create(par);
  free(par);
  if (!gen) return NULL;

  /* create hat and squeeze (setup procedure) */
  if ( _unur_utdr_hat(gen) )
    /* hat successfully created */
    return gen;

  else {  /* error */
    _unur_utdr_free(gen);
    return NULL;
  }

} /* end of _unur_utdr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_utdr_hat( struct unur_gen *gen )
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
  double fm;

  int setupok=1;
  double c,cfac,volc,volr,tly,tlys,try,trys,dl,dr,delta,delta1;

  /* check arguments */
  CHECK_NULL( gen, 0 );
  COOKIE_CHECK( gen,CK_UTDR_GEN, 0 );

  /* compute pdf at mode (if not given by user) */
  if (!(gen->set & UTDR_SET_PDFMODE)) {
    fm = PDF(DISTR.mode);
    /* fm must be positive */
    if (fm <= 0.) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"pdf(mode) <= 0.");
      return 0;
    }
    /* step 1.0 of algorithm UTDR */
    GEN.fm = fm;           /* pdf at mode  */
    GEN.hm = -1/sqrt(fm);  /* transformed pdf at mode  */
  }

  /** TODO: inititialisieren notwendig ?? **/
  try = 0.;
  trys = 0.;
  dl = 0.;
  dr = 0.;

  /* start of the set-up procedure */

  /* step 1.0 of algorithm UTDR */
  /* see above or in unur_utdr_set_pdfatmode() */

  do {

    /* 1.1 */
    cfac = (setupok) ? GEN.c_factor : 2.;     /* gibt es hier nur zwei varianten ?? ja*/
    c = cfac * DISTR.area/GEN.fm;
    setupok=1;         

    GEN.tlx = DISTR.mode - c;
    GEN.trx = DISTR.mode + c;

    /* 1.2 */
    if (GEN.il > -INFINITY && GEN.tlx < GEN.il) { 
      GEN.bl = GEN.il;
      GEN.al = 0.;
      GEN.voll = 0.;
      if (GEN.il < DISTR.mode) {
        GEN.tlx = DISTR.mode + (GEN.il - DISTR.mode) * 0.6;
        GEN.sal = (GEN.hm + 1./sqrt(PDF(GEN.tlx))) / (DISTR.mode - GEN.tlx);
      }  
    }
    else {
      tlys = PDF(GEN.tlx);
      if (tlys > 0.) 
	tlys = -1./sqrt(tlys);
      else {
	_unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"f(tlx)=0!!, Perhaps pdf or mode or domain wrong\n");
	goto error;
      }

      GEN.sal =  (GEN.hm - tlys) / (DISTR.mode - GEN.tlx);

      /* da kenn ich mich nicht aus! delta1 > 0; aber was ist mit delta, das is < 0. */
      /* Nein delta1> 0 da -tlys>0 und sal >0; im Fall sal==0 habe ich mich geirrt, da gehoert -tlys */
      delta = ( GEN.sal > 0. ) ? -tlys/GEN.sal : -tlys;
      delta1 = fabs(GEN.tlx);
      delta = GEN.delta_factor * ((delta1<=delta) ? delta : delta1);

      tly = -1./sqrt(PDF(GEN.tlx+delta));
      GEN.al = (tly-tlys)/delta;

      if (GEN.al <= 0.) 
	setupok = 0; /* break ?? soll wenn 1.3 erfolgreich ist, soll das dann
                        nochmals durchlaufen werden?? oder soll 1.3 gar nicht erst 
                        berechnet werden falls dieser Fall hier eintritt??
                        es braucht nichts mehr berechnet werden  */
/*ich mach e hier kein continue, da das zu Endlosschleife fuehren wuerde,sondern
  frage setupok spaeter ab*/
      else {
	GEN.bl = GEN.tlx + (GEN.hm - tly)/GEN.al;

      /* warum wir der rest ausgerechnet wenn setupok auf 0 gesetzt wurde ? */
/*nicht noetig!!*/

        dl = tly - GEN.al * GEN.tlx;
        GEN.voll = -1./(GEN.al * GEN.hm);
        GEN.col = GEN.voll;
        if (GEN.il > -INFINITY)
          GEN.voll += 1./(GEN.al * (GEN.al * GEN.il + dl));
      }
    }

    /* 1.3 */
    if(setupok) {
      if (GEN.ir < INFINITY && GEN.trx > GEN.ir) {
        GEN.br = GEN.ir;
        GEN.ar = 0.;
        volr = 0.;
        if (GEN.ir > DISTR.mode) {
          GEN.trx = DISTR.mode + (GEN.ir - DISTR.mode) * 0.6;
          GEN.sar = (GEN.hm + 1./sqrt(PDF(GEN.trx))) / (DISTR.mode - GEN.trx);
        } 
      }
      else {
        trys = sqrt(PDF(GEN.trx));
        if (trys>0.)
          trys= -1./trys;
        else {
          _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"f(trx)=0!!, Perhaps pdf or mode or domain wrong\n");
	  goto error;
        }
  
        /* siehe 1.2. !! */
        GEN.sar = (GEN.hm - trys) / (DISTR.mode - GEN.trx);
/*delta ist positiv, da trys<0 und sar <0 */
        delta = (GEN.sar<0.) ? trys/GEN.sar : -trys;
        delta1 = fabs(GEN.trx);
        delta = GEN.delta_factor * ((delta1<=delta) ? delta : delta1);
        try = -1./sqrt(PDF(GEN.trx-delta));
        GEN.ar = (trys - try)/delta;
        if (GEN.ar >= 0.) 
          setupok = 0;
        else{ 
          GEN.br = GEN.trx + (GEN.hm - try) / GEN.ar;
          dr = try - GEN.ar * GEN.trx;
          volr = 1./(GEN.ar * GEN.hm);
          GEN.cor = volr;
          if (GEN.ir<INFINITY)
            volr -= 1./(GEN.ar * (GEN.ar * GEN.ir + dr));
        }
      }
    }
    /* 1.4 */
    if(setupok) {
      volc = (GEN.br - GEN.bl) * GEN.fm;
      GEN.vollc = GEN.voll + volc;
      GEN.volcompl = GEN.vollc + volr;
      if (volc>0.) 
        GEN.brblvolc = (GEN.br - GEN.bl)/volc;
      if (GEN.ar!=0.) {
        GEN.drar = dr/GEN.ar;
        GEN.ooar2 = 1./(GEN.ar*GEN.ar);
      }
      if (GEN.al!=0.) {
        GEN.dlal = dl/GEN.al;
        GEN.ooal2 = 1./(GEN.al*GEN.al);
      }
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into log file */
    if (gen->debug) _unur_utdr_debug_init(gen,try,trys,cfac,setupok);
#endif

    if (cfac!=2.) {
      if(setupok)
        if (GEN.volcompl > 4. * DISTR.area || GEN.volcompl < 0.5 * DISTR.area)
        setupok=0;
    }
    else { 
      if (setupok==0 || GEN.volcompl > 8. * DISTR.area || GEN.volcompl < 0.5 * DISTR.area) {
        _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Area below hat too large! Perhaps pdf or mode wrong\n");
	goto error;
      }
    }

  } while (!setupok);

  return 1;

 error:
#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (gen->debug) _unur_utdr_debug_init(gen,try,trys,cfac,setupok);
#endif
  return 0;
} /* end of _unur_utdr_hat() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_utdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_UTDR_GEN);

  /* copy distribution object into generator object */
  memcpy( &(gen->distr), par->distr, sizeof( struct unur_distr ) );

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_utdr_sample;
  gen->destroy = _unur_utdr_free;
  gen->reinit = _unur_utdr_reinit;

  /* copy some parameters into generator object */
  GEN.il = DISTR.BD_LEFT;           /* left boundary of domain         */
  GEN.ir = DISTR.BD_RIGHT;          /* right boundary of domain        */
  GEN.fm = PAR.fm;                  /* pdf at mode                           */
  GEN.hm = PAR.hm;                  /* square root of pdf at mode            */
  GEN.c_factor = PAR.c_factor;
  GEN.delta_factor = PAR.delta_factor;

  /* mode must be in domain */
  DISTR.mode = max(DISTR.mode,GEN.il);
  DISTR.mode = min(DISTR.mode,GEN.ir);

  gen->method = par->method;        /* indicates method                      */
  gen->variant = par->variant;      /* indicates variant                     */
  gen->set = par->set;              /* indicates parameter settings          */
  gen->debug = par->debug;          /* debuging flags                        */
  gen->urng = par->urng;            /* pointer to urng                       */

  gen->urng_aux = NULL;             /* no auxilliary URNG required           */
  gen->gen_aux = NULL;              /* no auxilliary generator objects       */
  gen->gen_aux_2 = NULL;

  /* initialize parameters */
  /** TODO !!! **/
  /** ist das wirklich so noetig ?? **/
  GEN.vollc = 0.; 
  GEN.volcompl = 0.; 
  GEN.voll = 0.; 
  GEN.al = 0.; 
  GEN.ar = 0.; 
  GEN.col = 0.; 
  GEN.cor = 0.; 
  GEN.sal = 0.; 
  GEN.sar = 0.; 
  GEN.bl = 0.; 
  GEN.br = 0.; 
  GEN.tlx = 0.; 
  GEN.trx = 0.; 
  GEN.brblvolc = 0.; 
  GEN.drar = 0.; 
  GEN.dlal = 0.; 
  GEN.ooar2 = 0.; 
  GEN.ooal2 = 0.;
  /* constants of the hat and for generation*/

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_utdr_create() */

/*---------------------------------------------------------------------------*/

int
_unur_utdr_reinit( struct unur_gen *gen )
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
  _unur_check_gen_object( gen,UTDR );

  /* compute universal bounding rectangle */
  return _unur_utdr_hat( gen );
} /* end of _unur_utdr_reinit() */

/*****************************************************************************/

double
_unur_utdr_sample( struct unur_gen *gen )
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
  double u,v,x,help,linx;

  /* check arguments */
  CHECK_NULL(gen,0.);  COOKIE_CHECK(gen,CK_UTDR_GEN,0.);

  while (1) {
    /*2*/
    u = _unur_call_urng(gen) * GEN.volcompl;
    /*2.1*/
    if (u <= GEN.voll) {
      x = -GEN.dlal+GEN.ooal2/(u-GEN.col);
      help = GEN.al*(u-GEN.col);
      linx = help*help;
    }
    else {
      if (u <= GEN.vollc) {
	x = (u-GEN.voll) * GEN.brblvolc + GEN.bl;
	linx = GEN.fm;
      }
      else {
	x = - GEN.drar - GEN.ooar2 / (u-GEN.vollc - GEN.cor);
	help = GEN.ar * (u-GEN.vollc - GEN.cor);
	linx = help*help;
      }
    }
    /*2.2*/
    v = _unur_call_urng(gen) * linx;
    /*2.3*/
    if (x<DISTR.mode) {
      if (x >= GEN.tlx) {
	help = GEN.hm - (DISTR.mode - x) * GEN.sal;
	if (v * help * help <= 1.) return x;
      } 
    }
    else {
      if (x <= GEN.trx) {
	help = GEN.hm - (DISTR.mode - x) * GEN.sar;
	if (v * help * help <= 1.) return x; 
      }
    }
    if (v <= PDF(x)) return x; 
  }

} /* end of _unur_utdr_sample() */

/*****************************************************************************/

void
_unur_utdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_UTDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_UTDR_GEN,/*void*/);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  free(gen);

} /* end of _unur_utdr_free() */

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
_unur_utdr_debug_init( struct unur_gen *gen,
		       double try, double trys, double cfac, int setupok)
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
  CHECK_NULL(gen,/*void*/);  COOKIE_CHECK(gen,CK_UTDR_GEN,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = utdr (transformed density rejection with 3 points of contact)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  _unur_distr_cont_debug( &(gen->distr), gen->genid );

#if 0
  /** TODO **/
  fprintf(log,"%s: sampling routine = _unur_utdr_sample",gen->genid);
  if (gen->variant & UTDR_VARFLAG_VERIFY)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);
#endif

  /** TODO: c_factor, delta_factor **/

  fprintf(log,"%s: tlx %e bl %e mode %e br %e trx %e\n",gen->genid,GEN.tlx,GEN.bl,DISTR.mode,GEN.br,GEN.trx);
  fprintf(log,"%s: try %e trys %e ar %e \n",gen->genid,try,trys,GEN.ar);
  fprintf(log,"%s: cfac %e setupok %d volcompl %e pdf_area %e\n",gen->genid,cfac,setupok,GEN.volcompl,DISTR.area);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_utdr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
