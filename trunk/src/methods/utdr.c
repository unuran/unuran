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

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "UTDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_utdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_utdr_debug_init( struct unur_par *par, struct unur_gen *gen,
				   double try, double trys, double cfac, int setupok);
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.utdr
#define GEN     gen->data.utdr
#define SAMPLE  gen->sample.cont

#define PDF(x) ((*(GEN.pdf))((x),GEN.pdf_param,GEN.n_pdf_param))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_utdr_new( double (*pdf)(double x,double *pdf_param, int n_pdf_param), double mode )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdf  ... probability density function of the desired distribution  */
     /*   mode ... location of the mode of the distribution                  */
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

  /* allocate structure */
  par = _unur_malloc(sizeof(struct unur_par));
  COOKIE_SET(par,CK_UTDR_PAR);

  /* copy input */
  PAR.pdf          = pdf;
  PAR.mode         = mode;

  /* set default values */
  PAR.pdf_area     = 1.;         /* area below p.d.f.                        */
  PAR.il           = -INFINITY;  /* set left boundary of domain              */
  PAR.ir           = INFINITY;   /* set right boundary of domain             */
  PAR.pdf_param    = NULL;       /* no parameters for pdf                    */
  PAR.n_pdf_param  = 0;
  PAR.c_factor     = 0.664; 
          /* optimal value for the normal distribution, which is good for 
	     all bell-shaped densities. The minimax approach for that 
	     transformation has c_factor=2. */
  PAR.delta_factor = 0.00001;
          /* constant for choosing delta to replace the tangent.
	     default should not be changed if used doubles have at least
	     10 decimal digits precision. */

  par->method      = UNUR_METH_UTDR;  /* method and default variant          */
  par->set         = 0UL;          /* inidicate default parameters           */    
  par->urng        = unur_get_default_urng(); /* use default urng            */

  _unur_set_debugflag_default(par);  /* set default debugging flags          */
  _unur_set_genid(par,GENTYPE);      /* set generator identifier             */

  /* routine for starting generator */
  par->init = unur_utdr_init;

  return par;

} /* end of unur_utdr_new() */

/*****************************************************************************/

struct unur_gen *
unur_utdr_init( struct unur_par *par )
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

  int setupok=1;
  double c,cfac,volc,volr,tly,tlys,try,trys,dl,dr,delta,delta1;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_utdr_create(par);
  if (!gen) { free(par); return NULL; }

  /* start of the set-up procedure */

  /* step 1.0 of algorithm UTDR */
  GEN.fm = PDF(GEN.mode);       /* p.d.f. at mode of distribution */
  GEN.hm = -1./sqrt(GEN.fm);    /* transformed p.d.f. at mode     */

  do {

    /* 1.1 */
    cfac = (setupok) ? PAR.c_factor : 2.;     /* gibt es hier nur zwei varianten ?? ja*/
    c = cfac * PAR.pdf_area/GEN.fm;
    setupok=1;         

    GEN.tlx = GEN.mode - c;
    GEN.trx = GEN.mode + c;

    /* 1.2 */
    if (GEN.il > -INFINITY && GEN.tlx < GEN.il) { 
      GEN.bl = GEN.il;
      GEN.al = 0.;
      GEN.voll = 0.;
      if (GEN.il < GEN.mode) {
        GEN.tlx = GEN.mode + (GEN.il - GEN.mode) * 0.6;
        GEN.sal = (GEN.hm + 1./sqrt(PDF(GEN.tlx))) / (GEN.mode - GEN.tlx);
      }  
    }
    else {
      tlys = PDF(GEN.tlx);
      if (tlys > 0.) 
	tlys = -1./sqrt(tlys);
      else {
	_unur_warning(gen->genid,UNUR_ERR_INIT,"f(tlx)=0!!, Perhaps pdf or mode or domain wrong\n");
	free(par); unur_utdr_free(gen);
	return NULL; 
      }

      GEN.sal =  (GEN.hm - tlys) / (GEN.mode - GEN.tlx);

      /* da kenn ich mich nicht aus! delta1 > 0; aber was ist mit delta, das is < 0. */
      /* Nein delta1> 0 da -tlys>0 und sal >0; im Fall sal==0 habe ich mich geirrt, da gehoert -tlys */
      delta = ( GEN.sal > 0. ) ? -tlys/GEN.sal : -tlys;
      delta1 = fabs(GEN.tlx);
      delta = PAR.delta_factor * ((delta1<=delta) ? delta : delta1);

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
        if (GEN.ir > GEN.mode) {
          GEN.trx = GEN.mode + (GEN.ir - GEN.mode) * 0.6;
          GEN.sar = (GEN.hm + 1./sqrt(PDF(GEN.trx))) / (GEN.mode - GEN.trx);
        } 
      }
      else {
        trys = sqrt(PDF(GEN.trx));
        if (trys>0.)
          trys= -1./trys;
        else {
          _unur_warning(gen->genid,UNUR_ERR_INIT,"f(trx)=0!!, Perhaps pdf or mode or domain wrong\n");
          free(par); unur_utdr_free(gen);
          return NULL;
        }
  
        /* siehe 1.2. !! */
        GEN.sar = (GEN.hm - trys) / (GEN.mode - GEN.trx);
/*delta ist positiv, da trys<0 und sar <0 */
        delta = (GEN.sar<0.) ? trys/GEN.sar : -trys;
        delta1 = fabs(GEN.trx);
        delta = PAR.delta_factor * ((delta1<=delta) ? delta : delta1);
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

#if UNUR_DEBUG & UNUR_DB_INFO
    /* write info into log file */
    if (gen->debug) _unur_utdr_debug_init(par,gen,try,trys,cfac,setupok);
#endif

    if (cfac!=2.) {
      if(setupok)
        if (GEN.volcompl > 4. * PAR.pdf_area||GEN.volcompl<0.5*PAR.pdf_area)
        setupok=0;
    }
    else { 
      if (setupok==0 || GEN.volcompl > 8. * PAR.pdf_area||GEN.volcompl<0.5*PAR.pdf_area) {
        _unur_warning(gen->genid,UNUR_ERR_INIT,"Area below hat too large! Perhaps pdf or mode wrong\n");
        free(par); unur_utdr_free(gen);
        return NULL; 
      }
    }


  } while (!setupok);

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_utdr_init() */

/*****************************************************************************/

double
unur_utdr_sample( struct unur_gen *gen )
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
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_UTDR_GEN,0.);

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
    if (x<GEN.mode) {
      if (x >= GEN.tlx) {
	help = GEN.hm - (GEN.mode - x) * GEN.sal;
	if (v * help * help <= 1.) return x;
      } 
    }
    else {
      if (x <= GEN.trx) {
	help = GEN.hm - (GEN.mode - x) * GEN.sar;
	if (v * help * help <= 1.) return x; 
      }
    }
    if (v <= PDF(x)) return x; 
  }

} /* end of unur_utdr_sample() */

/*****************************************************************************/

void
unur_utdr_free( struct unur_gen *gen )
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

  /* magic cookies */
  COOKIE_CHECK(gen,CK_UTDR_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_free_genid(gen);
  if (GEN.pdf_param) free(GEN.pdf_param);
  free(gen);

} /* end of unur_utdr_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

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
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_UTDR_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_UTDR_GEN);

  /* routines for sampling and destroying generator */
  SAMPLE = unur_utdr_sample;
  gen->destroy = unur_utdr_free;

  /* set all pointers to NULL */
  GEN.pdf_param = NULL;
  GEN.n_pdf_param = 0;

  /* copy some parameters into generator object */
  _unur_copy_urng_pointer(par,gen);  /* pointer to urng into generator object*/
  _unur_copy_debugflag(par,gen);     /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);         /* copy generator identifier            */

  GEN.mode = PAR.mode;         /* mode of distribution     */
  GEN.il   = PAR.il;           /* left boundary of domain  */
  GEN.ir   = PAR.ir;           /* left boundary of domain  */
  GEN.pdf  = PAR.pdf;          /* p.d.f. of distribution   */
  gen->method = par->method;   /* indicates method and variant */

  /* allocation for array of pdf-parameters */
  if( PAR.n_pdf_param > 0 )
    GEN.pdf_param = _unur_malloc( PAR.n_pdf_param * sizeof(double) );
  
  /* copy parameters of distribution */
  GEN.n_pdf_param = PAR.n_pdf_param;
  if (PAR.pdf_param != NULL)
    for (i=0; i<PAR.n_pdf_param; i++)
      GEN.pdf_param[i] = PAR.pdf_param[i];

  /* return pointer to (almost empty) generator object */
  return(gen);
  
} /* end of _unur_utdr_create() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

static void
_unur_utdr_debug_init( struct unur_par *par, struct unur_gen *gen,
		       double try, double trys, double cfac, int setupok)
     /*----------------------------------------------------------------------*/
     /* write info about generator into logfile                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  log = unur_get_stream();

  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = utdr (transformed density rejection with 3 points of contact)\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: sampling routine = unur_utdr_sample",gen->genid);
  if (par->method & UNUR_MASK_SCHECK)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  fprintf(log,"%s:\n",gen->genid);

  fprintf(log,"%s: tlx %e bl %e mode %e br %e trx %e\n",gen->genid,GEN.tlx,GEN.bl,GEN.mode,GEN.br,GEN.trx);
  fprintf(log,"%s: try %e trys %e ar %e \n",gen->genid,try,trys,GEN.ar);
  fprintf(log,"%s: cfac %e setupok %d volcompl %e pdf_area %e\n",gen->genid,cfac,setupok,GEN.volcompl,PAR.pdf_area);
  fprintf(log,"%s:\n",gen->genid);

} /* end of _unur_utdr_debug_init() */

#endif

/*****************************************************************************/

