/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given p.d.f and .... of a T-concave distribution                     *
 *      produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:  pointer to the density, ....                                 *
 *                                                                           *
 *   PARAMETERS:                                                             *
 *      double *pdf_param    ... parameters of p.d.f.                        *
 *                               (default: NULL)                             *
 *      int     n_pdf_param  ... number of parameters of p.d.f.              *
 *                               (default: 0)                                *
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Ahrens J. H. (1993): Sampling from general distributions by         *
 *       suboptimal division of domains,                                     *
 *       Grazer Math. Berichte 319, 30pp.                                    *
 *                                                                           *
 *   [2] Ahrens J. H. (1995): An one-table method for sampling from          *
 *       continuous and discrete distributions,                              *
 *       Computing 54(2), pp. 127-146                                        *
 *                                                                           *
 *   SEE ALSO:                                                               *
 *   [3] Gilks, W. R. and Wild,  P. (1992):                                  *
 *       Adaptive rejection sampling for Gibbs sampling,                     *
 *       Applied Statistics 41, pp. 337-348                                  *
 *                                                                           *
 *   [4] Zaman, A. (1996), Generation of Random Numbers from an Arbitrary    *
 *       Unimodal Density by Cutting Corners, unpublished manuskript         *
 *       available at http://chenab.lums.edu.pk/~arifz/                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * ..... beschreibung ....                                                   *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Tue Sep 21 10:03:03 CEST 1999                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 1999 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *                                                                           *
 *   This library is free software; you can redistribute it and/or           *
 *   modify it under the terms of the GNU Library General Public             *
 *   License as published by the Free Software Foundation; either            *
 *   version 2 of the License, or (at your option) any later version.        *
 *                                                                           *
 *   This library is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 *   Library General Public License for more details.                        *
 *                                                                           *
 *   You should have received a copy of the GNU Library General Public       *
 *   License along with this library; if not, write to the Free              *
 *   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "TABL"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals.                                               */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_starting_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals, use given slopes                              */
/*---------------------------------------------------------------------------*/
static int _unur_tabl_get_starting_intervals_from_mode( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals, use mode and domain                           */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_split_a_starting_intervals( struct unur_par *par, struct unur_gen *gen, struct unur_tabl_interval *iv_slope );
/*---------------------------------------------------------------------------*/
/* split starting intervals according to [1]                                 */
/* SPLIT A (equal areas rule)                                                */
/*---------------------------------------------------------------------------*/

static int 
_unur_tabl_split_b_starting_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* split starting intervals according to [1]                                 */
/* SPLIT B, but instead of the iteration in [1] use "arcmean".               */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_split_interval( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			   double x, double fx, unsigned int split_mode );
/*---------------------------------------------------------------------------*/
/* split interval (replace old one by two new ones in same place)            */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *_unur_tabl_iv_stack_pop( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* pop an interval from the stack of free intervals.                         */
/*---------------------------------------------------------------------------*/

#if 0 /* we do not need this subroutine yet */
static void _unur_tabl_iv_stack_push( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* push the last popped interval back onto the stack.                         */
/*---------------------------------------------------------------------------*/
#endif

#if UNUR_DEBUG & UNUR_DB_INFO
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_intervals( struct unur_gen *gen, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals.                                                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR     par->data.tabl
#define GEN     gen->data.tabl
#define SAMPLE  gen->sample.cont

#define PDF(x) ((*(GEN.pdf))((x),GEN.pdf_param,GEN.n_pdf_param))

/*---------------------------------------------------------------------------*/

/* indicate variant of method */

#define TABL_MASK_SPLIT    0x003u   /* indicates how to split interval */
/*                                     split at        computation     convergence of hat */
#define TABL_SPLIT_POINT   0x000u   /* sampled point    none            slowest          */
#define TABL_SPLIT_MEAN    0x001u   /* mean point       slower          better           */
#define TABL_SPLIT_ARC     0x002u   /* "arcmean"        very slow       very good for almost unbounded domain */

#define TABL_MASK_STP      0x0f0u   /* indicates if starting intervals have to be split */
#define TABL_STP_SPLIT_A   0x010u   /* use equal area rule (SPLIT A in [1])   */
#define TABL_STP_SPLIT_B   0x020u   /* use main subdivisions (SPLIT B in [1]) */

/* Special debugging flags (do not use the first 3 bits) */
#define TABL_DB_IV         0x01u    /* show intervals                        */
#define TABL_DB_A_IV       0x02u    /* show intervals after split A, before split B */

/*---------------------------------------------------------------------------*/

#define min(x,y)   (((x)<(y)) ? (x) : (y))
#define max(x,y)   (((x)>(y)) ? (x) : (y))

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  User Interface                                                         **/
/*****************************************************************************/

struct unur_par *
unur_tabl_new( double (*pdf)(double x,double *pdf_param, int n_pdf_param) )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdf  ... probability density function of the desired distribution  */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* allocate structure */
  par = _unur_malloc( sizeof(struct unur_par) );
  COOKIE_SET(par,CK_TABL_PAR);

  /* copy input */
  PAR.pdf           = pdf;

  /* set default values */
  PAR.pdf_param     = NULL;      /* no parameters for pdf                    */
  PAR.n_pdf_param   = 0;         /* number of parameters                     */
  PAR.pdf_area      = 1.;        /* area below p.d.f.                        */

  PAR.slopes        = NULL;      /* pointer to slopes of p.d.f.              */
  PAR.n_slopes      = 0;         /* number of slopes                         */

  PAR.mode          = 0.;        /* mode of p.d.f.                           */
  PAR.bleft         = 0.;        /* left boundary of domain (no useful default) */
  PAR.bright        = 0.;        /* right boundary of domain (left = right --> cannot make hat) */
  PAR.mode          = 0.;        /* (exact!) location of mode                */

  PAR.n_starting_cpoints = 30;   /* number of starting points                */
  PAR.area_fract    = 0.25;      /* parameter for equal area rule (default from [1] ) */

  PAR.max_ivs       = 100;       /* maximum number of intervals              */
  PAR.max_ratio     = 0.95;      /* bound for ratio  Atotal / Asqueeze       */

  PAR.guide_factor  = 1.; /* guide table has same size as array of intervals */

  par->method       = UNUR_METH_TABL;         /* indicate method             */
  PAR.variant       = (TABL_SPLIT_ARC   |     /* variant: split at arc_mean  */
		       TABL_STP_SPLIT_A |     /* run SPLIT A on slopes       */
		       TABL_STP_SPLIT_B  );   /* run SPLIT B on slopes       */

  par->set          = 0UL;       /* inidicate default parameters             */    
  par->urng         = unur_get_default_urng(); /* use default urng           */

  _unur_set_debugflag_default(par); /* set default debugging flags           */
  _unur_set_genid(par,GENTYPE);     /* set generator identifier              */

  /* routine for starting generator */
  par->init = unur_tabl_init;

  return par;

} /* end of unur_tabl_new() */

/*****************************************************************************/

int 
unur_set_tabl_variant( struct unur_par *par, unsigned int variant )
     /*----------------------------------------------------------------------*/
     /* set variant of method                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   variant ... indicator for variant                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);
  COOKIE_CHECK(par,CK_TABL_PAR,0);

  /* changelog */
  par->set |= UNUR_SET_VARIANT;

  PAR.variant = variant;

  return 1;

} /* end if unur_set_tabl_variant() */

/*****************************************************************************/

struct unur_gen *
unur_tabl_init( struct unur_par *par )
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
  CHECK_NULL(par,NULL);
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_tabl_create(par);
  if (!gen) { free(par); return NULL; }

  /* get slopes for starting generator */
  if (!_unur_tabl_get_starting_intervals(par,gen)) {
    _unur_error(gen->genid,UNUR_ERR_INIT,"Cannot make hat function.");
    free(par); unur_tabl_free(gen);
    return NULL;
  }

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug) _unur_tabl_debug_init(par,gen);
  if (gen->debug & TABL_DB_A_IV)
    _unur_tabl_debug_intervals(gen,FALSE);
#endif

  /* split according to [1], run SPLIT B */
  if (PAR.variant & TABL_STP_SPLIT_B)
    while (GEN.n_ivs < PAR.n_starting_cpoints)
      if (!_unur_tabl_split_b_starting_intervals(par,gen))
	return NULL;
  
  /* we have to update the maximal number of intervals,
     if the user wants more starting points. */
  if (GEN.n_ivs > GEN.max_ivs) {
    /** TODO: do not allow too many intervals ?? **/
    _unur_warning(gen->genid,UNUR_ERR_INIT,"maximal number of intervals too small. increase.");
    GEN.max_ivs = GEN.n_ivs;
  }

  /* make initial guide table */
  _unur_tabl_make_guide_table(gen);

#if UNUR_DEBUG & UNUR_DB_INFO
  /* write info into log file */
  if (gen->debug)
    _unur_tabl_debug_intervals(gen,TRUE);
#endif

  /* free parameters */
  free(par);

  return gen;

} /* end of unur_tabl_init() */

/*****************************************************************************/

double
unur_tabl_sample_adaptive( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    do {
      /* since we do not upgrade the guide table every time we 
	 split an interval, there exists "black intervals".
	 The total area below the hat decreases, but the values 
	 of iv->Acum cannot be upgraded every time.
	 So there arise holes, i.e., these "black intervals",
	 that represent area above the hat function (cut off by
	 the splitting process).
	 if we hit such a hole, we have to reject the generated 
	 interval, and try again. 
      */

      /* sample from U(0,1) */
      u = _unur_call_urng(gen);
    
      /* look up in guide table and search for interval */
      iv =  GEN.guide[(int) (u * GEN.guide_size)];
      COOKIE_CHECK(iv,CK_TABL_IV,0.);
      u *= GEN.Atotal;
      while (iv->Acum < u) {
	if (iv->next == NULL)
	  /* we have hit an imaginary "black interval" at the end of the list. try again. */
	  break;
	iv = iv->next;
	COOKIE_CHECK(iv,CK_TABL_IV,0.);
      }
      
      /* reuse of uniform random number */
      u = iv->Acum - u;
    
    } while (u > iv->Ahat); /* check whether we have hit a "black interval" */

    /* generation w.r.t. squeeze should be inversion */
    if (iv->slope>0)
      u = iv->Ahat - u;

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
      /** TODO: possible overflow/underflow ?? **/
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (GEN.variant & TABL_MASK_SPLIT) );
  	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample_adaptive() */

/*****************************************************************************/

double
unur_tabl_sample( struct unur_gen *gen )
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for interval */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,0.);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      return( iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze ); 
      /** TODO: possible overflow/underflow ?? **/
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (GEN.variant & TABL_MASK_SPLIT) );
	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample() */

/*****************************************************************************/

double
unur_tabl_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
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
  struct unur_tabl_interval *iv;
  double u,x,fx;

  /* check arguments */
  CHECK_NULL(gen,0.);
  COOKIE_CHECK(gen,CK_TABL_GEN,0.);

  while(1) {

    /* sample from U(0,1) */
    u = _unur_call_urng(gen);

    /* look up in guide table and search for interval */
    iv =  GEN.guide[(int) (u * GEN.guide_size)];
    u *= GEN.Atotal;
    while (iv->Acum < u)
      iv = iv->next;

    COOKIE_CHECK(iv,CK_TABL_IV,0.);

    /* reuse of uniform random number
       (generation of squeeze should be inversion) */
    u = (iv->slope<0) ? (iv->Acum - u) : (iv->Ahat + u - iv->Acum);

    if( u <= iv->Asqueeze ) {
      /* below squeeze */
      x = iv->xmax + (iv->Asqueeze-u) * (iv->xmin - iv->xmax)/iv->Asqueeze;
      /** TODO: possible overflow/underflow ?? **/
      /* test whether p.d.f. is monotone */
      fx = PDF(x);
      if (fx > iv->fmax)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf > hat. pdf not monotone in interval");
      if (fx < iv->fmin)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf < squeeze. pdf not monotone in interval");
      /* at last return number */
      return x;
    }
    else {
      /* between spueeze and hat --> have to valuate p.d.f. */
      x = iv->xmax + (u-iv->Asqueeze) * (iv->xmin - iv->xmax)/(iv->Ahat - iv->Asqueeze);
      /** TODO: possible overflow/underflow ?? **/
      fx = PDF(x);
      /* test whether p.d.f. is monotone */
      if (fx > iv->fmax)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf > hat. pdf not monotone in interval");
      if (fx < iv->fmin)
	_unur_warning(gen->genid,UNUR_ERR_SAMPLE,"pdf < squeeze. pdf not monotone in interval");
      /* split interval */
      if (GEN.n_ivs < GEN.max_ivs && GEN.max_ratio * GEN.Atotal > GEN.Asqueeze) {
      	_unur_tabl_split_interval( gen, iv, x, fx, (GEN.variant & TABL_MASK_SPLIT) );
	_unur_tabl_make_guide_table(gen);
	/** TODO: it is not necessary to update the guide table every time. 
	    But then (1) some additional bookkeeping is required and
	    (2) the guide table method requires a acc./rej. step. **/
      }
  
      /* now accept or reject */
      u = _unur_call_urng(gen);
      if (fx >= u * (iv->fmax - iv->fmin) + iv->fmin)
	/** TODO: possible overflow/underflow ?? **/
	return x;
    }
  }

} /* end of unur_tabl_sample_check() */

/*****************************************************************************/

void
unur_tabl_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);
  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* write info into log file */
#if UNUR_DEBUG & UNUR_DB_INFO
  if (gen->debug) _unur_tabl_debug_free(gen);
#endif

  /* free linked list of intervals and others */
  _unur_free_mblocks(GEN.mblocks);

  /* free other memory */
  _unur_free_genid(gen);
  if (GEN.pdf_param) free(GEN.pdf_param);
  free(GEN.guide);
  free(gen);

} /* end of unur_tabl_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

static struct unur_gen *
_unur_tabl_create( struct unur_par *par )
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
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);

  /* allocate memory for generator object */
  gen = _unur_malloc( sizeof(struct unur_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_TABL_GEN);

  /* routines for sampling and destroying generator */
  SAMPLE = (par->method & UNUR_MASK_SCHECK) ? unur_tabl_sample_check : unur_tabl_sample;
  gen->destroy = unur_tabl_free;

  /* set all pointers to NULL */
  GEN.pdf_param   = NULL;
  GEN.n_pdf_param = 0;
  GEN.Atotal      = 0.;
  GEN.Asqueeze    = 0.;
  GEN.guide       = NULL;
  GEN.guide_size  = 0;
  GEN.iv          = NULL;
  GEN.n_ivs       = 0;
  GEN.iv_stack    = NULL;
  GEN.iv_free     = 0;
  GEN.mblocks     = NULL;

  /* copy some parameters into generator object */
  GEN.pdf         = PAR.pdf;           /* p.d.f. of distribution                */
  GEN.bleft       = PAR.bleft;         /* left boundary of domain               */
  GEN.bright      = PAR.bright;        /* right boundary of domain              */

  GEN.guide_factor = PAR.guide_factor; /* relative size of guide tables         */

  /* bounds for adding construction points  */
  GEN.max_ivs   = PAR.max_ivs;         /* maximum number of intervals           */
  GEN.max_ratio = PAR.max_ratio;       /* bound for ratio  Atotal / Asqueeze    */

  gen->method = par->method;           /* indicates method                      */
  GEN.variant = PAR.variant;           /* indicates variant                     */

  _unur_copy_urng_pointer(par,gen);    /* pointer to urng into generator object */
  _unur_copy_debugflag(par,gen);       /* copy debugging flags into generator object */
  _unur_copy_genid(par,gen);           /* copy generator identifier             */

  /* allocate memory for parameters of p.d.f. */
  if( PAR.n_pdf_param > 0 )
    GEN.pdf_param = _unur_malloc( PAR.n_pdf_param * sizeof(double) );

  /* copy parameters of distribution */
  GEN.n_pdf_param = PAR.n_pdf_param;
  if (PAR.pdf_param != NULL)
    for (i=0; i<PAR.n_pdf_param; i++)
      GEN.pdf_param[i] = PAR.pdf_param[i];

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_tabl_create() */

/*****************************************************************************/

static int
_unur_tabl_get_starting_intervals( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   pdf(a) >= pdf(b)                                                   */
     /*----------------------------------------------------------------------*/
{

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* we have two cases: 
     (1) we are given slopes --> check these, compute domain if necessary
     (2) we are given domain and mode --> compute slopes */

  if (PAR.n_slopes > 0 )
    /* slopes are given */
    return _unur_tabl_get_starting_intervals_from_slopes(par,gen);

  if ( (par->set & UNUR_SET_DOMAIN) && (par->set & UNUR_SET_MODE) )
    /* no slopes given. need domain and mode */
    /* compute slopes */
    return _unur_tabl_get_starting_intervals_from_mode(par,gen);

  /* else */
  _unur_error(gen->genid,UNUR_ERR_INIT,"number of slopes <= 0, domain or mode not given.");
  return 0;

} /* end of _unur_tabl_get_starting_intervals() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_get_starting_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals, slopes are given by user.                */
     /* estimate domain when not given.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*                                                                      */
     /* comment:                                                             */
     /*   a slope <a,b> is an interval [a,b] or [b,a] such that              */
     /*   pdf(a) >= pdf(b)                                                   */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  int i;

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* init counter of intervals */
  GEN.n_ivs = 0;
  iv = NULL;

  /* compute initial intervals */
  for ( i=0; i < 2*PAR.n_slopes; i+=2 ) {
    /* get a new interval and link into list */
    if (i==0)
      iv = GEN.iv = _unur_tabl_iv_stack_pop(gen);    /* the first interval */
    else
      iv = iv->next = _unur_tabl_iv_stack_pop(gen);  /* all the other intervals */
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    /* max and min of p.d.f. in interval */
    iv->xmax = PAR.slopes[i];      
    iv->fmax = PDF(iv->xmax);
    iv->xmin = PAR.slopes[i+1];    
    iv->fmin = PDF(iv->xmin);

    /* check slopes */
    if (iv->fmax < iv->fmin) {
      _unur_error(gen->genid,UNUR_ERR_INIT,"slopes non-decreasing");
      return 0;
    }

    /* area of slope and sign of slope (increasing/decreasing) */
    iv->slope = (iv->xmax > iv->xmin) ? 1 : -1;
    iv->Ahat = iv->slope * (iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = iv->slope * (iv->xmax - iv->xmin) * iv->fmin;
    /** TODO: possible overflow/underflow ?? **/
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* estimate domain */
    if (!(par->set & UNUR_SET_DOMAIN)) {
      if (iv->slope > 0) {
	GEN.bleft = min(GEN.bleft,iv->xmin);
	GEN.bright = max(GEN.bright,iv->xmax);
      }
      else {
	GEN.bleft = min(GEN.bleft,iv->xmax);
	GEN.bright = max(GEN.bright,iv->xmin);
      }
    }

    /* split interval following [1], split A */
    if (PAR.variant & TABL_STP_SPLIT_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return 0;
    }
  }

  /* terminate list */
  iv->next = NULL;

  /* o.k. */
  return 1;

} /* end of _unur_tabl_get_starting_intervals_from_slopes() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_get_starting_intervals_from_mode( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute starting intervals                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par          ... pointer to parameter list                         */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;

  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* init linked list of intervals */
  GEN.n_ivs = 0;

  /* compute initial intervals */
  while (1) {
    /* the first interval */
    iv = GEN.iv = _unur_tabl_iv_stack_pop(gen);
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    if (PAR.mode <= PAR.bleft) {
      /* only one ascending interval <a,b> = [a,b] */
      iv->xmax = PAR.mode;
      iv->xmin = PAR.bright;
      break;
    }

    if (PAR.mode >= PAR.bright) {
      /* only one descending interval <a,b> = [b,a] */
      iv->xmax = PAR.mode;
      iv->xmin = PAR.bleft;
      break;
    }

    /* one descending and one ascending interval */
    iv->xmax = PAR.mode;
    iv->xmin = PAR.bleft;

    /* the second interval */
    iv = iv->next = _unur_tabl_iv_stack_pop(gen);  /* all the other intervals */
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    iv->xmax = PAR.mode;
    iv->xmin = PAR.bright;
    break;
  }

  /* terminate list */
  iv->next = NULL;

  /* compute parameters */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);

    /* max and min of p.d.f. in interval */
    iv->fmax = PDF(iv->xmax);
    iv->fmin = PDF(iv->xmin);

    /* area of slope and sign of slope (increasing/decreasing) */
    iv->slope = (iv->xmax > iv->xmin) ? 1 : -1;
    iv->Ahat = iv->slope * (iv->xmax - iv->xmin) * iv->fmax;
    /** TODO: possible overflow/underflow ?? **/
    iv->Asqueeze = iv->slope * (iv->xmax - iv->xmin) * iv->fmin;
    /** TODO: possible overflow/underflow ?? **/
    /* avoid strange (possible) floating point execption on non IEEE754 architecture */
    iv->Acum = 0.;

    /* split interval following [1], split A */
    if (PAR.variant & TABL_STP_SPLIT_A) {
      iv = _unur_tabl_split_a_starting_intervals( par, gen, iv );
      if (iv == NULL) return 0;
    }

  }

  /* o.k. */
  return 1;

} /* end of _unur_tabl_get_starting_intervals_from_mode() */

/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_split_a_starting_intervals( struct unur_par *par, 
				       struct unur_gen *gen, 
				       struct unur_tabl_interval *iv_slope )
     /*----------------------------------------------------------------------*/
     /* split starting intervals according to [1]                            */
     /* SPLIT A (equal areas rule)                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter list                            */
     /*   gen       ... pointer to generator object                          */
     /*   iv_slope  ... pointer to interval of slope                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to last interval in list of splitted slope                 */
     /*   NULL on error                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv, *iv_last, *iv_new;
  double bar_area, x;
  
  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,NULL);
  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  COOKIE_CHECK(iv_slope,CK_TABL_IV,NULL);
  if (iv_slope->slope != 1 && iv_slope->slope != -1 ) {
    _unur_warning( gen->genid, UNUR_ERR_INIT, "invalid slope.");
    return iv_slope;
  }

  iv = iv_slope;        /* pointer to actual interval */
  iv_last = iv_slope;   /* pointer to last interval in list */
  /* (maximal) area of bar (= hat in one interval) */
  bar_area = PAR.pdf_area * PAR.area_fract;

  while (iv->Ahat > bar_area) {

    switch (iv->slope) {
    case +1:
      /* move from right to left */
      x = iv->xmax - bar_area / iv->fmax;
      iv_new = _unur_tabl_split_interval( gen, iv, x, PDF(x), TABL_SPLIT_POINT );
      if (iv_last == iv_slope)
	iv_last = iv_new;
      break;
    case -1:
      /* move from left to right */
      x = iv->xmax + bar_area / iv->fmax;
      iv = _unur_tabl_split_interval( gen, iv, x, PDF(x), TABL_SPLIT_POINT );
      break;
    }
  }

  /* pointer to last interval */
  return ((iv->slope < 0) ? iv : iv_last);

} /* end of _unur_tabl_split_a_starting_intervals() */

/*---------------------------------------------------------------------------*/

static int
_unur_tabl_split_b_starting_intervals( struct unur_par *par, 
				       struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* split starting intervals according to [1]                            */
     /* SPLIT B, but instead of the iteration in [1] use "arcmean".          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter list                                  */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... error                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Amean;  /* mean area between hat and squeeze in slope */
  
  /* check arguments */
  COOKIE_CHECK(par,CK_TABL_PAR,0);
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* compute mean area between squeeze and hat for each interval */
  Amean = 0.;
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    Amean += iv->Ahat - iv->Asqueeze;
  }
  Amean /= GEN.n_ivs;

  /* now split intervals */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    if ((iv->Ahat - iv->Asqueeze) >= Amean) {
      /* new point instead of the interation of [1] we use "arcmean" */
      iv = _unur_tabl_split_interval( gen, iv, 0., 0., TABL_SPLIT_ARC );
      if (GEN.n_ivs >= PAR.n_starting_cpoints)
	/* no more intervals, yet */
	break;
    }
  }

  /* o.k. */
  return 1;

} /* end of _unur_tabl_split_b_starting_intervals() */

/*****************************************************************************/

static struct unur_tabl_interval *
_unur_tabl_split_interval( struct unur_gen *gen,
				struct unur_tabl_interval *iv_old, 
				double x, double fx, 
				unsigned int split_mode )
     /*----------------------------------------------------------------------*/
     /* split interval (replace old one by two new ones in same place)       */
     /* new interval is inserted immedately after old one.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval that has to be split                   */
     /*   x   ... splitting point                                            */
     /*   fx  ... value of p.d.f. at splitting point                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new (right) interval                                    */
     /*   NULL on error                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv_new;

  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);
  COOKIE_CHECK(iv_old,CK_TABL_IV,NULL);

  /* There are three possibilities for the splitting point:
     (1) use x and avoid computation of pdf(x). 
     (2) use middle of interval. converges faster in many cases.
     (3) use "arc_mean" of interval. 
         converges faster when domain is almost unbounded. */
  switch( split_mode ) {
  case TABL_SPLIT_POINT:    /* (1) */
    /* nothing to do (default) */
    break;
  case TABL_SPLIT_MEAN:     /* (2) */
    x = 0.5 * (iv_old->xmin + iv_old->xmax); 
    fx = PDF(x);
    break;
  case TABL_SPLIT_ARC:      /* (3) */
    x = _unur_arcmean(iv_old->xmin, iv_old->xmax); 
    fx = PDF(x);
    break;
  default: 
    /* this should not happen */
    _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"Invalid variant, use default");
    break;
  }

  /* we need a new interval */
  iv_new = _unur_tabl_iv_stack_pop(gen);
  COOKIE_CHECK(iv_new,CK_TABL_IV,0);

  /* iv_new has the same slope as iv_old */
  iv_new->slope = iv_old->slope;

  /* we have to distinguish between two cases:
     pdf is increasing (slope = +1) or
     pdf is decreasing (slope = -1). */

  switch (iv_old->slope) {
  case -1:
    /* (x) The iv_new inherits the minimum of iv_old.
           iv_old keeps the maximum.
       (x) The splitting point is the maximum of iv_new and
           the minimum of iv_old.
    */
    iv_new->xmin  = iv_old->xmin;  
    iv_new->fmin = iv_old->fmin;
    iv_old->xmin  = iv_new->xmax = x; 
    iv_old->fmin = iv_new->fmax = fx; 
    break;
  case +1: /* the other way round */
    /* (x) The iv_new inherits the maximum of iv_old.
           iv_old keeps the minimum.
       (x) The splitting point is the minimum of iv_new and
           the maximum of iv_old.
    */
    iv_new->xmax  = iv_old->xmax;  
    iv_new->fmax = iv_old->fmax;
    iv_old->xmax  = iv_new->xmin = x; 
    iv_old->fmax = iv_new->fmin = fx; 
    break;
  default: 
    /* this should not happen */
    _unur_warning(gen->genid,UNUR_ERR_SAMPLE,"Invalid slope. Cannot split interval.");
    return 0;
  }

  /* compute the areas in both intervals */
  iv_old->Acum -= iv_old->Ahat;

  /** TODO: possible overflow/underflow ?? **/
  iv_new->Ahat     = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmax;
  iv_new->Asqueeze = iv_new->slope * (iv_new->xmax - iv_new->xmin) * iv_new->fmin;
  iv_old->Ahat     = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmax;
  iv_old->Asqueeze = iv_old->slope * (iv_old->xmax - iv_old->xmin) * iv_old->fmin;

  /* update cumulated areas */
  iv_old->Acum += iv_old->Ahat;
  iv_new->Acum = iv_old->Acum + iv_new->Ahat;


  /* insert iv_new into linked list of intervals.
     iv_old is stored on the left hand side of iv_new. */
  iv_new->next = iv_old->next;
  iv_old->next = iv_new;

  return iv_new;

} /* end of _unur_tabl_split_interval() */

/*****************************************************************************/

static int
_unur_tabl_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 (--> successful)                                                 */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0.                                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_tabl_interval *iv;
  double Acum, Asqueezecum, Astep;
  int j;

  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,0);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN.guide) {
    int max_guide_size = (GEN.guide_factor > 0.) ? (GEN.max_ivs * GEN.guide_factor) : 1;
    GEN.guide = _unur_malloc( max_guide_size * sizeof(struct unur_tabl_interval*) );
  }

  /* first we need the cumulated areas of rectangles */
  Acum = 0.;            /* area below hat */
  Asqueezecum = 0.;     /* area below squeeze */
  for (iv = GEN.iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    Acum += iv->Ahat;
    Asqueezecum += iv->Asqueeze;
    iv->Acum = Acum;
  }
    
  /* total area below hat */
  GEN.Atotal = Acum;
  GEN.Asqueeze = Asqueezecum;

  /* actual size of guide table */
  GEN.guide_size = GEN.n_ivs;

  /* make table (use variant 2; see dis.c) */
  Astep = GEN.Atotal / GEN.guide_size;
  Acum=0.;
  for( j=0, iv=GEN.iv; j < GEN.guide_size; j++ ) {
    COOKIE_CHECK(iv,CK_TABL_IV,0);
    while( iv->Acum < Acum )
      if( iv->next != NULL )    /* skip to next segment if it exists */
        iv = iv->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_INIT,"roundoff error while making guide table!");
	break;
      }
    GEN.guide[j] = iv;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN.guide_size ;j++ )
    GEN.guide[j] = iv;

  return 1;
} /* end of _unur_tabl_make_guide_table() */

/*****************************************************************************/

static struct unur_tabl_interval *
_unur_tabl_iv_stack_pop( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* pop free interval from stack; allocate memory block if necessary.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to interval                                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,NULL);

  /* look for an unused segment */
  if( ! GEN.iv_free ) {
    /* no more unused segments. make some. */
    GEN.iv_stack = _unur_malloc( UNUR_MALLOC_SIZE * sizeof(struct unur_tabl_interval) );

    /* reset counter */
    GEN.iv_free = UNUR_MALLOC_SIZE;

    /* set cookies */
    COOKIE_SET_ARRAY( GEN.iv_stack, CK_TABL_IV, UNUR_MALLOC_SIZE);

    /* add to list of allocated blocks */
    _unur_add_mblocks( &(GEN.mblocks), GEN.iv_stack ); 
  }

  /* update ....                                   */
  --(GEN.iv_free);   /* pointer to free segments  */
  ++(GEN.n_ivs);     /* counter for used segments */

  /* return pointer to segment */
  return (GEN.iv_stack + GEN.iv_free);

} /* end of _unur_tabl_iv_stack_pop() */

/*---------------------------------------------------------------------------*/

#if 0 /* we do not need this subroutine yet */
static void
_unur_tabl_iv_stack_push( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* push the last popped interval back onto the stack.                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  /* update counters and pointers */
  --(GEN.n_ivs);
  ++(GEN.iv_free);
} /* end of _unur_tabl_iv_stack_push() */
#endif

/*-----------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

#if UNUR_DEBUG & UNUR_DB_INFO

/*---------------------------------------------------------------------------*/

#define empty_line() fprintf(log,"%s:\n",gen->genid);

/*---------------------------------------------------------------------------*/

static void
_unur_tabl_debug_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into logfile                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  int i;

  log = unur_get_log();

  empty_line();
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  empty_line();

  fprintf(log,"%s: sampling routine = unur_tabl_sample",gen->genid);
  if (par->method & UNUR_MASK_SCHECK)
    fprintf(log,"_check()\n");
  else
    fprintf(log,"()\n");
  empty_line();

  fprintf(log,"%s: p.d.f with %d arguments\n",gen->genid,PAR.n_pdf_param);
  if (PAR.n_pdf_param)
    for( i=0; i<PAR.n_pdf_param; i++ )
      fprintf(log,"%s:\tparam[%d] = %g\n",gen->genid,i,PAR.pdf_param[i]);
  empty_line();

  if (par->method & UNUR_MASK_MODE )
    fprintf(log,"%s: mode = %g\n",gen->genid,PAR.mode);
  else
    fprintf(log,"%s: mode unknown\n",gen->genid);
  fprintf(log,"%s: domain = (%g, %g)",gen->genid,GEN.bleft,GEN.bright);
  if (!(par->set & UNUR_SET_DOMAIN))
      fprintf(log,"   (computed from given slopes)");
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: area below p.d.f. ",gen->genid);
  if (par->set & UNUR_SET_AREA)
    fprintf(log,"= %g\n",PAR.pdf_area);
  else
    fprintf(log,"not given. assume 1.\n");
  fprintf(log,"%s: area fraction for equal area rule = %g ",gen->genid,PAR.area_fract);
  _unur_print_if_default(par,UNUR_SET_TABL_C);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: maximum number of intervals        = %d",gen->genid,PAR.max_ivs);
  _unur_print_if_default(par,UNUR_SET_MAX_IVS);
  fprintf(log,"\n%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,PAR.max_ratio*100.);
  _unur_print_if_default(par,UNUR_SET_MAX_RATIO);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*PAR.guide_factor);
  _unur_print_if_default(par,UNUR_SET_FACTOR);
  fprintf(log,"\n");
  empty_line();

  fprintf(log,"%s: split intervals at ",gen->genid);
  switch( GEN.variant & TABL_MASK_SPLIT ) {
  case TABL_SPLIT_MEAN:
    fprintf(log,"mean point");
    break;
  case TABL_SPLIT_ARC:
    fprintf(log,"\"arcmean\" point");
    break;
  case TABL_SPLIT_POINT:
  default: 
    fprintf(log,"sample point");
    break;
  }
  fprintf(log," when using adaptive sampling.\n");
  empty_line();

  if (par->set & UNUR_SET_SLOPES) {
    fprintf(log,"%s: slopes = %d\n",gen->genid,PAR.n_slopes);
    for (i=0; i<PAR.n_slopes; i++) {
      if ( PAR.slopes[2*i] > PAR.slopes[2*i+1] )
	fprintf(log,"%s:   (+)  ",gen->genid);
      else
	fprintf(log,"%s:   (-)  ",gen->genid);
      fprintf(log,"< %#g, %#g >\n", PAR.slopes[2*i], PAR.slopes[2*i+1] );
    }
  }
  else
    fprintf(log,"%s: no slopes given. compute from domain and mode.\n",gen->genid);

  if (PAR.variant & TABL_STP_SPLIT_A)
    fprintf(log,"%s: split slopes by equal area rule (SPLIT A).\n",gen->genid);
  if (PAR.variant & TABL_STP_SPLIT_B)
    fprintf(log,"%s: split slopes by main subdivision rule (SPLIT B).\n",gen->genid);
  empty_line();

  fprintf(log,"%s: number of starting intervals (at least) = %d",gen->genid,PAR.n_starting_cpoints);
  _unur_print_if_default(par,UNUR_SET_N_STP);
  fprintf(log,"\n");
  empty_line();

} /* end of _unur_tabl_debug_init() */

/*****************************************************************************/

static void
_unur_tabl_debug_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into logfile            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  log = unur_get_log();

  empty_line();
  fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
  empty_line();
  _unur_tabl_debug_intervals(gen,TRUE);
  empty_line();

  fflush(log);

} /* end of _unur_tabl_debug_free() */

/*****************************************************************************/

static void
_unur_tabl_debug_intervals( struct unur_gen *gen, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into logfile                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *log;
  struct unur_tabl_interval *iv;
  double sAsqueeze;
  int i;

  /* check arguments */
  CHECK_NULL(gen,/*void*/);
  COOKIE_CHECK(gen,CK_TABL_GEN,/*void*/);

  log = unur_get_log();

  fprintf(log,"%s: intervals = %d\n",gen->genid,GEN.n_ivs);
  if (gen->debug & TABL_DB_IV) {
    fprintf(log,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    for (iv = GEN.iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,/*void*/);
      fprintf(log,"%s:[%3d]: (",gen->genid,i);
      switch (iv->slope) {
      case 1:  fprintf(log,"+"); break;
      case 0:  fprintf(log,"0"); break;
      case -1: fprintf(log,"-"); break;
      default:
	_unur_warning(gen->genid,UNUR_ERR_GENERIC,"invalid value for iv->slope.");
      }
      fprintf(log,")   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }

  if (!print_areas) return;

  if (GEN.Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    GEN.Atotal = -1.;   /* to avoid floating point exceptions */
  }

  /* print and sum areas below squeeze and hat */
  if (gen->debug & TABL_DB_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN.iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,/*void*/); 
      sAsqueeze += iv->Asqueeze;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / GEN.Atotal,
	      iv->Ahat, iv->Ahat * 100. / GEN.Atotal, 
	      iv->Acum, iv->Acum * 100. / GEN.Atotal);
    }
    fprintf(log,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / GEN.Atotal, GEN.Atotal);
    empty_line();
  }
    
  /* summary of areas */
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Asqueeze, GEN.Asqueeze * 100./GEN.Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN.Atotal - GEN.Asqueeze, (GEN.Atotal - GEN.Asqueeze) * 100./GEN.Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, GEN.Atotal);

  empty_line();

} /* end of _unur_tabl_debug_intervals */

/*****************************************************************************/

#endif
