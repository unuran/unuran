/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distr_corder.c                                               *
 *                                                                           *
 *   manipulate univariate continuous order statistics                       *
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

#include <source_distributions.h>

/*---------------------------------------------------------------------------*/

static const char distr_name[] = "order statistics";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont    /* underlying (base) distribution          */
#define OS    os->data.cont       /* order statistics                        */

/*---------------------------------------------------------------------------*/

/* function prototypes                                                       */
static double _unur_pdf_corder( double x, UNUR_DISTR *os );
static double _unur_dpdf_corder( double x, UNUR_DISTR *os );
static double _unur_cdf_corder( double x, UNUR_DISTR *os );

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous order statistics                                  **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_corder_new( struct unur_distr *distr, int n, int k )
     /*----------------------------------------------------------------------*/
     /* Create an object for order statistics of for a                       */
     /* sample size n and rank k.                                            */
     /* `distr' must be a pointer to a univariate continuous distribution.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to  univariate continuous distribution.          */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_distr *os;

  /* check arguments */
  _unur_check_NULL( distr_name,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  /* check parameters n and k */
  if (n < 2 || k < 1 || k > n) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,"n < 2 or k < 1 or k > n");
    return NULL;
  }

  /* allocate structure */
  os = _unur_malloc( sizeof(struct unur_distr) );
  if (!os) return NULL;

  /* set magic cookie */
  COOKIE_SET(os,CK_DISTR_CONT);

  /* set type of distribution */
  os->type = UNUR_DISTR_CONT;

  /* set id to generic distribution */
  os->id = UNUR_DISTR_CORDER;

  /* name of distribution */
  os->name = distr_name;

  /* this is a derived distribution */
  /* allocate memory ... */
  os->base = _unur_malloc( sizeof(struct unur_distr) );
  if (!os->base) return NULL;
  /* ... and copy distribution object */
  memcpy( os->base, distr, sizeof(struct unur_distr) );

  /* set parameters for order statistics */
  OS.n_params = 2;                 /* two parameters: n and k                */
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;
  /* there is no need (?) to set the other parameters to 0 */

  /* copy data */
  OS.area = DISTR.area;            /* area below p.d.f. (same as for distr)  */
  OS.domain[0] = DISTR.domain[0];  /* left boundary of domain                */
  OS.domain[1] = DISTR.domain[1];  /* right boundary of domain               */

  OS.pdf  = (DISTR.pdf)  ? _unur_pdf_corder  : NULL;  /* pointer to p.d.f.   */
  OS.dpdf = (DISTR.dpdf) ? _unur_dpdf_corder : NULL;  /* derivative of p.d.f.*/
  OS.cdf  = (DISTR.cdf)  ? _unur_cdf_corder  : NULL;  /* pointer to c.d.f.   */

  /* set defaults                                                            */
  OS.mode         = INFINITY;      /* location of mode (default: not known)  */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  /* there is no function for computing the mode of the order statistics     */
  DISTR.upd_mode  = NULL;

  /* there is no necessity for a function that computes the area below pdf   */
  DISTR.upd_area  = NULL;

  /* parameters set */
  os->set = distr->set & ~UNUR_DISTR_SET_MODE; /* mode not derived from distr */

  /* return pointer to object */
  return os;

} /* end of unur_distr_corder_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_corder_get_distribution( struct unur_distr *os )
     /*----------------------------------------------------------------------*/
     /* get pointer to distribution object for underlying distribution       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os ... pointer to order statistics                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to underlying distribution                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, NULL );
  _unur_check_distr_object( os, CONT, NULL );

  /* check distribution */
  if (os->type != UNUR_DISTR_CORDER) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"");
    return NULL; 
  }

  return os->base;
} /* end of unur_distr_corder_get_distribution() */

/*---------------------------------------------------------------------------*/

int
unur_distr_corder_set_rank( struct unur_distr *os, int n, int k )
     /*----------------------------------------------------------------------*/
     /* change sample size n and rank k of order statistics.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, 0 );
  _unur_check_distr_object( os, CONT, 0 );

  /* check distribution */
  if (os->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return 0; }
  COOKIE_CHECK(os,CK_DISTR_CONT,0);

  /* check parameters n and k */
  if (n < 2 || k < 1 || k > n) {
    _unur_warning(distr_name,UNUR_ERR_DISTR_INVALID,"n < 2 or k < 1 or k > n");
    return 0;
  }

  /* changelog */
  os->set &= ~UNUR_DISTR_SET_MODE; /* mode unknown */

  /* copy parameters */
  OS.params[0] = (double) n;
  OS.params[1] = (double) k;

  /* o.k. */
  return 1;
} /* end of unur_distr_corder_set_rank() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_corder_get_rank( struct unur_distr *os, int *n, int *k )
     /*----------------------------------------------------------------------*/
     /* get sample size n and rank k of order statistics.                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   n     ... sample size                                              */
     /*   k     ... rank                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( distr_name, os, 0 );
  _unur_check_distr_object( os, CONT, 0 );

  /* check distribution */
  if (os->type != UNUR_DISTR_CONT) {
    _unur_error(distr_name,UNUR_ERR_DISTR_INVALID,""); return 0; }
  COOKIE_CHECK(os,CK_DISTR_CONT,0);

  /* copy parameters */
  *n = (int)(OS.params[0] + 0.5);
  *k = (int)(OS.params[1] + 0.5);

  /* o.k. */
  return 1;
} /* end of unur_distr_corder_get_rank() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** p.d.f., its derivative and c.d.f. of order statistics                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double
_unur_pdf_corder( double x, struct unur_distr *os )
     /* 
	pdf(x) = b(F(x)) * f(x)

	b(.) ... pdf of beta(k,n-k+1) distribution
	f(.) ... pdf of underlying distribution
	F(.) ... cdf of underlying distribution
     */
{ 
  double Fx;    /* cdf of underlying distribution at x */
  double fx;    /* pdf of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */

  double x_LOGNORMCONSTANT;

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  CHECK_NULL( os->base->data.cont.pdf, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  x_LOGNORMCONSTANT = _unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q);

  /* pdf(x) = b(F(x)) * f(x) */
  return exp(log(fx) + (p-1.)*log(Fx) + (q-1.)*log(1.-Fx) - x_LOGNORMCONSTANT);

} /* end of _unur_pdf_corder() */

/*---------------------------------------------------------------------------*/

double
_unur_dpdf_corder( double x, struct unur_distr *os )
     /* 
	pdf'(x) = b'(F(x)) * f(x)^2 + b(F(x)) * f'(x)

	b(.) ... pdf of beta(k,n-k+1) distribution
	f(.) ... pdf of underlying distribution
	F(.) ... cdf of underlying distribution
     */
{
  double Fx;    /* cdf of underlying distribution at x */
  double fx;    /* pdf of underlying distribution at x */
  double dfx;   /* derivative of pdf of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */
  double dpdf;  /* derivative of pdf of order statistics */
  double lFx, lFy;

  double x_LOGNORMCONSTANT;

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  CHECK_NULL( os->base->data.cont.dpdf, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);
  fx = (*(os->base->data.cont.pdf)) (x, os->base);
  dfx = (*(os->base->data.cont.dpdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  x_LOGNORMCONSTANT = _unur_gammaln(p) + _unur_gammaln(q) - _unur_gammaln(p+q);

  lFx = log(Fx);
  lFy = log(1.-Fx);

  /* pdf'(x) = b'(F(x)) * f(x)^2 + b(F(x)) * f'(x) */
  dpdf = ( exp(2.*log(fx) + (p-2.)*lFx + (q-2.)*lFy - x_LOGNORMCONSTANT)
	   * ( (p-1.)*(1.-Fx) - (q-1.)*Fx ));
  dpdf += exp((p-1.)*lFx + (q-1.)*lFy - x_LOGNORMCONSTANT) * dfx;

  return dpdf;
} /* end of _unur_dpdf_corder() */

/*---------------------------------------------------------------------------*/

double
_unur_cdf_corder( double x, struct unur_distr *os ) 
     /* 
	cdf(x) = B(F(x))

	B(.) ... cdf of beta(k,n-k+1) distribution
	F(.) ... cdf of underlying distribution
     */
{
  double Fx;    /* cdf of underlying distribution at x */
  double p,q;   /* shape parameters for beta distribution */

  /* check arguments */
  _unur_check_NULL( NULL, os, INFINITY );
  CHECK_NULL( os->base, INFINITY );
  CHECK_NULL( os->base->data.cont.cdf, INFINITY );
  _unur_check_distr_object( os, CONT, INFINITY );

  Fx = (*(os->base->data.cont.cdf)) (x, os->base);

  p = OS.params[1];                       /* k     */
  q = OS.params[0] - OS.params[1] + 1.;   /* n-k+1 */

  /* cdf(x) = B(F(x)) */
  return _unur_incbeta(Fx,p,q);

} /* end of _unur_cdf_corder() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** debug                                                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

void
_unur_distr_corder_debug( struct unur_distr *os, char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   os    ... pointer to order statistics                              */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
  FILE *log;

  /* check arguments */
  CHECK_NULL(os,/*void*/);
  COOKIE_CHECK(os,CK_DISTR_CONT,/*void*/);

  log = unur_get_stream();

  /* print data about order statistics */
  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = order statistics of continuous univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,os->name);
  fprintf(log,"%s:\tsample size\tn = %d\n",genid,(int)(OS.params[0]+0.5));
  fprintf(log,"%s:\trank\t\tk = %d\n",genid,(int)(OS.params[1]+0.5));
  fprintf(log,"%s:\n",genid);

  if (os->set & UNUR_DISTR_SET_MODE)
    fprintf(log,"%s:\tmode = %g\n",genid,OS.mode);
  else
    fprintf(log,"%s:\tmode unknown\n",genid);

  fprintf(log,"%s:\tdomain = (%g, %g)",genid,OS.domain[0],OS.domain[1]);
  _unur_print_if_default(os,UNUR_DISTR_SET_DOMAIN);

  fprintf(log,"\n%s:\tarea below p.d.f. = %g",genid,OS.area);
  _unur_print_if_default(os,UNUR_DISTR_SET_PDFAREA);
  fprintf(log,"\n%s:\n",genid);

  /* print data about underlying distribution */
  fprintf(log,"%s: Underlying distribution:\n",genid);
  _unur_distr_cont_debug(os->base, genid);

} /* end of _unur_distr_corder_debug() */

/*---------------------------------------------------------------------------*/

