/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      set.h                                                        *
 *                                                                           *
 *   set parameters for (new) generators                                     *
 *                                                                           *
 *   PARAMETER: struct unur_par *                                            *
 *                                                                           *
 *   return:                                                                 *
 *     1 ... on success                                                      *
 *     0 ... on error                                                        *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************
 *                                                                           *
 *   author: Josef.Leydold @ statistik.wu-wien.ac.at                         *
 *                                                                           *
 *   last modification: Mon Oct 18 16:59:52 CEST 1999                        *
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

#define _unur_warning_notrequired(par,param) \
   _unur_db_warning((par)->genid,UNUR_ERR_SET_NOTREQU,__FILE__,__LINE__,\
                    " - %s",param)

#define _unur_warning_invalid(par,param) \
   _unur_db_warning((par)->genid,UNUR_ERR_SET_INVALID,__FILE__,__LINE__,\
                    " - %s",param)

#define _unur_warning_set(par,param) \
   _unur_db_warning((par)->genid,UNUR_ERR_SET,__FILE__,__LINE__,\
                    "%s",param)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for the distribution and its p.d.f.                         **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_domain( struct unur_par *par, double left, double right )
/*---------------------------------------------------------------------------*/
/* set left and right boundary of domain of p.d.f.                           */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   left  ... left boundary point                                           */
/*   right ... right boundary point                                          */
/*                                                                           */
/* comment:                                                                  */
/*   the new boundary points may be +/- INFINITY                             */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning_invalid(par,"domain, left >= right");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.bleft = left;
    par->data.arou.bright = right;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    if (left <= -INFINITY || right >= INFINITY) {
      _unur_warning_invalid(par,"domain, +/- INFINITY not allowed");
      return 0;
    }
    par->data.tabl.bleft = left;
    par->data.tabl.bright = right;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.bleft = left;
    par->data.tdr.bright = right;
    break;

  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    par->data.utdr.il = left;
    par->data.utdr.ir = right;
    break;

  default:
    _unur_warning_notrequired(par,"domain");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_DOMAIN;

  /* o.k. */
  return 1;

} /* end of unur_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_set_domain_vec( struct unur_par *par, double **domain )
/*---------------------------------------------------------------------------*/
/* set coordinates for domain boundary                                       */
/*                                                                           */
/* parameters:                                                               */
/*   par    ... pointer to parameter for building generator object           */
/*   domain ... coordinates of domain boundary                               */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_RECT:
    COOKIE_CHECK(par,CK_RECT_PAR,0);
    par->data.rect.domain = domain;
    break;

  default:
    _unur_warning_notrequired(par,"domain");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_DOMAIN;

  /* o.k. */
  return 1;

} /* end of unur_set_domain_vec() */

/*---------------------------------------------------------------------------*/

int
unur_set_pdf_param(  struct unur_par *par, double *pdf_params, int n_params )
/*---------------------------------------------------------------------------*/
/* set array of parameters for p.d.f.                                        */
/*                                                                           */
/* parameters:                                                               */
/*   par        ... pointer to parameter for building generator object       */
/*   pdf_params ... list of arguments                                        */
/*   n_params   ... number of arguments                                      */
/*                                                                           */
/* comment:                                                                  */
/*   must be called unless p.d.f. has does not require parameters            */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (pdf_params == NULL) {
    _unur_warning_invalid(par,"p.d.f. params is NULL pointer");
    return 0;
  }
  if (n_params < 0) {
    _unur_warning_invalid(par,"number of p.d.f. params is negative");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.pdf_param = pdf_params;
    par->data.arou.n_pdf_param = n_params;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.pdf_param = pdf_params;
    par->data.tabl.n_pdf_param = n_params;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.pdf_param = pdf_params;
    par->data.tdr.n_pdf_param = n_params;
    break;

  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    par->data.utdr.pdf_param = pdf_params;
    par->data.utdr.n_pdf_param = n_params;
    break;

  default:
    _unur_warning_notrequired(par,"p.d.f. parameters");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_PDFPARAM;

  /* o.k. */
  return 1;

} /* end of unur_set_pdf_param() */

/*---------------------------------------------------------------------------*/

int
unur_set_mode( struct unur_par *par, double mode )
/*---------------------------------------------------------------------------*/
/* set mode of p.d.f.                                                        */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   mode  ... mode of p.d.f.                                                */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.mode = mode;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.mode = mode;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.mode = mode;
    break;

  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    par->data.utdr.mode = mode;
    break;

  default:
    _unur_warning_notrequired(par,"mode of p.d.f.");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_MODE;

  /* indicate that mode is known (we use a bit in method) */
  par->method |= UNUR_MASK_MODE;

  /* o.k. */
  return 1;

} /* end of unur_set_mode() */

/*---------------------------------------------------------------------------*/

int
unur_set_usemode( struct unur_par *par, int usemode )
/*---------------------------------------------------------------------------*/
/* set flag for using mode of p.d.f.                                         */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   mode  ... mode of p.d.f.                                                */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    break;

  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    break;

  default:
    _unur_warning_notrequired(par,"mode of p.d.f.");
    return 0;
  }
  
  /* indicate that mode should be used (we use a bit in method) */
  par->method = (usemode) ? (par->method | UNUR_MASK_MODE) : (par->method & (~UNUR_MASK_MODE));

  /* o.k. */
  return 1;

} /* end of unur_set_usemode() */

/*---------------------------------------------------------------------------*/

int
unur_set_pdf_area( struct unur_par *par, double area )
/*---------------------------------------------------------------------------*/
/* set the (approximate) area below p.d.f.                                   */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   area  ... area                                                          */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (area <= 0.) {
    _unur_warning_invalid(par,"pdf area <= 0");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.pdf_area = area;
    break;

  case UNUR_METH_UTDR:
    COOKIE_CHECK(par,CK_UTDR_PAR,0);
    par->data.utdr.pdf_area = area;
    break;

  default:
    _unur_warning_notrequired(par,"area below p.d.f.");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_AREA;

  /* o.k. */
  return 1;

} /* end of unur_set_pdf_area() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of discrete distributions                    **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_factor( struct unur_par *par, double factor )
/*---------------------------------------------------------------------------*/
/* set factor for relative size of (search|guide|alias) table                */
/*                                                                           */
/* parameters:                                                               */
/*   par    ... pointer to parameter for building generator object           */
/*   factor ... relative size for table                                      */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning_invalid(par,"relative table size < 0");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_DAU:
    COOKIE_CHECK(par,CK_DAU_PAR,0);
    if (factor < 1.) {
      _unur_warning_invalid(par,"relative urn size < 1");
      return 0;
    }
    par->data.dau.urn_factor = factor;
    break;

  case UNUR_METH_DIS:
    COOKIE_CHECK(par,CK_DIS_PAR,0);
    par->data.dis.guide_factor = factor;
    break;

  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.guide_factor = factor;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.guide_factor = factor;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.guide_factor = factor;
    break;

  default:
    _unur_warning_notrequired(par,"size of search table");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_FACTOR;

  return 1;

} /* end of unur_set_factor() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for generators of continuous distributions                  **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_cpoints( struct unur_par *par, int n_stp, double *stp )
/*---------------------------------------------------------------------------*/
/* set construction points for hat and/or its number for initialization      */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   n_stp ... number of starting points                                     */
/*   stp   ... pointer to array of starting points                           */
/*             (NULL for changing only the number of default points)         */
/*---------------------------------------------------------------------------*/
{
  int i; 

  /* check arguments */
  CHECK_NULL(par,0);

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning_invalid(par,"number of starting points < 0");
    return 0;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning_invalid(par,"starting points not strictly monotonically increasing");
	return 0;
      }


  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.starting_cpoints = stp;
    par->data.arou.n_starting_cpoints = n_stp;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    _unur_warning_notrequired(par,"starting points for hat");
    par->data.tabl.n_starting_cpoints = n_stp;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.starting_cpoints = stp;
    par->data.tdr.n_starting_cpoints = n_stp;
    break;

  default:
    _unur_warning_notrequired(par,"construction points for hat");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_N_STP | ((stp) ? UNUR_SET_STP : 0);

  /* o.k. */
  return 1;

} /* end of unur_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_set_max_ratio( struct unur_par *par, double max_ratio )
/*---------------------------------------------------------------------------*/
/* set bound for ratio A(squeeze) / A(hat)                                   */
/*                                                                           */
/* parameters:                                                               */
/*   par       ... pointer to parameter for building generator object        */
/*   max_ratio ... upper bound for ratio to add a new construction point     */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning_invalid(par,"ratio Atotal / Asqueeze not in [0,1]");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.max_ratio = max_ratio;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.max_ratio = max_ratio;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.max_ratio = max_ratio;
    break;

  default:
    _unur_warning_notrequired(par,"bound for ratio A(squeeze) / A(hat)");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_MAX_RATIO;

  /* o.k. */
  return 1;

} /* end of unur_set_max_ratio() */

/*---------------------------------------------------------------------------*/

int
unur_set_max_intervals( struct unur_par *par, int max_ivs )
/*---------------------------------------------------------------------------*/
/* set maximum number of intervals or segments                               */
/*                                                                           */
/* parameters:                                                               */
/*   par     ... pointer to parameter for building generator object          */
/*   max_ivs ... maximum number of intervals                                 */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning_invalid(par,"maximum number of intervals/segments < 1");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    par->data.arou.max_segs = max_ivs;
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.max_ivs = max_ivs;
    break;

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    par->data.tdr.max_ivs = max_ivs;
    break;

  default:
    _unur_warning_notrequired(par,"maximal number of intersvals/segments");
    return 0;
  }
  
  /* changelog */
  par->set |= UNUR_SET_MAX_IVS;

  /* o.k. */
  return 1;

} /* end of unur_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int
unur_set_tdr_c( struct unur_par *par, double c )
/*---------------------------------------------------------------------------*/
/* set parameter c for transformation T_c                                    */           
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*   c   ... parameter c                                                     */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    if (c > 0.) {
      _unur_warning_invalid(par,"c > 0");
      return 0;
    }
    if (c <= -1.) {
      _unur_warning_invalid(par,"c <= -1 only if domain is bounded. Use `TABL' method then.");
      return 0;
    }
    if (c != 0 && c > -0.5) {
      _unur_warning_invalid(par,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
      c = -0.5;
    }
    
    par->data.tdr.c_T = c;
    break;

  default:
    _unur_warning_notrequired(par,"TDR c");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_TDR_C;

  return 1;

} /* end of unur_set_tdr_c() */

/*---------------------------------------------------------------------------*/

int
unur_set_tabl_c( struct unur_par *par, double c )
/*---------------------------------------------------------------------------*/
/* set parameter for equal area rule (each bar has size c * area(pdf))       */           
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*   c   ... area factor                                                     */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    if (c < 0.) {
      _unur_warning_invalid(par,"area factor < 0");
      return 0;
    }
    par->data.tabl.area_fract = c;
    break;

  default:
    _unur_warning_notrequired(par,"TABL c");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_TABL_C;

  return 1;

} /* end of unur_set_tabl_c() */

/*---------------------------------------------------------------------------*/

int
unur_set_slopes( struct unur_par *par, double *slopes, int n_slopes )
/*---------------------------------------------------------------------------*/
/* set slopes of p.d.f.                                                      */
/*                                                                           */
/* parameters:                                                               */
/*   par      ... pointer to parameter for building generator object         */
/*   slopes   ... pointer to list of slopes                                  */
/*   n_slopes ... number of slopes                                           */
/*                                                                           */
/* comment:                                                                  */
/*   a slope <a,b> is an interval [a,b] or [b,a] such that pdf(a) >= pdf(b)  */
/*                                                                           */
/*   slopes must be decreasing, non-overlapping and sorted                   */
/*---------------------------------------------------------------------------*/
{
  int i;
  double al, bl;

  /* check arguments */
  CHECK_NULL(par,0);

  /* check new parameter for generator */
  if( n_slopes <= 0 ) {
    _unur_warning_invalid(par,"number of slopes <= 0");
    return 0;
  }

  /* check slopes */
  al = slopes[0];
  bl = slopes[1];
  for( i=1; i<n_slopes; i++ ) {
    /* we do not check here if f(a) >= f(b), since we make no calculations heres */
    if( al > slopes[2*i] || bl > slopes[2*i+1] ) {
      _unur_warning_invalid(par,"slopes overlapping or not in ascending order");
      return 0;
    }
    al = slopes[2*i];
    bl = slopes[2*i+1];
  }

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    par->data.tabl.slopes = slopes;
    par->data.tabl.n_slopes = n_slopes;
    break;

  default:
    _unur_warning_notrequired(par,"slopes");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_SLOPES;

  return 1;

} /* end of unur_set_slopes() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Parameters for all generators                                          **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_set_urng( struct unur_par *par, UNUR_URNG_TYPE urng )
/*---------------------------------------------------------------------------*/
/* set uniform random number generator                                       */
/*                                                                           */
/* parameters:                                                               */
/*   par     ... pointer to parameter for building generator object          */
/*   urng    ... pointer to uniform random number generator                  */
/*                                                                           */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);
  CHECK_NULL(urng,0);

  par->urng = urng;

  return 1;
} /* end of unur_set_urng() */

/*---------------------------------------------------------------------------*/

int
unur_set_variant( struct unur_par *par, unsigned long variant )
/*---------------------------------------------------------------------------*/
/* set variant of method                                                     */
/*                                                                           */
/* parameters:                                                               */
/*   par     ... pointer to parameter for building generator object          */
/*   variant ... indicator for variant                                       */
/*                                                                           */
/* comment:                                                                  */
/*   see method´s file for descreption of variants.                          */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* test size of variant */
  if (variant & ~UNUR_MASK_VARIANT) {
    _unur_warning_invalid(par,"variant, too many bits");
    return 0;
  }

  switch (par->method & UNUR_MASK_METHOD) {

  case UNUR_METH_DIS:
    COOKIE_CHECK(par,CK_DIS_PAR,0);
    /* check new parameter for generator */
    if (variant < 0 || variant > 2) {
      _unur_warning_invalid(par,"variant");
      return 0;
    }
    break;

  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    break;

  default:
    _unur_warning_set(par," - no variants for method");
    return 0;
  }

  /* changelog */
  par->set |= UNUR_SET_VARIANT;

  par->method = (par->method & (~UNUR_MASK_VARIANT)) | variant;

  return 1;

} /* end of unur_set_variant() */

/*---------------------------------------------------------------------------*/

int
unur_set_check(  struct unur_par *par, int check )
/*---------------------------------------------------------------------------*/
/* turn testing of sampling on/off                                           */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   check ... 0  = no check                                                 */
/*             !0 = check                                                    */
/*                                                                           */
/* comment:                                                                  */
/*   no checking is default                                                  */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  switch (par->method & UNUR_MASK_METHOD) {
  case UNUR_METH_AROU:
    COOKIE_CHECK(par,CK_AROU_PAR,0);
    break;
  case UNUR_METH_TABL:
    COOKIE_CHECK(par,CK_TABL_PAR,0);
    break;
  case UNUR_METH_TDR:
    COOKIE_CHECK(par,CK_TDR_PAR,0);
    break;

  default:
    _unur_warning_set(par," - no checking mode for method");
    return 0;
  }

  /* we use a bit in method */
  par->method = (check) ? (par->method | UNUR_MASK_SCHECK) : (par->method & (~UNUR_MASK_SCHECK));
  
  /* o.k. */
  return 1;

} /* end of unur_set_check() */

/*---------------------------------------------------------------------------*/

int
unur_set_copyall(  struct unur_par *par, int copy )
/*---------------------------------------------------------------------------*/
/* turn copaing of all inputs into generator object on/off                   */
/* (for some methods, some input parameters are not necessary to run         */
/* the generator but is a usefull debugging information                      */
/*                                                                           */
/* parameters:                                                               */
/*   par   ... pointer to parameter for building generator object            */
/*   copy  ... 0  = do not copy unnecessary parameters                       */
/*             !0 = copy all                                                 */
/*                                                                           */
/* comment:                                                                  */
/*   no copying is default                                                   */
/*---------------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(par,0);

  /* we use a bit in method */
  par->method = (copy) ? (par->method | UNUR_MASK_COPYALL) : (par->method & (~UNUR_MASK_COPYALL));
  
  /* o.k. */
  return 1;

} /* end of unur_set_copyall() */

/*---------------------------------------------------------------------------*/

int
unur_set_debug( struct unur_par *par, unsigned long db )
/*---------------------------------------------------------------------------*/
/* set debugging flag for generator                                          */
/*                                                                           */
/* parameters:                                                               */
/*   par ... pointer to parameter for building generator object              */
/*   db  ... debugging flag                                                  */
/*---------------------------------------------------------------------------*/
{
#if UNUR_DEBUG & UNUR_DB_INFO
  par->debug = db;
  return 1;
#else
  _unur_warning_set(par,"Debugging not available. Recompile library.");
  return 0;
#endif
} /* end of unur_set_debug() */
  
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Auxilliary routines                                                    **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

