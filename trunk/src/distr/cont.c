/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cont.c                                                       *
 *                                                                           *
 *   manipulate univariate continuous distribution objects                   *
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
#include <utils/fminmax_source.h>	
#include "distr_source.h"
#include "distr.h"
#include "cont.h"

/*---------------------------------------------------------------------------*/

static const char unknown_distr_name[] = "unknown";

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.cont

/* for derived distributions (e.g. order statistics):
   data of underlying distributions */
#define BASE  distr->base->data.cont

/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for PDF.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for dPDF.                                          */
/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for CDF.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for HR.                                            */
/*---------------------------------------------------------------------------*/

void _unur_distr_cont_free( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* destroy distribution object.                                              */
/*---------------------------------------------------------------------------*/

static int _unur_distr_cont_find_mode( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* find mode of unimodal univariate PDF numerically                          */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate continuous distributions                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_cont_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate continuous with given p.d.f.                        */
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
  int i;

  /* allocate structure */
  distr = _unur_malloc( sizeof(struct unur_distr) );
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_CONT);

  /* set type of distribution */
  distr->type = UNUR_DISTR_CONT;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* name of distribution */
  distr->name = unknown_distr_name;
  distr->name_str = NULL;

  /* this is not a derived distribution */
  distr->base = NULL;

  /* destructor */
  distr->destroy = _unur_distr_cont_free;

  /* clone */
  distr->clone = _unur_distr_cont_clone;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to PDF                         */
  DISTR.dpdf      = NULL;          /* pointer to derivative of PDF           */
  DISTR.cdf       = NULL;          /* pointer to CDF                         */
  DISTR.hr        = NULL;          /* pointer to HR                          */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.n_params  = 0;             /* number of parameters of the pdf        */
  /* initialize parameters of the p.d.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for p.d.f.
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.mode       = INFINITY;     /* location of mode (default: not known)  */
  DISTR.area       = 1.;           /* area below PDF (default: not known)    */

  DISTR.trunc[0] = DISTR.domain[0] = -INFINITY; /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INFINITY;  /* right boundary of domain  */

  DISTR.set_params = NULL;         /* funct for setting parameters and domain*/
  DISTR.upd_mode   = _unur_distr_cont_find_mode; /* funct for computing mode */
  DISTR.upd_area   = NULL;         /* funct for computing area               */

  DISTR.pdftree    = NULL;         /* pointer to function tree for PDF       */
  DISTR.dpdftree   = NULL;         /* pointer to function tree for dPDF      */
  DISTR.cdftree    = NULL;         /* pointer to function tree for CDF       */
  DISTR.hrtree     = NULL;         /* pointer to function tree for HR        */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

} /* end of unur_distr_cont_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_cont_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.cont

  struct unur_distr *clone;
  int len;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  /* allocate memory */
  clone = _unur_malloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy function trees into generator object (when there is one) */
  CLONE.pdftree  = (DISTR.pdftree)  ? _unur_fstr_dup_tree(DISTR.pdftree)  : NULL;
  CLONE.dpdftree = (DISTR.dpdftree) ? _unur_fstr_dup_tree(DISTR.dpdftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree)  ? _unur_fstr_dup_tree(DISTR.cdftree)  : NULL;
  CLONE.hrtree   = (DISTR.hrtree)   ? _unur_fstr_dup_tree(DISTR.hrtree)   : NULL;

  /* copy user name for distribution */
  if (distr->name_str) {
    len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_malloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  /* for a derived distribution we also have to copy the underlying */
  /* distribution object                                            */
  if (distr->base != NULL) {
    clone->base = _unur_distr_cont_clone( distr->base);
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_cont_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_cont_free( struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, RETURN_VOID );

  /* function trees */
  if (DISTR.pdftree)  _unur_fstr_free(DISTR.pdftree);
  if (DISTR.dpdftree) _unur_fstr_free(DISTR.dpdftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);
  if (DISTR.hrtree)   _unur_fstr_free(DISTR.hrtree);

  /* derived distribution */
  if (distr->base) _unur_distr_free(distr->base);

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);
  
  COOKIE_CLEAR(distr);
  free( distr );
} /* end of _unur_distr_cont_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdf( struct unur_distr *distr, UNUR_FUNCT_CONT *pdf )
     /*----------------------------------------------------------------------*/
     /* set PDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to PDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* we do not allow overwriting a pdf */
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.pdf = pdf;
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_dpdf( struct unur_distr *distr, UNUR_FUNCT_CONT *dpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of PDF of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   dpdf   ... pointer to derivative of PDF                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, dpdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a dpdf */
  if (DISTR.dpdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of dPDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.dpdf = dpdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_dpdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdf( struct unur_distr *distr, UNUR_FUNCT_CONT *cdf )
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
  _unur_check_NULL( distr->name, cdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a cdf */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_hr( struct unur_distr *distr, UNUR_FUNCT_CONT *hr )
     /*----------------------------------------------------------------------*/
     /* set HR of distribution                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   hr    ... pointer to HR                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, hr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting a cdf */
  if (DISTR.hr != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of HR not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.hr = hr;
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_hr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfstr( struct unur_distr *distr, const char *pdfstr )
     /*----------------------------------------------------------------------*/
     /* set PDF and its derivative of distribution via a string interface    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   pdfstr ... string that describes function term of PDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a PDF */
  if (DISTR.pdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse PDF string */
  if ( (DISTR.pdftree = _unur_fstr2tree(pdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pdf  = _unur_distr_cont_eval_pdf_tree;

  /* make derivative */
  if ( (DISTR.dpdftree = _unur_fstr_make_derivative(DISTR.pdftree)) == NULL )
    return UNUR_ERR_DISTR_DATA;
  DISTR.dpdf = _unur_distr_cont_eval_dpdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_pdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
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
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }

  /* set evaluation function */
  DISTR.cdf  = _unur_distr_cont_eval_cdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_cdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_hrstr( struct unur_distr *distr, const char *hrstr )
     /*----------------------------------------------------------------------*/
     /* set HR of distribution via a string interface                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   hrstr  ... string that describes function term of CDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, hrstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.hr != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.hrtree = _unur_fstr2tree(hrstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }

  /* set evaluation function */
  DISTR.hr  = _unur_distr_cont_eval_hr_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_hrstr() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_pdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for PDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for PDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.pdftree) ? _unur_fstr_eval_tree(DISTR.pdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_pdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_dpdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for derivative of PDF.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for derivative of PDF                           */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   derivative of PDF at x                                             */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.dpdftree) ? _unur_fstr_eval_tree(DISTR.dpdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_dpdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_cdf_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for CDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_cdf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_cont_eval_hr_tree( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for HR.                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at x                                                           */
     /*----------------------------------------------------------------------*/
{
  return ((DISTR.hrtree) ? _unur_fstr_eval_tree(DISTR.hrtree,x) : INFINITY);
} /* end of _unur_distr_cont_eval_hr_tree() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_pdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get PDF string that is given via the string interface                */
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
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.pdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.pdftree,"x","PDF",TRUE);
} /* end of unur_distr_cont_get_pdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_dpdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get string for derivative of PDF that is given via string interface  */
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
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.dpdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.dpdftree,"x","dPDF",TRUE);
} /* end of unur_distr_cont_get_dpdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_cdfstr( const struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} /* end of unur_distr_cont_get_cdfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_cont_get_hrstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get HR string that is given via the string interface                 */
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
  _unur_check_distr_object( distr, CONT, NULL );
  _unur_check_NULL( NULL, DISTR.hrtree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.hrtree,"x","HR",TRUE);
} /* end of unur_distr_cont_get_hrstr() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_pdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to PDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to PDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.pdf;
} /* end of unur_distr_cont_get_pdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_dpdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to derivative of PDF of distribution                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to derivative of PDF                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.dpdf;
} /* end of unur_distr_cont_get_dpdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_cdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to CDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to CDF                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.cdf;
} /* end of unur_distr_cont_get_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_CONT *
unur_distr_cont_get_hr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to HR of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to HR                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, CONT, NULL );

  return DISTR.hr;
} /* end of unur_distr_cont_get_hr() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_pdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PDF of distribution at x                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for pdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pdf(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.pdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_PDF(x,distr);
} /* end of unur_distr_cont_eval_pdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_dpdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate derivative of PDF of distribution at x                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for dpdf                                        */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   (pdf(x))'                                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.dpdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_dPDF(x,distr);
} /* end of unur_distr_cont_eval_dpdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_cdf( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of distribution at x                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   cdf(x)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.cdf == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_CDF(x,distr);
} /* end of unur_distr_cont_eval_cdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_cont_eval_hr( double x, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate HR of distribution at x                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   x     ... argument for cdf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   hr(x)                                                              */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  if (DISTR.hr == NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_cont_HR(x,distr);
} /* end of unur_distr_cont_eval_hr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfparams( struct unur_distr *distr, const double *params, int n_params )
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
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name,params,UNUR_ERR_NULL);

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* even if the set routine fails, the derived parameters are
     marked as unknown. but this is o.k. since in this case something
     has been wrong. */

  /* use special routine for setting parameters
     (if there is one) */

  if (distr->base && BASE.set_params) 
    return (BASE.set_params(distr->base,params,n_params));

  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));

  /* otherwise simply copy parameters */

  /* but first check number of new parameter for the distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  if (distr->base) {
    BASE.n_params = n_params;
    if (n_params) memcpy( BASE.params, params, n_params*sizeof(double) );
  }

  else {
    DISTR.n_params = n_params;
    if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_pdfparams( const struct unur_distr *distr, const double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pdf parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number of pdf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, CONT, 0 );

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       the parameters for the underlying distributions are returned */
    *params = (BASE.n_params) ? BASE.params : NULL;
    return BASE.n_params;
  }
  else {
    *params = (DISTR.n_params) ? DISTR.params : NULL;
    return DISTR.n_params;
  }

} /* end of unur_distr_cont_get_pdfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_domain( struct unur_distr *distr, double left, double right )
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
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (left >= right) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* we have to deal with the mode */
  if ( (distr->set & UNUR_DISTR_SET_MODE) &&
       (left  >= DISTR.domain[0]) &&
       (right <= DISTR.domain[1]) ) {
    /* mode has been set and new domain is subset of old domain */
    /* we assume the density is unimodal and thus monotone on
       either side of the mode!! */
    if ( DISTR.mode < left)       DISTR.mode = left;
    else if ( DISTR.mode > right) DISTR.mode = right;
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
  if (distr->set & UNUR_DISTR_SET_MODE) {
    distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		    UNUR_DISTR_SET_TRUNCATED | 
		    UNUR_DISTR_SET_MASK_DERIVED );
    distr->set |= UNUR_DISTR_SET_MODE;
  }
  else {
    distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		    UNUR_DISTR_SET_TRUNCATED | 
		    UNUR_DISTR_SET_MASK_DERIVED );
  }

  if (distr->base) {
    /* for derived distributions (e.g. order statistics)
       we also set the domain for the underlying distribution */
    BASE.trunc[0] = BASE.domain[0] = left;
    BASE.trunc[1] = BASE.domain[1] = right;
    distr->base->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
			  UNUR_DISTR_SET_TRUNCATED | 
			  UNUR_DISTR_SET_MASK_DERIVED );
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_domain( const struct unur_distr *distr, double *left, double *right )
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
     /*   if no boundaries have been set +/- INFINITY is returned.           */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = -INFINITY;
  *right = INFINITY;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_get_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_get_truncated( const struct unur_distr *distr, double *left, double *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the truncated domain of the        */
     /* distribution                                                         */
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
     /*   if no boundaries have been set +/- INFINITY is returned.           */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = -INFINITY;
  *right = INFINITY;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[0] : DISTR.domain[0];
  *right = (distr->set & UNUR_DISTR_SET_TRUNCATED) ? DISTR.trunc[1] : DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_get_truncated() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_mode( struct unur_distr *distr, double mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of p.d.f.                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_cont_set_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_mode( struct unur_distr *distr )
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
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

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

} /* end of unur_distr_cont_upd_mode() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_mode( struct unur_distr *distr )
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
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INFINITY;
    }
    else {
      /* compute mode */
      if (unur_distr_cont_upd_mode(distr)!=UNUR_SUCCESS) {
	/* finding mode not successfully */
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INFINITY;
      }
    }
  }

  return DISTR.mode;

} /* end of unur_distr_cont_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_cont_set_pdfarea( struct unur_distr *distr, double area )
     /*----------------------------------------------------------------------*/
     /* set area below p.d.f.                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   area  ... area below p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (area <= 0.) {
    _unur_error(NULL,UNUR_ERR_DISTR_SET,"pdf area <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  DISTR.area = area;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_cont_set_pdfarea() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_cont_upd_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute area below p.d.f. of distribution (if possible)        */
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
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );

  if (DISTR.upd_area == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute area */
  if (((DISTR.upd_area)(distr)!=UNUR_SUCCESS) || DISTR.area <= 0.) {
    /* computing of area failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"upd area <= 0");
    DISTR.area = 1.;   /* avoid possible floating point exceptions */
    distr->set &= ~UNUR_DISTR_SET_PDFAREA;
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  return UNUR_SUCCESS;
} /* end of unur_distr_cont_upd_pdfarea() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_cont_get_pdfarea( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get area below p.d.f. of distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   area below p.d.f. of distribution                                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, CONT, INFINITY );

  /* area known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PDFAREA) ) {
    /* try to compute area */
    if (DISTR.upd_area == NULL) {
      /* no function to compute area available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"area");
      return INFINITY;
    }
    else {
      /* compute area */
      unur_distr_cont_upd_pdfarea( distr );
    }
  }

  return DISTR.area;

} /* end of unur_distr_cont_get_pdfarea() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_cont_debug( const struct unur_distr *distr, const char *genid )
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
  COOKIE_CHECK(distr,CK_DISTR_CONT,RETURN_VOID);

  log = unur_get_stream();

  /* is this a derived distribution */
  if (distr->base) {
    switch (distr->id) {
    case UNUR_DISTR_CORDER:
      _unur_distr_corder_debug(distr,genid);
      return;
    default:
      _unur_warning(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return;
    }
  }

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tPDF with %d argument(s)\n",genid,DISTR.n_params);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(log,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(log,"%s:\tmode = %g\n",genid,DISTR.mode);
  else
    fprintf(log,"%s:\tmode unknown\n",genid);

  fprintf(log,"%s:\tdomain = (%g, %g)",genid,DISTR.domain[0],DISTR.domain[1]);
  _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);

  fprintf(log,"\n%s:\tarea below p.d.f. = %g",genid,DISTR.area);
  _unur_print_if_default(distr,UNUR_DISTR_SET_PDFAREA);
  fprintf(log,"\n%s:\n",genid);

} /* end of _unur_distr_cont_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_distr_cont_find_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* find mode of univariate continuous distribution numerically          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  struct unur_funct_generic fs;
  double mode;

  /* check arguments */
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, CONT, UNUR_ERR_DISTR_INVALID );
  if (DISTR.pdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PDF required for finding mode numerically"); 
    return UNUR_ERR_DISTR_DATA;
  }

  fs.f = (UNUR_FUNCT_GENERIC) DISTR.pdf;  
 
  mode = _unur_util_find_max(fs, DISTR.domain[0], DISTR.domain[1], DISTR.mode);

  if (!(_unur_FP_is_infinity( mode )) ){
    /* mode successfully computed */
    DISTR.mode = mode;
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE; 
  }
  else {
    /* computing mode did not work */
    /* (we do not change mode entry in distribution object) */
    return UNUR_ERR_DISTR_DATA;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_distr_cont_find_mode() */

/*---------------------------------------------------------------------------*/


#undef MAXIT
#undef SQRT_EPSILON

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
