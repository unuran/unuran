/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      distribution.c                                               *
 *                                                                           *
 *   manipulate distribution obejct                                          *
 *                                                                           *
 *   PARAMETER: struct unur_distr *                                          *
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

#include <unur_methods.h>

#include <unur_cookies.h>
#include <unur_errno.h>
#include <unur_math.h>
#include <unur_utils.h>

/*---------------------------------------------------------------------------*/

static char unknown_distr_name[] = "Unknown";

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_new_cont( void )
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
#define DISTR distr->data.cont
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

  /* name of distribution */
  distr->name = unknown_distr_name;

  /* set defaults                                                            */
  DISTR.pdf       = NULL;          /* pointer to p.d.f.                      */
  DISTR.dpdf      = NULL;          /* pointer to derivative of p.d.f.        */
  DISTR.cdf       = NULL;          /* pointer to c.d.f.                      */

  DISTR.n_params  = 0;             /* number of parameters of the pdf        */
  /* initialize parameters of the p.d.f.                                     */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.mode      = 0.;            /* location of mode                       */
  DISTR.area      = 1.;            /* area below p.d.f.                      */
  DISTR.domain[0] = -INFINITY;     /* left boundary of domain                */
  DISTR.domain[1] = INFINITY;      /* right boundary of domain               */

  distr->set = 0u;                 /* no parameters set                      */
  
  /* return pointer to object */
  return distr;

#undef DISTR
} /* end of unur_distr_new_cont() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_pdf( struct unur_distr *distr, void *pdf )
     /*----------------------------------------------------------------------*/
     /* set p.d.f. of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);
  CHECK_NULL(pdf,0);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.pdf = pdf;
    return 1;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }

} /* end of unur_distr_set_pdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_dpdf( struct unur_distr *distr, void *dpdf )
     /*----------------------------------------------------------------------*/
     /* set derivative of p.d.f. of distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to derivative of p.d.f.                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);
  CHECK_NULL(dpdf,0);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.dpdf = dpdf;
    return 1;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }
} /* end of unur_distr_set_dpdf() */


/*---------------------------------------------------------------------------*/

int
unur_distr_set_cdf( struct unur_distr *distr, void *cdf )
     /*----------------------------------------------------------------------*/
     /* set p.d.f. of distribution                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pdf   ... pointer to p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);
  CHECK_NULL(cdf,0);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.cdf = cdf;
    return 1;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }
} /* end of unur_distr_set_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_params( struct unur_distr *distr, double *params, int n_params )
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
  int i;

  /* check arguments */
  CHECK_NULL(distr,0);

  /* check new parameter for generator */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_warning(NULL,UNUR_ERR_GENERIC,"invalid number parameter");
    return 0;
  }

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.n_params = n_params;
    for (i=0; i < n_params; i++)
      distr->data.cont.params[i] = params[i];
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PARAMS;

  /* o.k. */
  return 1;
} /* end of unur_distr_set_param() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_mode( struct unur_distr *distr, double mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of p.d.f.                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.mode = mode;
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return 1;
} /* end of unur_distr_set_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_pdfarea( struct unur_distr *distr, double area )
     /*----------------------------------------------------------------------*/
     /* set area below p.d.f.                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   area  ... area below p.d.f.                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.area = area;
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }
  /* changelog */
  distr->set |= UNUR_DISTR_SET_PDFAREA;

  /* o.k. */
  return 1;

} /* end of unur_distr_set_pdfarea() */

/*---------------------------------------------------------------------------*/

int
unur_distr_set_domain( struct unur_distr *distr, double left, double right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,0);

  /* check new parameter for generator */
  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_SET_INVALID,"domain, left >= right");
    return 0;
  }

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,0);
    distr->data.cont.domain[0] = left;
    distr->data.cont.domain[1] = right;
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
    return 0;
  }

  /* changelog */
  distr->set |= UNUR_DISTR_SET_DOMAIN;

  /* o.k. */
  return 1;

} /* end of unur_distr_set_domain() */

/*---------------------------------------------------------------------------*/

void
unur_distr_copy( struct unur_distr *distr1, struct unur_distr *distr2 )
     /*----------------------------------------------------------------------*/
     /* copy distribution object distr2 into distr1                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr1 ... pointer to target object                                */
     /*   distr2 ... pointer to source object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr1,/*void*/);
  CHECK_NULL(distr2,/*void*/);

  /* copy main structure */
  memcpy( distr1, distr2, sizeof( struct unur_distr ) );

} /* end of unur_distr_copy() */

/*---------------------------------------------------------------------------*/

void
unur_distr_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(distr,/*void*/);

  switch (distr->type) {
  case UNUR_DISTR_CONT:
    COOKIE_CHECK(distr,CK_DISTR_CONT,/*void*/);
    break;
  default:
    _unur_warning(NULL,UNUR_ERR_UNKNOWNDISTR,"");
  }

  free( distr );

} /* end of unur_distr_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

void
_unur_distr_debug_cont( struct unur_distr *distr, char *genid )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into logfile                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*----------------------------------------------------------------------*/
{
#define DISTR distr->data.cont

  FILE *log;
  int i;

  /* check arguments */
  CHECK_NULL(distr,/*void*/);
  COOKIE_CHECK(distr,CK_DISTR_CONT,/*void*/);

  log = unur_get_stream();

  fprintf(log,"%s: distribution:\n",genid);
  fprintf(log,"%s:\ttype = continuous univariate distribution\n",genid);
  fprintf(log,"%s:\tname = %s\n",genid,distr->name);

  fprintf(log,"%s:\tp.d.f with %d argument(s)\n",genid,DISTR.n_params);
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

#undef DISTR
} /* end of _unur_distr_debug_cont() */

/*---------------------------------------------------------------------------*/
