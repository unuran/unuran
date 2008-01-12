/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: distr_info.c                                                      *
 *                                                                           *
 *   create info strings for distribution objects                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
#include "distr.h"
#include "distr_source.h"

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_distr_info_typename( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* type and name of the given generator object.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;

  /* name of distribution */
  _unur_string_append(info,"   name      = %s\n", distr->name);

  /* type of distribution */
  _unur_string_append(info,"   type      = ");
  switch (distr->type) {
  case UNUR_DISTR_CONT:
    _unur_string_append(info,"continuous univariate distribution\n");
    return;
  case UNUR_DISTR_CEMP:
    _unur_string_append(info,"continuous empirical univariate distribution\n");
    return;
  case UNUR_DISTR_CVEC:
    _unur_string_append(info,"continuous multivariate distribution\n");
    return;
  case UNUR_DISTR_CVEMP:
    _unur_string_append(info,"continuous empirical multivariate distribution\n");
    return;
  case UNUR_DISTR_DISCR:
    _unur_string_append(info,"discrete univariate distribution\n");
    return;
  case UNUR_DISTR_MATR:
    _unur_string_append(info,"matrix distribution\n");
    return;
  default:
    _unur_error(distr->name,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  }

} /* end of unur_distr_info_typename() */


/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
