/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: testcounter.c                                                     *
 *                                                                           *
 *   Count function evaluations                                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
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

#include "testunuran.h"

#include <unur_source.h>
#include <distr/distr_source.h>

/*---------------------------------------------------------------------------*/
/* count number of function evaluations                                      */

/* counter for evaluations of PDF and similar functions */
static int counter_pdf = 0;

/* pointer to original functions */
static double (*cont_pdf_to_use)(double x, const struct unur_distr *distr);

/* wrapper for functions */
static double cont_pdf_with_counter( double x, const struct unur_distr *distr ) {
  ++counter_pdf; return cont_pdf_to_use(x,distr); }

/* set and start counter for PDF in parameter object */
int start_counter_pdf( UNUR_PAR *par )
{
  struct unur_distr *distr;

/*   distr = (struct unur_distr *) par->distr; */

  /* make a copy (clone) of the distribution object */
  distr = _unur_distr_clone(par->distr);

  /* set pointer to distribution object in parameter object to cloned one */
  par->distr = distr;
  par->distr_is_privatecopy = FALSE;

  /* exchange pointer to PDF etc. with counting wrapper */
  switch (distr->type) {
  case UNUR_DISTR_CONT:
    /* this only works for continuous distributions */
    if (distr->data.cont.pdf == cont_pdf_with_counter) {
      /* we must avoid this situation as it causes infinite recursion */
      return UNUR_FAILURE;
    } 

    cont_pdf_to_use = distr->data.cont.pdf; 
    distr->data.cont.pdf = cont_pdf_with_counter;
    break;

  default:
    return UNUR_FAILURE;
  }

  /* reset counter */
  counter_pdf = 0;

  return UNUR_SUCCESS;
}

/* get number of PDF evaluations */
int get_counter_pdf(void)
{
  return counter_pdf;
}

/* reset counter to 0 */
void reset_counter_pdf(void)
{
  counter_pdf = 0;
}

/*---------------------------------------------------------------------------*/
