/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      codegen_UNURAN.c                                             *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ACG (Automatic Code Generator)                                       *
 *      (UNURAN version)                                                     *
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

#include <stdarg.h>
#include "codegen_source.h"

/*---------------------------------------------------------------------------*/

int
unur_acg_UNURAN( struct unur_gen *gen, FILE *out, const char *distr_name )
     /*----------------------------------------------------------------------*/
     /* Automatic code generator (C version)                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   out        ... output stream                                       */
     /*   distr_name ... name of distribution                                */
     /*                  (used to name routines, if NULL the UNURAN          */
     /*                   build-in name is used.)                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  char *rand_name;                /* names for function */
  int return_code;                /* exit status of routine */

  /* check arguments */
  _unur_check_NULL("unur_acg", gen, 0 );

  /* make name of PDF function and sampling routine */
  if (distr_name == NULL) 
    distr_name = unur_distr_get_name( &(gen->distr) );

  rand_name = _unur_malloc((6+strlen(distr_name)) * sizeof(char));
  sprintf(rand_name,"rand_%s",distr_name);

  /* make code */
  switch (gen->method) {
  case UNUR_METH_TDR:
    return_code =
      _unur_acg_UNURAN_header( &(gen->distr), out, rand_name ) &&
      _unur_acg_UNURAN_tdr_ps( gen, out, rand_name );
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return_code = 0;
  }

  /* clear memory */
  free(rand_name);

  /* make error message in source file if generation failed */
  if (return_code == 0)
    fprintf(out,"\n#error Sorry. Could not make generator code!!\n\n");

  /* end */
  _unur_acg_C_print_section_title( out, "End of Generator" );

  return return_code;

} /* end of unur_acg_UNURAN() */

/*---------------------------------------------------------------------------*/
