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
 *      (This is only a version for testing ACG. It should not be used.)     *
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
unur_acg_UNURAN( struct unur_gen *gen, FILE *out, const char *distr_name, int with_main )
     /*----------------------------------------------------------------------*/
     /* Automatic code generator (C version)                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   out        ... output stream                                       */
     /*   distr_name ... name of distribution                                */
     /*                  (used to name routines, if NULL the UNURAN          */
     /*                   build-in name is used.)                            */
     /*   with_main  ... whether to include main into source                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  char *pdf_name, *rand_name;     /* names for function */
  int return_code;                /* exit status of routine */

  int n_cpoints = 30;

  /* check arguments */
  _unur_check_NULL("unur_acg", gen, 0 );

  /* make name of PDF function and sampling routine */
  if (distr_name == NULL) 
    distr_name = unur_distr_get_name( &(gen->distr) );

  pdf_name = _unur_malloc((5+strlen(distr_name)) * sizeof(char));
  sprintf(pdf_name,"pdf_%s",distr_name);

  rand_name = _unur_malloc((6+strlen(distr_name)) * sizeof(char));
  sprintf(rand_name,"rand_%s",distr_name);

  /* make code */
  switch (gen->method) {
  case UNUR_METH_TDR:
    return_code =
      _unur_acg_UNURAN_header( &(gen->distr), out, rand_name ) &&
      _unur_acg_UNURAN_PDF( &(gen->distr), out, pdf_name ) &&
      _unur_acg_UNURAN_tdr_ps( gen, out, rand_name, n_cpoints );
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return_code = 0;
  }

  /* clear memory */
  free(pdf_name);
  free(rand_name);

  /* make error message in source file if generation failed */
  if (return_code == 0)
    fprintf(out,"\n#error Sorry. Could not make generator code!!\n\n");

  /* end */
  _unur_acg_C_print_section_title( out, "End of Generator" );

  return return_code;

} /* end of unur_acg_UNURAN() */

/*---------------------------------------------------------------------------*/

int 
_unur_acg_UNURAN_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf_name)
     /*----------------------------------------------------------------------*/
     /* Code generator for PDFs of UNURAN build-in standard distributions.   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr      ... pointer to distribution object                      */
     /*   out        ... output stream                                       */
     /*   pdf_name   ... name of pdf function                                */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  double *fpar;
  int n_fpar;
  int i;

  /* make section header for code file */
  fprintf(out,"\n");
  _unur_acg_C_print_section_rule(out);
  _unur_acg_C_print_section_line(out," PDF for %s distribution.",distr->name);
  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  /* include header files */
  fprintf(out,"#include <unuran.h>\n\n");

  /* PDF routine */
  fprintf(out,"double %s (double x)\n{\n",pdf_name);

  /* declare variables */
  fprintf(out,"\tstatic UNUR_DISTR *distr = NULL;\n"); 
  fprintf(out,"\n");

  /* set up */
  fprintf(out,"\tif (distr == NULL) {\n"); 

  /* make distribution object */
  n_fpar = unur_distr_cont_get_pdfparams(distr,&fpar);
  if (n_fpar == 0)
    fprintf(out,"\t\tdouble *fpar = NULL;\n");
  else {
    fprintf(out,"\t\tdouble fpar[] = {\n");
    for (i=0; i<n_fpar-1; i++)
      fprintf(out,"\t\t\t%.20e,\n",fpar[i]);
    fprintf(out,"\t\t\t%.20e };\n",fpar[n_fpar-1]);
  }

  fprintf(out,"\t\tdistr = unur_distr_%s(fpar,%d);\n",distr->name,n_fpar);

  /* emergency exit */
  fprintf(out,"\t\tif (distr == NULL) {\n");
  fprintf(out,"\t\t\tfprintf(stderr, \"Cannot create distribution object.\\n\");\n");
  fprintf(out,"\t\t\texit (EXIT_FAILURE);\n");
  fprintf(out,"\t\t}\n\n");

  fprintf(out,"\t}\n\n");

  /* evaluate PDF */
  fprintf(out,"\treturn unur_distr_cont_eval_pdf( x, distr );\n");

  fprintf(out,"}\n\n");

  return 1;
} /* end of _unur_acg_UNURAN_PDF() */

/*---------------------------------------------------------------------------*/
