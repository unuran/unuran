/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      codegen_C.c                                                  *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ACG (Automatic Code Generator)                                       *
 *      (C version)                                                          *
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
unur_acg_C( struct unur_gen *gen, FILE *out, const char *distr_name )
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
  char *pdf_name, *rand_name;     /* names for function */
  int return_code;                /* exit status of routine */

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
      _unur_acg_C_demo_urng( out ) &&
      _unur_acg_C_PDF( &(gen->distr), out, pdf_name ) &&
      _unur_acg_C_tdr_ps( gen, out, rand_name, pdf_name );
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
  _unur_acg_C_print_sectionheader( out, 1, "End of Generator" );

  return return_code;

} /* end of unur_acg_C() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_C_demo_urng( FILE *out )
     /*----------------------------------------------------------------------*/
     /* uniform random number generator (for demo mode only)                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  _unur_acg_C_print_sectionheader
    ( out, 1, 
      "Uniform (pseudo-)random number generator"
      );
      
  fprintf(out,"/* Define the uniform (pseudo-)random number generator              */\n");
  fprintf(out,"/* of your choice here.                                             */\n\n");
  fprintf(out,"/* #define uniform()   your_uniform_rng() */\n");
  fprintf(out,"\n");

  _unur_acg_C_print_sectionheader
    ( out, 9,
      "LCG (Linear Congruential Generator) by G. Marsaglia (1972).",
      "  x_(n+1) = 69069 * x_n + 1 mod 2^32",
      "",
      "WARNING! This short build-in uniform random number generator",
      "is not state-of-the-art and should not be used for simulations.",
      "It should be replaced by a modern generator of your choice",
      "with a (much) longer period (see above).",
      "E.g. Mersenne Twister by Makoto Matsumoto and Takuji Nishimura,",
      "see http://www.math.keio.ac.jp/~matumoto/emt.html" 
      );

  /*************************************************************
   * Marsaglia G. (1972), m = 2^32, a = 69069, c = 1           *
   * The structure of linear congruential sequences, in:       *
   * Applications of Number Theory to Numerical Analysis, S.K. *
   * Zaremba, ed., Academic Press, New York 1972               *
   *************************************************************/
  
  fprintf(out,"#ifndef uniform\n");
  fprintf(out,"#define uniform() _uniform_demo()\n");
  fprintf(out,"\n");
  fprintf(out,"static double _uniform_demo (void)\n");
  fprintf(out,"{\n");
  fprintf(out,"\tstatic unsigned long int x = %d;   /* seed  */\n",rand());
  fprintf(out,"\tx = 69069 * x + 1;\n");
  fprintf(out,"\treturn (x * %.20e + %.20e);\n", pow(2.,-32), pow(2.,-33));
  fprintf(out,"}\n");
  fprintf(out,"#endif\n\n");

  return 1;
} /* end of _unur_acg_C_demo_urng() */

/*---------------------------------------------------------------------------*/

void
_unur_acg_C_print_sectionheader( FILE *out, int n_lines, ... )
     /*----------------------------------------------------------------------*/
     /* print a section header with n_lines lines to output stream           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   n_lines   ... number of lines                                      */
     /*   ...       ... (optional) arguments to be be printed                */
     /*----------------------------------------------------------------------*/
{
  char *string;     /* string to printed in section header   */
  va_list ap;       /* pointer to variable list of arguments */
  const char hrule[] = "/* ---------------------------------------------------------------- */\n";

  /* start of variable parameter list */
  va_start(ap, n_lines);

  /* write into output stream */
  fprintf (out, "\n");
  fprintf (out, hrule);
  for (; n_lines>0; n_lines--) {
    string = va_arg( ap, char* );
    fprintf(out,"/* %-64.64s */\n",string);
  }
  fprintf (out, hrule);
  fprintf (out,"\n");
        
  /* end of variable parameter list */
  va_end(ap);

} /* end of _unur_acg_C_print_sectionheader() */

/*---------------------------------------------------------------------------*/
