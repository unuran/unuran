/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      codegen_JAVA.c                                               *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ACG (Automatic Code Generator)                                       *
 *      (JAVA Version)                                                       *
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
unur_acg_JAVA( struct unur_gen *gen, FILE *out, const char *distr_name, int with_main )
     /*----------------------------------------------------------------------*/
     /* Automatic code generator (JAVA version)                              */
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
  char pdf_name[] = "pdf";        /* name for pdf function */
  char rand_name[] = "sample";    /* name for sampling routine */
  char *class_name;               /* name for generator class */
  int return_code;                /* exit status of routine */

  /* check arguments */
  _unur_check_NULL("unur_acg", gen, 0 );

  /* make name for generator class */
  if (distr_name == NULL) 
    distr_name = unur_distr_get_name( &(gen->distr) );

  class_name = _unur_malloc((15+strlen(distr_name)) * sizeof(char));
  sprintf(class_name,"Generator_%s",distr_name);

  /* make code */
  switch (gen->method) {
  case UNUR_METH_TDR:
    return_code =
      _unur_acg_JAVA_header       ( out, &(gen->distr), class_name ) &&
      _unur_acg_JAVA_begin_class  ( out, gen, class_name ) &&
      _unur_acg_JAVA_urng         ( out ) &&
      _unur_acg_JAVA_tdr_class_IV ( out, gen ) &&
      _unur_acg_JAVA_PDF          ( out, &(gen->distr), pdf_name ) &&
      _unur_acg_JAVA_tdr_ps       ( out, gen, rand_name, pdf_name ) &&
      _unur_acg_JAVA_end_class    ( out, gen ) &&
      _unur_acg_JAVA_print_section_title( out, "End of Generator" );
    if (with_main && return_code)
      _unur_acg_JAVA_main( out, class_name );
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return_code = 0;
  }

  /* clear memory */
  free(class_name);

  /* make error message in source file if generation failed */
  if (return_code == 0)
    fprintf(out,"\n#error Sorry. Could not make generator code!!\n\n");

  return return_code;

} /* end of unur_acg_JAVA() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_urng( FILE *out )
     /*----------------------------------------------------------------------*/
     /* uniform random number generator                                      */
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
  _unur_acg_JAVA_print_section_title( out, "Uniform (pseudo-)random number generator");
      
  fprintf(out,"\t// Define the uniform (pseudo-)random number generator\n");
  fprintf(out,"\t// of your choice here by operator overload.\n");
  fprintf(out,"\t//\n");
  fprintf(out,"\t// Otherwise the built-in generator is used which might\n");
  fprintf(out,"\t// NOT be state-of-the-art. A good choice is e.g. the\n");
  fprintf(out,"\t// Mersenne Twister by Makoto Matsumoto and Takuji Nishimura,\n");
  fprintf(out,"\t// see http://www.math.keio.ac.jp/~matumoto/emt.html.\n");
  fprintf(out,"\t//\n");
  fprintf(out,"\t// double Random ()\n");
  fprintf(out,"\t// {\n");
  fprintf(out,"\t//\t...\n");
  fprintf(out,"\t//\treturn U;\n");
  fprintf(out,"\t// }\n");

  return 1;
} /* end of _unur_acg_JAVA_urng() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_print_section_title( FILE *out, const char *title )
     /*----------------------------------------------------------------------*/
     /* print a section header with title to output stream                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   title     ... section title                                        */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"\n");
  _unur_acg_JAVA_print_section_rule(out);
  _unur_acg_JAVA_print_section_line(out,title);
  _unur_acg_JAVA_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_JAVA_print_section_title() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_print_section_rule( FILE *out )
     /*----------------------------------------------------------------------*/
     /* print a rule for section header                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*----------------------------------------------------------------------*/
{
  const char hrule[] = "/* ---------------------------------------------------------------- */\n";
  fprintf (out, hrule);

  return 1;
} /* end of _unur_acg_JAVA_print_section_rule() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_print_section_line( FILE *out, const char *format, ... )
     /*----------------------------------------------------------------------*/
     /* print a section header with n_lines lines to output stream           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   format    ... format for output                                    */
     /*   ...       ... (optional) arguments to be be printed                */
     /*----------------------------------------------------------------------*/
{
  char buffer[256];
  va_list ap;       /* pointer to variable list of arguments */

  /* start of variable parameter list */
  va_start(ap, format);

  /* write into output stream */
  vsprintf(buffer,format,ap);
  fprintf(out,"/* %-64.64s */\n",buffer);
        
  /* end of variable parameter list */
  va_end(ap);

  return 1;
} /* end of _unur_acg_JAVA_print_section_line() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_main( FILE *out, const char *class_name )
     /*----------------------------------------------------------------------*/
     /* Print main (JAVA version)                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out        ... output stream                                       */
     /*   class_name ... name of generator class                             */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  _unur_acg_JAVA_print_section_title( out, "Main");

  fprintf(out,"public class %s_demo {\n",class_name);
  fprintf(out,"\tpublic static void main(String[] args) throws Exception {\n\n");
  fprintf(out,"\t\tdouble x;\n\n");
  fprintf(out,"\t\tfor (int i=0; i<10; i++) {\n");
  fprintf(out,"\t\t\tx = %s.sample();\n",class_name);
  fprintf(out,"\t\t\tSystem.out.println(x);\n");
  fprintf(out,"\t\t}\n");
  fprintf(out,"\t}\n");
  fprintf(out,"}\n\n");


  _unur_acg_JAVA_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_JAVA_main() */

/*---------------------------------------------------------------------------*/
