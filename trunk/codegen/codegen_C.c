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
#include <methods/unur_methods_source.h>

/*---------------------------------------------------------------------------*/

int
unur_acg_C( const UNUR_GEN *gen, FILE *out, const char *distr_name, int with_main )
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

  /* check arguments */
  _unur_check_NULL("unur_acg", gen, 0 );

  /* make name of PDF function and sampling routine */
  if (distr_name == NULL) 
    distr_name = unur_distr_get_name( gen->distr );

  pdf_name = _unur_malloc((5+strlen(distr_name)) * sizeof(char));
  sprintf(pdf_name,"pdf_%s",distr_name);

  rand_name = _unur_malloc((6+strlen(distr_name)) * sizeof(char));
  sprintf(rand_name,"rand_%s",distr_name);

  /* make code */
  switch (gen->method) {
  case UNUR_METH_TDR:
    return_code =
      _unur_acg_C_header    ( out, gen->distr, rand_name ) &&
      _unur_acg_C_demo_urng ( out ) &&
      _unur_acg_C_PDF       ( out, gen->distr, pdf_name ) &&
      _unur_acg_C_tdr_ps    ( out, gen, rand_name, pdf_name ) &&
      _unur_acg_C_print_section_title( out, "End of Generator" );
    if (with_main && return_code)
      _unur_acg_C_main( out, rand_name );
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
  _unur_acg_C_print_section_title( out, "Uniform (pseudo-)random number generator");
      
  fprintf(out,"/* Define the uniform (pseudo-)random number generator              */\n");
  fprintf(out,"/* of your choice here.                                             */\n\n");
  fprintf(out,"/* #define uniform()   your_uniform_rng() */\n");
  fprintf(out,"\n");

  fprintf(out,"\n");
  _unur_acg_C_print_section_rule(out);
  _unur_acg_C_print_section_line(out,"LCG (Linear Congruential Generator)      by Park & Miller (1988)");
  _unur_acg_C_print_section_line(out,"  x_(n+1) = 16807 * x_n mod (2^32 - 1)        [Minimal Standard]");
  _unur_acg_C_print_section_line(out,"");
  _unur_acg_C_print_section_line(out,"WARNING! This short build-in uniform random number generator");
  _unur_acg_C_print_section_line(out,"is not state-of-the-art and should not be used for simulations.");
  _unur_acg_C_print_section_line(out,"It should be replaced by a modern generator of your choice");
  _unur_acg_C_print_section_line(out,"with a (much) longer period (see above).");
  _unur_acg_C_print_section_line(out,"E.g. Mersenne Twister by Makoto Matsumoto and Takuji Nishimura,");
  _unur_acg_C_print_section_line(out,"see http://www.math.keio.ac.jp/~matumoto/emt.html");
  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  /*************************************************************
   * Park and Miller (1988).                                   *
   * Random number generators: good ones are hard to find.     *
   * Comm. ACM 31, pp. 1192--1201.                             *
   *************************************************************/
  
  fprintf(out,"#ifndef uniform\n");
  fprintf(out,"#define uniform() _uniform_demo()\n");
  fprintf(out,"\n");
  fprintf(out,"static double _uniform_demo (void)\n");
  fprintf(out,"{\n");
  fprintf(out,"\tstatic unsigned int x = %d;  /* seed  */\n",rand());
  fprintf(out,"\tint hi, lo, test;\n\n");

  fprintf(out,"#\tdefine a 16807       /* multiplicator */\n");
  fprintf(out,"#\tdefine m 2147483647  /* modulus */\n");
  fprintf(out,"#\tdefine q 127773      /* m / a */\n");
  fprintf(out,"#\tdefine r 2836        /* m %% a */\n\n");

  fprintf(out,"\thi = x / q;\n");
  fprintf(out,"\tlo = x %% q;\n");
  fprintf(out,"\ttest = a * lo - r * hi;\n");
  fprintf(out,"\tx = (test > 0 ) ? test : test + m;\n");
  fprintf(out,"\treturn (x * 4.656612875245796924105750827e-10);\n\n");

  fprintf(out,"#\tundef a\n");
  fprintf(out,"#\tundef m\n");
  fprintf(out,"#\tundef q\n");
  fprintf(out,"#\tundef r\n");
  fprintf(out,"}\n");
  fprintf(out,"#endif /* uniform */\n\n");

  return 1;
} /* end of _unur_acg_C_demo_urng() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_C_print_section_title( FILE *out, const char *title )
     /*----------------------------------------------------------------------*/
     /* print a section header with title to output stream                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   title     ... section title                                        */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"\n");
  _unur_acg_C_print_section_rule(out);
  _unur_acg_C_print_section_line(out,title);
  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_C_print_section_title() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_C_print_section_rule( FILE *out )
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
} /* end of _unur_acg_C_print_section_rule() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_C_print_section_line( FILE *out, const char *format, ... )
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
} /* end of _unur_acg_C_print_section_line() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_C_main( FILE *out, const char *rand_name )
     /*----------------------------------------------------------------------*/
     /* Print main (C version)                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   rand_name ... name of sampling routine                             */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  _unur_acg_C_print_section_title( out, "Main");

  fprintf(out,"\n");
  fprintf(out,"#include<stdio.h>\n");
  fprintf(out,"\n");
      
  fprintf(out,"int main()\n{\n");
  fprintf(out,"\tint i;\n");
  fprintf(out,"\tdouble x;\n\n");
  fprintf(out,"\tfor (i=0; i<10; i++) {\n");
  fprintf(out,"\t\tx = %s();\n",rand_name);
  fprintf(out,"\t\tprintf(\"%%g\\n\",x);\n");
  fprintf(out,"\t}\n");
  fprintf(out,"\treturn 1;\n");
  fprintf(out,"}\n");

  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_C_main() */

/*---------------------------------------------------------------------------*/

