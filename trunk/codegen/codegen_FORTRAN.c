/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      codegen_FORTRAN.c                                            *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      ACG (Automatic Code Generator)                                       *
 *      (FORTRAN version)                                                    *
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
unur_acg_FORTRAN( const UNUR_GEN *gen, FILE *out, const char *distr_name, int with_main )
     /*----------------------------------------------------------------------*/
     /* Automatic code generator (FORTRAN version)                           */
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
    distr_name = unur_distr_get_name( &(gen->distr) );

  pdf_name = _unur_malloc((5+strlen(distr_name)) * sizeof(char));
  sprintf(pdf_name,"f%.5s",distr_name);

  rand_name = _unur_malloc((6+strlen(distr_name)) * sizeof(char));
  sprintf(rand_name,"r%.5s",distr_name);

  /* make code */
  switch (gen->method) {
  case UNUR_METH_TDR:
    return_code =
      _unur_acg_FORTRAN_header    ( out, &(gen->distr), rand_name ) &&
      _unur_acg_FORTRAN_demo_urng ( out ) &&
      _unur_acg_FORTRAN_PDF       ( out, &(gen->distr), pdf_name ) &&
      _unur_acg_FORTRAN_tdr_ps    ( out, gen, rand_name, pdf_name ) &&
      _unur_acg_FORTRAN_print_section_title( out, "End of Generator" );
    if (with_main && return_code)
      _unur_acg_FORTRAN_main( out, rand_name );
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

} /* end of unur_acg_FORTRAN() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_demo_urng( FILE *out )
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
  _unur_acg_FORTRAN_print_section_line( out, "Uniform (pseudo-)random number generator");
      
  fprintf(out,"C  Replace the uniform (pseudo-)random number generator 'urand()'\n");
  fprintf(out,"C  by a generator of your choice.\n\n");
  fprintf(out,"\n");

  fprintf(out,"\n");
  _unur_acg_FORTRAN_print_section_rule(out);
  _unur_acg_FORTRAN_print_section_line(out,"LCG (Linear Congruential Generator)        by Park & Miller (1988)");
  _unur_acg_FORTRAN_print_section_line(out,"  x_(n+1) = 16807 * x_n mod (2^32 - 1)          [Minimal Standard]");
  _unur_acg_FORTRAN_print_section_line(out,"");
  _unur_acg_FORTRAN_print_section_line(out,"WARNING! This short build-in uniform random number generator");
  _unur_acg_FORTRAN_print_section_line(out,"is not state-of-the-art and should not be used for simulations.");
  _unur_acg_FORTRAN_print_section_line(out,"It should be replaced by a modern generator of your choice");
  _unur_acg_FORTRAN_print_section_line(out,"with a (much) longer period (see above).");
  _unur_acg_FORTRAN_print_section_line(out,"E.g. Mersenne Twister by Makoto Matsumoto and Takuji Nishimura,");
  _unur_acg_FORTRAN_print_section_line(out,"see http://www.math.keio.ac.jp/~matumoto/emt.html");
  _unur_acg_FORTRAN_print_section_rule(out);
  fprintf(out,"\n");

  /*************************************************************
   * Park and Miller (1988).                                   *
   * Random number generators: good ones are hard to find.     *
   * Comm. ACM 31, pp. 1192--1201.                             *
   *************************************************************/
  
  fprintf(out,"      DOUBLE PRECISION FUNCTION urand()\n\n");

  fprintf(out,"      INTEGER a, m, q, r, x, hi, lo, test\n");
  fprintf(out,"      PARAMETER (a = 16807)\n");
  fprintf(out,"      PARAMETER (m = 2147483647)\n");
  fprintf(out,"      PARAMETER (q = 127773)\n");
  fprintf(out,"      PARAMETER (r = 2836)\n\n");

  fprintf(out,"C     state variable\n");
  fprintf(out,"      DATA x/173923/\n");
  fprintf(out,"      SAVE x\n\n");

  fprintf(out,"      hi = x / q\n");
  fprintf(out,"      lo = MOD(x,q)\n");
  fprintf(out,"      test = a * lo - r * hi\n");
  fprintf(out,"      IF (test .gt. 0) THEN\n");
  fprintf(out,"          x = test\n");
  fprintf(out,"      ELSE\n");
  fprintf(out,"          x = test + m\n");
  fprintf(out,"      END IF\n\n");

  fprintf(out,"      urand = x * 4.656612875245796924105750827D-10\n\n");

  fprintf(out,"      END\n\n");


  return 1;
} /* end of _unur_acg_FORTRAN_demo_urng() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_print_section_title( FILE *out, const char *title )
     /*----------------------------------------------------------------------*/
     /* print a section header with title to output stream                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   title     ... section title                                        */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"\n");
  _unur_acg_FORTRAN_print_section_rule(out);
  _unur_acg_FORTRAN_print_section_line(out,title);
  _unur_acg_FORTRAN_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_FORTRAN_print_section_title() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_print_section_rule( FILE *out )
     /*----------------------------------------------------------------------*/
     /* print a rule for section header                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*----------------------------------------------------------------------*/
{
  const char hrule[] = "* ------------------------------------------------------------------ *\n";
  fprintf (out, hrule);

  return 1;
} /* end of _unur_acg_FORTRAN_print_section_rule() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_print_section_line( FILE *out, const char *format, ... )
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
  fprintf(out,"* %-66.66s *\n",buffer);
        
  /* end of variable parameter list */
  va_end(ap);

  return 1;
} /* end of _unur_acg_FORTRAN_print_section_line() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_print_double( FILE *out, double x )
     /*----------------------------------------------------------------------*/
     /* print a double number in FORTRAN format                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out ... output stream                                              */
     /*   x   ... double precision number                                    */
     /*----------------------------------------------------------------------*/
{
  char buf[256];
  char *here_is_e;

  /* C format */
  sprintf(buf,"%.20e",x);
  
  /* replace `e' by `D' */
  here_is_e = strchr(buf, 'e');
  if (here_is_e)
    *here_is_e = 'D';
  
  /* print */
  fprintf(out,"%s",buf);

  return 1;
} /* end of _unur_acg_FORTRAN_print_double() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_main( FILE *out, const char *rand_name )
     /*----------------------------------------------------------------------*/
     /* Print main (FORTRAN version)                                         */
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
  _unur_acg_FORTRAN_print_section_title( out, "Main");

  fprintf(out,"      PROGRAM demo\n\n");
  fprintf(out,"      IMPLICIT DOUBLE PRECISION (a-h,o-z)\n\n");
  fprintf(out,"      DO 1 i=1,10\n");
  fprintf(out,"         x = %s()\n",rand_name);
  fprintf(out,"         PRINT *, x\n");
  fprintf(out,"1     CONTINUE\n\n");
  fprintf(out,"      END\n\n");
  _unur_acg_FORTRAN_print_section_rule(out);
  fprintf(out,"\n");

  return 1;
} /* end of _unur_acg_FORTRAN_main() */

/*---------------------------------------------------------------------------*/
