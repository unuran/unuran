/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_codegen_JAVA.c                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Make generator code for TDR.                                         *
 *      (JAVA version)                                                       *
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

#include "codegen_source.h"
#include <tdr_source.h>

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_tdr_ps( FILE *out, 
		       struct unur_gen *gen, 
		       const char *rand_name, 
		       const char *pdf_name )
     /*----------------------------------------------------------------------*/
     /* code generator for method TDR variant PS (proportional squeeze)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   gen       ... pointer to generator object                          */
     /*   rand_name ... name of sampling routine                             */
     /*   pdf_name  ... name of PDF                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{

  struct unur_tdr_interval *iv;
  int i,j;

  /* check arguments */
  _unur_check_NULL("ACG",gen, 0);
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* check variant */
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
  case TDR_VARIANT_IA:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return 0;
  case TDR_VARIANT_PS:
    /* it should work! */
    break;
  }

  /* make section header for code file */
  fprintf(out,"\n");
  _unur_acg_JAVA_print_section_rule(out);
  _unur_acg_JAVA_print_section_line(out,"Sampling from %.35s distribution.",gen->distr.name);
  _unur_acg_JAVA_print_section_line(out,"Method: TDR - PS (Transformed Density Rejection / prop. squeeze)");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_acg_JAVA_print_section_line(out,"        Transformation = log(x) ... c = 0");         break;
  case TDR_VAR_T_SQRT:
    _unur_acg_JAVA_print_section_line(out,"        Transformation = -1/sqrt(x)  ... c = -1/2"); break;
  case TDR_VAR_T_POW:
    _unur_acg_JAVA_print_section_line(out,"        Transformation = -x^c  ... c = %g",GEN.c_T); break;
  }
  _unur_acg_JAVA_print_section_line(out,"        hat / squeeze ratio = %g",GEN.Atotal / GEN.Asqueeze);
  _unur_acg_JAVA_print_section_rule(out);
  fprintf(out,"\n");

  /* sampling routine */
  fprintf(out, "\tstatic double %s ()\n\t{\n",rand_name);

  /* constants */
  fprintf(out, "\t\t/* data */\n");

  /* guide table */
  fprintf(out, "\t\tfinal int guide_size = %d;\n", GEN.guide_size);
  fprintf(out, "\t\tfinal int guide[] = {\n\t\t\t");
  iv = GEN.iv;
  j = 0;
  for (i=0; i<GEN.guide_size; i++) {
    while (GEN.guide[i] != iv) {
      iv=iv->next;
      j++;
    }
    if (iv == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 0;
    }
    if (i==0)
      fprintf(out,"%d",j);
    else
      fprintf(out,(i%16)?", %d":",\n\t\t\t%d",j);
  }      
  fprintf(out," };\n");

  /* hat function */
  fprintf(out, "\t\tfinal double Atotal = %.20e;\n", GEN.Atotal);


  fprintf(out, "\t\tfinal IV[] iv = {\n");
  for ( iv=GEN.iv , i=0; i<GEN.n_ivs; iv=iv->next, i++){
    fprintf(out, "\t\t\tnew IV(\t"); /* zeilenklammer auf */
    fprintf(out, "%.20e, ",iv->x);             /* x         */
    switch (gen->variant & TDR_VARMASK_T) {    /* fx or Tfx */
    case TDR_VAR_T_LOG:
      fprintf(out, "%.20e, ",iv->fx);
      break;
    case TDR_VAR_T_SQRT:
      fprintf(out, "%.20e, ",iv->Tfx);
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 0;
    }
    fprintf(out, "\n\t\t\t\t");
    fprintf(out, "%.20e, ",iv->dTfx);           /* dTfx   */
    fprintf(out, "%.20e, ",iv->sq);             /* sq     */
    fprintf(out, "\n\t\t\t\t");
    fprintf(out, "%.20e, ",iv->Acum);           /* Acum   */
    fprintf(out, "%.20e ",iv->Ahatr);           /* Ahatr  */
    fprintf(out, (i<GEN.n_ivs-1)? "),\n" : ")\n"); /* zeilenklammer zu */
  }
  fprintf(out, "\t\t};\n\n");


  fprintf(out, "\t\t/* code */\n");

  /* declare variables */
  fprintf(out, "\t\tint I;\n");          /* random interval */
  fprintf(out, "\t\tdouble U;\n");       /* uniform random number */
  fprintf(out, "\t\tdouble V;\n");       /* random number for rejection */
  fprintf(out, "\t\tdouble X;\n");       /* non-uniform random variate */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out, "\t\tdouble t;\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out, "\t\tdouble Thx;\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  } /* end switch */
  fprintf(out, "\t\n");

  /* main loop */
  fprintf(out, "\t\twhile (1==1) {\n");

  /* sample from U(0,1) */
  fprintf(out, "\t\t\tU = Math.random();\n");

  /* look up in guide table and search for interval */
  fprintf(out, "\t\t\tI =  guide[(int) (U * guide_size)];\n");
  fprintf(out, "\t\t\tU *= Atotal;\n");
  fprintf(out, "\t\t\twhile (iv[I].Acum < U) { I++; }\n");

  /* reuse of uniform random number */
  fprintf(out, "\t\t\tU -= iv[I].Acum - iv[I].Ahatr;\n");
  /* result: U in (-A_hatl, A_hatr) */

  /* generate from hat distribution */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out, "\t\t\tt = iv[I].dTfx * U / iv[I].fx;\n");
    fprintf(out, "\t\t\tif (Math.abs(t) > 1.e-8) {\n");
    fprintf(out, "\t\t\t\tX = iv[I].x + Math.log(t + 1.) * U / (iv[I].fx * t); }\n");
    fprintf(out, "\t\t\telse {\n");
    fprintf(out, "\t\t\t\tX = iv[I].x + U / iv[I].fx * (1 - t/2.); }\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out, "\t\t\tX = iv[I].x + (U * iv[I].Tfx * iv[I].Tfx) / (1.-iv[I].Tfx*iv[I].dTfx*U);\n");
    break;
  } /* end switch */

  /* accept or reject */
  fprintf(out, "\t\t\tV = Math.random();\n");

  /* squeeze acceptance */
  fprintf(out, "\t\t\tif (V <= iv[I].sq) { return X; }\n");

  /* evaluate hat at X:
     get uniform random number between 0 and hat(X) */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out, "\t\t\tV *= iv[I].fx * Math.exp(iv[I].dTfx * (X - iv[I].x));\n");
    break;
  case TDR_VAR_T_SQRT:
     /* transformed hat at X */ 
    fprintf(out, "\t\t\tThx = iv[I].Tfx + iv[I].dTfx * (X - iv[I].x);\n");
    fprintf(out, "\t\t\tV /= Thx*Thx;\n");
    break;
  } /* end switch */

  /* main rejection */
  fprintf(out, "\t\t\tif (V <= %s(X)) { return X; }\n",pdf_name);

  /* end of loop */
  fprintf(out, "\t\t}\n");

  /* end of function */
  fprintf(out, "\t}\n");

  /* o.k. */
  return 1;

} /* end of _unur_acg_JAVA_tdr_ps() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_tdr_class_IV( FILE *out, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* code generator for method TDR variant PS -- class                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   gen       ... pointer to generator object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{

  /* check arguments */
  _unur_check_NULL("ACG",gen, 0);
  COOKIE_CHECK(gen,CK_TDR_GEN,0);

  /* check variant */
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
  case TDR_VARIANT_IA:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return 0;
  case TDR_VARIANT_PS:
    /* it should work! */
    break;
  }

  _unur_acg_JAVA_print_section_title(out, "Class for intervals.");

  /* class IV */
  fprintf(out,"\tprivate static class IV {\n");
  fprintf(out,"\t\tprivate double x, fx, Tfx, dTfx, sq, Acum, Ahatr;\n\n");

  fprintf(out, "\t\t/* Constructor */\n");
  fprintf(out, "\t\tpublic IV(double x, double fx_or_Tfx, double dTfx,\n");
  fprintf(out, "\t\t\t\tdouble sq, double Acum, double Ahatr)\n");
  fprintf(out, "\t\t{\n");

  fprintf(out,"\t\t\tthis.x\t\t= x;\n");
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\t\tthis.fx\t= fx_or_Tfx;\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\t\t\tthis.Tfx\t= fx_or_Tfx;\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  fprintf(out,"\t\t\tthis.dTfx\t= dTfx;\n");
  fprintf(out,"\t\t\tthis.sq\t\t= sq;\n");
  fprintf(out,"\t\t\tthis.Acum\t= Acum;\n");
  fprintf(out,"\t\t\tthis.Ahatr\t= Ahatr;\n");

  fprintf(out, "\t\t}\n\t}\n");

  /* o.k. */
  return 1;
} /* _unur_acg_JAVA_tdr_class_IV() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_begin_class( FILE *out, struct unur_gen *gen, const char *class_name )
     /*----------------------------------------------------------------------*/
     /* begin of class containing PDF (as member fuction) and the tdr method */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   gen       ... pointer to generator object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{

  _unur_acg_JAVA_print_section_title(out, "Class for generator.");
  fprintf(out,"public class %s {\n", class_name);

  /* o.k. */
  return 1;
} /* end of _unur_acg_JAVA_begin_class() */

/*---------------------------------------------------------------------------*/

int
_unur_acg_JAVA_end_class( FILE *out, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* end of class containing PDF (as member fuction) and the tdr method   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   gen       ... pointer to generator object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... on success                                                   */
     /*   0 ... on error                                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{

  fprintf(out, "\n} /* end of class Generator */\n");

  /* o.k. */
  return 1;
} /* end of _unur_acg_JAVA_end_class() */

/*---------------------------------------------------------------------------*/
