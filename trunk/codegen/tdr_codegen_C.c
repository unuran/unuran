/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_codegen_C.c                                              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Make generator code for TDR.                                         *
 *      (C Version)                                                          *
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
_unur_acg_C_tdr_ps( struct unur_gen *gen, 
		    FILE *out, 
		    const char *rand_name, 
		    const char *pdf_name )
     /*----------------------------------------------------------------------*/
     /* code generator for method TDR variant PS (proportional squeeze)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   out       ... output stream                                        */
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
  char buffer[80], transf[80];

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
  sprintf(buffer,"Sampling from %.35s distribution.",gen->distr.name);
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    sprintf(transf,"        Transformtaion = log(x) ... c = 0");         break;
  case TDR_VAR_T_SQRT:
    sprintf(transf,"        Transformtaion = -1/sqrt(x)  ... c = -1/2"); break;
  case TDR_VAR_T_POW:
    sprintf(transf,"        Transformtaion = -x^c  ... c = %g",GEN.c_T); break;
  }
  _unur_acg_C_print_sectionheader
    ( out, 3, 
      buffer,
      "Method: TDR - PS (Transformed Density Rejection / prop. squeeze)",
      transf
      );

  /* sampling routine */
  fprintf(out,"double %s (void)\n{\n",rand_name);

  /* constants */
  fprintf(out,"\t/* data */\n");

  /* guide table */
  fprintf(out,"\tconst int guide_size = %d;\n", GEN.guide_size);
  fprintf(out,"\tconst int guide[%d] = {\n\t\t", GEN.guide_size);
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
      fprintf(out,(i%16)?", %d":",\n\t\t%d",j);
  }      
  fprintf(out," };\n");

  /* hat function */
  fprintf(out,"\tconst double Atotal = %.20e;\n", GEN.Atotal);

  fprintf(out,"\tconst struct {\n");
  fprintf(out,"\t\tdouble x;\n");
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\tdouble fx;\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\t\tdouble Tfx;\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  fprintf(out,"\t\tdouble dTfx;\n");
  fprintf(out,"\t\tdouble sq;\n");
  fprintf(out,"\t\tdouble Acum;\n");
  fprintf(out,"\t\tdouble Ahatr;\n");
  fprintf(out,"\t} iv[%d] = {\n",GEN.n_ivs);
  
  for (iv=GEN.iv; iv->next!=NULL; iv=iv->next) {
    fprintf(out,(iv==GEN.iv)?"\t\t{ ":",\n\t\t{ ");
    fprintf(out,"%.20e, ",iv->x);
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      fprintf(out,"%.20e, ",iv->fx);
      break;
    case TDR_VAR_T_SQRT:
      fprintf(out,"%.20e, ",iv->Tfx);
      break;
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return 0;
    }
    fprintf(out,"\n\t\t  ");
    fprintf(out,"%.20e, ",iv->dTfx);
    fprintf(out,"%.20e, ",iv->sq);
    fprintf(out,"\n\t\t  ");
    fprintf(out,"%.20e, ",iv->Acum);
    fprintf(out,"%.20e ",iv->Ahatr);
    fprintf(out,"}");
  }
  fprintf(out," };\n");

  fprintf(out,"\n");

  fprintf(out,"\t/* code */\n");

  /* declare variables */
  fprintf(out,"\tint I;\n");          /* random interval */
  fprintf(out,"\tdouble U;\n");       /* uniform random number */
  fprintf(out,"\tdouble V;\n");       /* random number for rejection */
  fprintf(out,"\tdouble X;\n");       /* non-uniform random variate */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\tdouble t;\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\tdouble Thx;\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  } /* end switch */
  fprintf(out,"\n");

  /* main loop */
  fprintf(out,"\twhile (1) {\n");

  /* sample from U(0,1) */
  fprintf(out,"\t\tU = uniform();\n");

  /* look up in guide table and search for interval */
  fprintf(out,"\t\tI =  guide[(int) (U * guide_size)];\n");
  fprintf(out,"\t\tU *= Atotal;\n");
  fprintf(out,"\t\twhile (iv[I].Acum < U) I++;\n");

  /* reuse of uniform random number */
  fprintf(out,"\t\tU -= iv[I].Acum - iv[I].Ahatr;\n");
  /* result: U in (-A_hatl, A_hatr) */

  /* generate from hat distribution */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\tt = iv[I].dTfx * U / iv[I].fx;\n");
    fprintf(out,"\t\tif (fabs(t) > 1.e-8)\n");
    fprintf(out,"\t\t\tX = iv[I].x + log(t + 1.) * U / (iv[I].fx * t);\n");
    fprintf(out,"\t\telse\n");
    fprintf(out,"\t\t\tX = iv[I].x + U / iv[I].fx * (1 - t/2.);\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"\t\tX = iv[I].x + (U * iv[I].Tfx * iv[I].Tfx) / (1.-iv[I].Tfx*iv[I].dTfx*U);\n");
    break;
  } /* end switch */

  /* accept or reject */
  fprintf(out,"\t\tV = uniform();\n");

  /* squeeze acceptance */
  fprintf(out,"\t\tif (V <= iv[I].sq) return X;\n");

  /* evaluate hat at X:
     get uniform random number between 0 and hat(X) */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"\t\tV *= iv[I].fx * exp(iv[I].dTfx * (X - iv[I].x));\n");
    break;
  case TDR_VAR_T_SQRT:
     /* transformed hat at X */ 
    fprintf(out,"\t\tThx = iv[I].Tfx + iv[I].dTfx * (X - iv[I].x);\n");
    fprintf(out,"\t\tV /= Thx*Thx;\n");
    break;
  } /* end switch */

  /* main rejection */
  fprintf(out,"\t\tif (V <= %s(X)) return X;\n",pdf_name);

  /* end of loop */
  fprintf(out,"\t}\n");

  /* end of function */
  fprintf(out,"}\n");

  /* o.k. */
  return 1;

} /* end of _unur_acg_C_tdr_ps() */

/*---------------------------------------------------------------------------*/
