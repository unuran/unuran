/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_codegen_FORTRAN.c                                        *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Make generator code for TDR.                                         *
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

#include "codegen_source.h"
#include <methods/tdr_source.h>

/*---------------------------------------------------------------------------*/

int
_unur_acg_FORTRAN_tdr_ps( FILE *out, 
			  const UNUR_GEN *gen, 
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
#define print_data(i,x)  \
  do { \
    fprintf(out,((i)==0)?"     *   ":(((i)%2)?", ":",\n     *   ")); \
    _unur_acg_FORTRAN_print_double(out,(x)); \
  } while (0)


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
  _unur_acg_FORTRAN_print_section_rule(out);
  _unur_acg_FORTRAN_print_section_line(out,"Sampling from %.35s distribution.",gen->distr->name);
  _unur_acg_FORTRAN_print_section_line(out,"Method: TDR - PS (Transformed Density Rejection / prop. squeeze)");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_acg_FORTRAN_print_section_line(out,"        Transformation = log(x) ... c = 0");         break;
  case TDR_VAR_T_SQRT:
    _unur_acg_FORTRAN_print_section_line(out,"        Transformation = -1/sqrt(x)  ... c = -1/2"); break;
  case TDR_VAR_T_POW:
    _unur_acg_FORTRAN_print_section_line(out,"        Transformation = -x^c  ... c = %g",GEN.c_T); break;
  }
  _unur_acg_FORTRAN_print_section_line(out,"        hat / squeeze ratio = %g",GEN.Atotal / GEN.Asqueeze);
  _unur_acg_FORTRAN_print_section_rule(out);
  fprintf(out,"\n");

  /* sampling routine */
  fprintf(out,"      DOUBLE PRECISION FUNCTION %s()\n",rand_name);
  fprintf(out,"\n");

  /* constants */
  fprintf(out,"C\n");
  fprintf(out,"C     data\n");
  fprintf(out,"C\n");

  fprintf(out,"      IMPLICIT DOUBLE PRECISION (A-H,O-Z)\n");
  fprintf(out,"      INTEGER GSIZE, GUIDE\n");
  fprintf(out,"      PARAMETER (gsize=%d)\n",GEN.guide_size);
  fprintf(out,"      PARAMETER (Atotal=");
  _unur_acg_FORTRAN_print_double(out,GEN.Atotal);
  fprintf(out,")\n");
  fprintf(out,"      DIMENSION GUIDE(0:gsize-1)\n");
  fprintf(out,"      DIMENSION x(0:%d)\n",GEN.n_ivs-1);
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"      DIMENSION fx(0:%d)\n",GEN.n_ivs-1);
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"      DIMENSION Tfx(0:%d)\n",GEN.n_ivs-1);
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }
  fprintf(out,"      DIMENSION dTfx(0:%d)\n",GEN.n_ivs-1);
  fprintf(out,"      DIMENSION sq(0:%d)\n",GEN.n_ivs-1);
  fprintf(out,"      DIMENSION Acum(0:%d)\n",GEN.n_ivs-1);
  fprintf(out,"      DIMENSION Ahatr(0:%d)\n",GEN.n_ivs-1);
  fprintf(out,"      \n");

  /* guide table */
  fprintf(out,"      DATA guide/\n");
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
      fprintf(out,"     *   %d",j);
    else
      fprintf(out,(i%12)?",%d":",\n     *   %d",j);
  }      
  fprintf(out,"/\n");

  /* hat function */
  fprintf(out,"      DATA x/\n");
  for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
    print_data(i,iv->x);
  fprintf(out,"/\n");

  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"      DATA fx/\n");
    for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
      print_data(i,iv->fx);
    fprintf(out,"/\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"      DATA Tfx/\n");
    for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
      print_data(i,iv->Tfx);
    fprintf(out,"/\n");
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return 0;
  }

  fprintf(out,"      DATA dTfx/\n");
  for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
    print_data(i,iv->dTfx);
  fprintf(out,"/\n");

  fprintf(out,"      DATA sq/\n");
  for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
    print_data(i,iv->sq);
  fprintf(out,"/\n");

  fprintf(out,"      DATA Acum/\n");
  for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
    print_data(i,iv->Acum);
  fprintf(out,"/\n");

  fprintf(out,"      DATA Ahatr/\n");
  for (i=0, iv=GEN.iv; iv->next!=NULL; i++, iv=iv->next)
    print_data(i,iv->Ahatr);
  fprintf(out,"/\n\n");

  fprintf(out,"C\n");
  fprintf(out,"C     code\n");
  fprintf(out,"C\n\n");

  /* main loop */
  fprintf(out,"1     CONTINUE\n");

  /* sample from U(0,1) */
  fprintf(out,"         U = urand()\n");
  fprintf(out,"         W = U * Atotal\n");

  /* look up in guide table and search for interval */
  fprintf(out,"         DO 3 I = guide(INT(U*gsize)), %d\n",GEN.n_ivs-1);
  fprintf(out,"            IF (Acum(I).GE.W) GOTO 5\n");
  fprintf(out,"3        CONTINUE\n");

  /* reuse of uniform random number */
  fprintf(out,"5        W = W - Acum(I) + Ahatr(I)\n");
  /* result: U in (-A_hatl, A_hatr) */

  /* generate from hat distribution */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"         t = dTfx(I) * W / fx(I)\n");
    fprintf(out,"         IF (DABS(t) > 1.d-8) THEN\n");
    fprintf(out,"            %s = x(I) + DLOG(t+1.d0) * W / (fx(I) * t) \n",rand_name);
    fprintf(out,"         ELSE\n");
    fprintf(out,"            %s = x(I) + W / (fx(I) * (1.d0 - t/2.d0)) \n",rand_name);
    fprintf(out,"         ENDIF\n");
    break;
  case TDR_VAR_T_SQRT:
    fprintf(out,"         %s = x(I)+(W*Tfx(I)*Tfx(I))/(1.d0-Tfx(I)*dTfx(I)*W)\n",rand_name);
    break;
  } /* end switch */

  /* accept or reject */
  fprintf(out,"         V = urand()\n");

  /* squeeze acceptance */
  fprintf(out,"         IF (V .LE. sq(I)) RETURN\n");

  /* evaluate hat at X:
     get uniform random number between 0 and hat(X) */
  switch (gen->variant & TDR_VARMASK_T) {
  case TDR_VAR_T_LOG:
    fprintf(out,"         V = V * fx(I) * DEXP(dTfx(I) * (%s - x(I)))\n",rand_name);
    break;
  case TDR_VAR_T_SQRT:
     /* transformed hat at X */ 
    fprintf(out,"         Thx = Tfx(I) + dTfx(I) * (%s - x(I))\n",rand_name);
    fprintf(out,"         V = V / (Thx*Thx)\n");
    break;
  } /* end switch */

  /* main rejection */
  fprintf(out,"      IF (V .GT. %s(%s)) GOTO 1\n",pdf_name,rand_name);

  /* end of function */
  fprintf(out,"\n");
  fprintf(out,"      END\n");

  /* o.k. */
  return 1;

#undef print_data

} /* end of _unur_acg_FORTRAN_tdr_ps() */

/*---------------------------------------------------------------------------*/
