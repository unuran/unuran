/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_codegen_UNURAN.c                                         *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Make generator code for TDR.                                         *
 *      (UNURAN Version)                                                     *
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
_unur_acg_UNURAN_tdr_ps( FILE *out, 
			 struct unur_gen *gen, 
			 const char *rand_name,
			 int n_cpoints ) 
     /*----------------------------------------------------------------------*/
     /* code generator for method TDR variant PS (proportional squeeze)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out       ... output stream                                        */
     /*   gen       ... pointer to generator object                          */
     /*   rand_name ... name of sampling routine                             */
     /*   n_cpoints ... number of construction points                        */
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
  _unur_acg_C_print_section_rule(out);
  _unur_acg_C_print_section_line(out,"Sampling from %.35s distribution.",gen->distr.name);
  _unur_acg_C_print_section_line(out,"Method: TDR - PS (Transformed Density Rejection / prop. squeeze)");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_acg_C_print_section_line(out,"        Transformation = log(x)  ... c = 0");        break;
  case TDR_VAR_T_SQRT:
    _unur_acg_C_print_section_line(out,"        Transformation = -1/sqrt(x)  ... c = -1/2"); break;
  case TDR_VAR_T_POW:
    _unur_acg_C_print_section_line(out,"        Transformation = -x^c  ... c = %g",GEN.c_T); break;
  }
  _unur_acg_C_print_section_line(out,"        hat / squeeze ratio = %g",GEN.Atotal / GEN.Asqueeze);
  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  /* include header files */
  fprintf(out,"#include <unuran.h>\n\n");

  /* sampling routine */
  fprintf(out,"double %s (void)\n{\n",rand_name);

  /* declare variables */
  fprintf(out,"\tUNUR_DISTR *distr;\n");           /* distribution object */
  fprintf(out,"\tUNUR_PAR *par;\n");               /* parameter object */
  fprintf(out,"\tstatic UNUR_GEN *gen = NULL;\n"); /* generator object */
  fprintf(out,"\n");

  /* set up */

  /* copy PDF parameters */
  fprintf(out,"\tif (gen == NULL) {\n\n"); 
  n_fpar = unur_distr_cont_get_pdfparams(unur_get_distr(gen),&fpar);
  if (n_fpar == 0)
    fprintf(out,"\t\tdouble *fpar = NULL;\n");
  else {
    fprintf(out,"\t\tdouble fpar[] = {\n");
    for (i=0; i<n_fpar-1; i++)
      fprintf(out,"\t\t\t%.20e,\n",fpar[i]);
    fprintf(out,"\t\t\t%.20e };\n",fpar[n_fpar-1]);
  }

  /* make distribution object */
  fprintf(out,"\t\tdistr = unur_distr_%s(fpar,%d);\n",
	  unur_distr_get_name(unur_get_distr(gen)),n_fpar );

  /* set domain (if changed) */
  if( !(gen->distr.set & UNUR_DISTR_SET_STDDOMAIN) ) {
    fprintf(out,"\t\tunur_distr_cont_set_domain(distr, ");
    if (gen->distr.data.cont.domain[0] <= -UNUR_INFINITY) 
      fprintf(out,"-UNUR_INFINITY, ");
    else
      fprintf(out,"%.20e, ", gen->distr.data.cont.domain[0] );
    if (gen->distr.data.cont.domain[1] >= UNUR_INFINITY) 
      fprintf(out,"UNUR_INFINITY );\n");
    else
      fprintf(out,"%.20e );\n", gen->distr.data.cont.domain[1] );
  }
  fprintf(out,"\n");

  /* make parameter object (using TDR) */
  fprintf(out,"\t\tpar = unur_tdr_new(distr);\n");

  /* set necessary parameters for transformation method */
  fprintf(out,"\t\tunur_tdr_set_variant_ps(par);\n");
  fprintf(out,"\t\tunur_tdr_set_cpoints(par, %d, NULL);\n",n_cpoints);
  fprintf(out,"\t\tunur_tdr_set_c(par, %g);\n",gen->data.tdr.c_T);

  fprintf(out,"\n");
  _unur_acg_C_print_section_rule(out);
  _unur_acg_C_print_section_line(out, "Uniform (pseudo-)random number generator");
  _unur_acg_C_print_section_line(out,"");
  _unur_acg_C_print_section_line(out,"The UNURAN default uniform random number generator is used.");
  _unur_acg_C_print_section_line(out,"It can be replaced by a urng of your choice using the");
  _unur_acg_C_print_section_line(out,"unur_set_urng() call (see UNURAN user manual for details).");
  _unur_acg_C_print_section_line(out,"Eg, suppose `my_urng' is a pointer to your urng then uncomment");
  _unur_acg_C_print_section_line(out,"");
  _unur_acg_C_print_section_line(out,"             unur_set_urng(par,my_urng);");
  _unur_acg_C_print_section_line(out,"");
  _unur_acg_C_print_section_line(out,"");
  _unur_acg_C_print_section_rule(out);
  fprintf(out,"\n");

  /* initialize generator */
  fprintf(out,"\t\tgen = unur_init(par);\n");

  /* emergency exit */
  fprintf(out,"\t\tif (gen == NULL) {\n");
  fprintf(out,"\t\t\tfprintf(stderr, \"Cannot create generator object.\\n\");\n");
  fprintf(out,"\t\t\texit (EXIT_FAILURE);\n");
  fprintf(out,"\t\t}\n\n");

  /* destroy distribution object */
  fprintf(out,"\t\tunur_distr_free(distr);\n");

  fprintf(out,"\t}\n");  
  fprintf(out,"\n");

  /* sample */
  fprintf(out,"\treturn unur_sample_cont(gen);\n");

  /* end of function */
  fprintf(out,"}\n");

  /* o.k. */
  return 1;

} /* end of _unur_acg_UNURAN_tdr_ps() */

/*---------------------------------------------------------------------------*/
