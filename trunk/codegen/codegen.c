/*---------------------------------------------------------------------------*/

#include <stdarg.h>
#include "codegen_source.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_acg( struct unur_gen *gen, FILE *out, const char *distr_name )
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


  switch (gen->method) {
  case UNUR_METH_TDR:
    if (! _unur_acg_C_PDF(&(gen->distr),out,pdf_name)) {
      return_code = 0;
    }
    else
      return_code = _unur_tdr_ps_codegen( gen, out, rand_name, pdf_name );
    break;
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return_code = 0;
  }

  /* clear memory */
  free(pdf_name);
  free(rand_name);

  /* end */
  return return_code;

} /* end of unur_acg() */

/*---------------------------------------------------------------------------*/

int
unur_xxx_default_urng( FILE *out, const char *name)
{

  return 1;
} /* end of unur_xxx_default_urng() */

/*---------------------------------------------------------------------------*/

void
_unur_acg_print_sectionheader( FILE *out, int n_lines, ... )
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

} /* end of _unur_acg_print_sectionheader() */

/*---------------------------------------------------------------------------*/
