#include <source_unuran.h>
#include "PDFgen_source.h"


/*---------------------------------------------------------------------------*/

int
unurgen (struct unur_gen *gen, FILE *out, const char *distr_name)
{
  _unur_check_NULL("unurgen", gen, 0 );

  switch (gen->method) {
  case UNUR_METH_TDR:
    return _unur_tdr_ps_codegen( gen, out, distr_name );
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return 0;
  }

} /* end of unurgen() */

/*---------------------------------------------------------------------------*/
