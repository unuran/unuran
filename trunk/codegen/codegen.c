/*---------------------------------------------------------------------------*/

#include "codegen_source.h"

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

int
unur_acg( struct unur_gen *gen, FILE *out, const char *distr_name )
{
  _unur_check_NULL("unur_acg", gen, 0 );

  switch (gen->method) {
  case UNUR_METH_TDR:
    return _unur_tdr_ps_codegen( gen, out, distr_name );
  default:
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot make generator code");
    return 0;
  }

} /* end of unur_acg() */

/*---------------------------------------------------------------------------*/

int
unur_xxx_default_urng( FILE *out, const char *name)
{

  return 1;
} /* end of unur_xxx_default_urng() */

/*---------------------------------------------------------------------------*/
