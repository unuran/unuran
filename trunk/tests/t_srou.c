/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Tests                                                                    *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/
#ifdef T_SROU
/*---------------------------------------------------------------------------*/

static const char method_name[] = "srou";

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static int test_ok = TRUE;        /* all tests ok (boolean)                  */
static FILE *TESTLOG = NULL;      /* test log file                           */
static FILE *UNURANLOG = NULL;    /* unuran log file                         */

/*---------------------------------------------------------------------------*/

int main()
{ 
  char filename[128];

  /* open test log file */
  sprintf(filename,"%s_test.log",method_name);
  TESTLOG = fopen(filename,"w");
  if (TESTLOG == NULL) exit (-1);

  /* set output stream for unuran messages */
  sprintf(filename,"%s_unuran.log",method_name);
  UNURANLOG = fopen(filename,"w");
  if (UNURANLOG == NULL) exit (-1);
  unur_set_stream( UNURANLOG );
  
  /* close log files and exit */
  fclose(TESTLOG);
  fclose(UNURANLOG);
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#else
/*---------------------------------------------------------------------------*/
int main() { exit(77); } /* ignore test */
/*---------------------------------------------------------------------------*/
#endif  /* T_SROU */
/*---------------------------------------------------------------------------*/


#if 0
UNUR_PAR *unur_srou_new( UNUR_DISTR *distribution );
int unur_srou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
int unur_srou_set_pdfatmode( UNUR_PAR *parameters, double fmode );
int unur_srou_set_verify( UNUR_PAR *parameters, int verify );
int unur_srou_set_usesqueeze( UNUR_PAR *parameters, int usesqueeze );
int unur_srou_set_usemirror( UNUR_PAR *parameters, int usemirror );
int unur_srou_chg_pdfparams( UNUR_GEN *generator, double *params, int n_params );
int unur_srou_chg_domain( UNUR_GEN *generator, double left, double right );
int unur_srou_chg_mode( UNUR_GEN *generator, double mode );
int unur_srou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
int unur_srou_chg_pdfatmode( UNUR_GEN *generator, double fmode );
int unur_srou_chg_pdfarea( UNUR_GEN *generator, double area );
#endif

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/
