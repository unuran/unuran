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

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

int test_ok = TRUE;               /* all tests ok (boolean)                  */
int test_ok_local = TRUE;         /* running test ok (boolean)               */
FILE *TESTLOG = NULL;             /* test log file                           */

static FILE *UNURANLOG = NULL;    /* unuran log file                         */

/*---------------------------------------------------------------------------*/

void test_srou_new( void );

/*---------------------------------------------------------------------------*/

int main()
{ 
  char filename[128];
  /* open test log file */
  sprintf(filename,"%s_test.log",__FILE__);
  TESTLOG = fopen(filename,"w");
  if (TESTLOG == NULL) exit (-1);

  /* set output stream for unuran messages */
  sprintf(filename,"%s_unuran.log",__FILE__);
  UNURANLOG = fopen(filename,"w");
  if (UNURANLOG == NULL) exit (-1);
  unur_set_stream( UNURANLOG );

  /* start test */
  printf("%s: ",__FILE__);

  /* run tests */
  test_srou_new();

  /* test finished */
  printf("\n");

  /* close log files and exit */
  fclose(TESTLOG);
  fclose(UNURANLOG);
  exit( (test_ok) ? 0 : -1 );
}

/*---------------------------------------------------------------------------*/

void test_srou_new( void )
{
  // UNUR_PAR *par;

  /* start test */
  printf("[new ");
  fprintf(TESTLOG,"[new]\n");

  test_ok_local = TRUE;

  /* check error handling */
  check_expected_NULL( (unur_srou_new(NULL)) );
     
  /* test finished */
  test_ok &= test_ok_local;
  (test_ok_local) ? printf("... ok] ") : printf("... failed] ");
}

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
