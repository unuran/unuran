/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Common test routines                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include "t_unuran.h"

/*---------------------------------------------------------------------------*/

/* compare error code */
void do_compare_errno( int line, unsigned errno )
{
  if (unur_errno != errno) {
    test_ok_local = FALSE;
    fprintf(TESTLOG,"line %d: Wrong error code ... Failed\n",line);
  }
}

/* check for expected NULL pointer */
void do_check_expected_NULL( int line, void *ptr )
{
  if (ptr != NULL) { 
    test_ok_local = FALSE;
    fprintf(TESTLOG,"line %d: NULL pointer expected ... Failed\n",line);
  }
  do_compare_errno(line,UNUR_ERR_NULL);
}

/*---------------------------------------------------------------------------*/
