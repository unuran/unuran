/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
     $Id$ 
 *****************************************************************************
 *                                                                           *
 *  Compare CDF, PDF and derivatives of PDF of varios distributions          *
 *  with Mathematica output.                                                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

/*  #define DEBUG 1 */

/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unuran.h>
#include <config.h>

#ifdef WITH_DMALLOC
#  include <dmalloc.h>
#endif

/*---------------------------------------------------------------------------*/

int test_cdf_pdf( UNUR_DISTR *distr, char *datafile, double max_diff );
/* compare CDF, PDF and dPDF for continuous univariate distributions         */

/*---------------------------------------------------------------------------*/
