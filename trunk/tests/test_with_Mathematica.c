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

#include "test_with_Mathematica.h"

/*---------------------------------------------------------------------------*/

int test_cont_cdf_pdf( UNUR_DISTR *distr, char *datafile, double max_diff )
     /*----------------------------------------------------------------------*/
     /* test CDF, PDF and derivative of PDF by comparing to data             */
     /* read from file created by Mathematica.                               */
     /*                                                                      */
     /* the ordering in this file is:                                        */
     /*   n   par[1] par[2] ... par[n]   x   CDF   PDF   dPDF                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   datafile ... name of file where data are stored                    */
     /*   max_diff ... maximal allowed difference                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if test was successful                                       */
     /*   0 ... difference too large                                         */
     /*----------------------------------------------------------------------*/
{
#define BUFSIZE      1024       /* size of line buffer */
#define MAX_FPARAMS  10         /* maximal number of parameters for distribution */

  FILE *fp;                     /* file handle for input file */
  int lineno;                   /* line number in file */
  char buffer[BUFSIZE];         /* buffer for reading line */
  char *ptr_buffer;             /* pointer into buffer */

  int n_fparams;                /* number of parameters for distribution */
  double fparams[10];           /* array for parameters for distribution */
  double x;                     /* argument */
  double CDF_e, PDF_e, dPDF_e;  /* expected values for CDF, PDF and derivative of PDF at x */
  double CDF_o, PDF_o, dPDF_o;  /* observed values for CDF, PDF and derivative of PDF at x */
  double CDF_d, PDF_d, dPDF_d;  /* differences between observed and expected values */

  double CDF_md, PDF_md, dPDF_md;  /* maximal difference (absolute) */

  int i;

  /* initialize */
  CDF_md = PDF_md = dPDF_md = 0.;

  /* print info on screen */
  printf("%s ...\n",datafile);

  /* open file for reading */
  if ( (fp = fopen(datafile,"r")) == NULL ) {
    printf("ERROR: could not open file %s \n",datafile);
    return(-1.);
  }

  /* read file */
  for (lineno = 0;; ++lineno) {

    /* read next line */
    fgets(buffer, BUFSIZE-1, fp);
    if (feof(fp)) break;
    ptr_buffer = buffer;

    /* get number of parameters */
    n_fparams = strtol(buffer, &ptr_buffer, 10);
    if (n_fparams < 0 || n_fparams >= MAX_FPARAMS) {
      printf("ERROR: invalid number of parameters for distribution: %d \n",n_fparams);
      return(-1.);
    }

    /* read parameters */
    for (i=0; i<n_fparams; i++) {
      fparams[i] = strtod( ptr_buffer, &ptr_buffer );
    }

    /* read argument */
    x = strtod( ptr_buffer, &ptr_buffer );

    /* read CDF */
    CDF_e = strtod( ptr_buffer, &ptr_buffer );

    /* read PDF */
    PDF_e = strtod( ptr_buffer, &ptr_buffer );

    /* read dPDF */
    dPDF_e = strtod( ptr_buffer, &ptr_buffer );

#ifdef DEBUG
    /* print data */
    printf("n_fparams = %d\n", n_fparams);

    /* print parameters */
    for (i=0; i<n_fparams; i++)
      printf("\t%d: %g\n",i,fparams[i]);

    /* print argument x */
    printf("expected x = %g:\t",x);

    /* print CDF, PDF and derivative of PDF at x */
    printf("%g, %g, %g\n",CDF_e,PDF_e,dPDF_e);

#endif

    /* set parameters for distribution */
    unur_distr_cont_set_pdfparams(distr,fparams,n_fparams);
    unur_distr_cont_upd_pdfarea(distr); 

    /* compute CDF, PDF and derivative of PDF at x */
    CDF_o = unur_distr_cont_eval_cdf(x, distr);
    PDF_o = unur_distr_cont_eval_pdf(x, distr);
    dPDF_o = unur_distr_cont_eval_dpdf(x, distr);

    /* compute differnces */
    CDF_d = CDF_o - CDF_e;
    PDF_d = PDF_o - PDF_e;
    dPDF_d = dPDF_o - dPDF_e;

#ifdef DEBUG
    /* print argument x */
    printf("observed x = %g:\t",x);

    /* print CDF, PDF and derivative of PDF at x */
    printf("%g, %g, %g\n",CDF_o,PDF_o,dPDF_o);

    /* print differences */
    printf("diff     x = %g:\t",x);
    printf("%g, %g, %g\n\n",CDF_d, PDF_d, dPDF_d);
#endif
    
    /* compare differences */
    CDF_d = fabs(CDF_d);
    PDF_d = fabs(PDF_d);
    dPDF_d = fabs(dPDF_d);
    dPDF_e = fabs(dPDF_e);

    if (CDF_d > CDF_md)   CDF_md = CDF_d;
    if (PDF_d > PDF_md)   PDF_md = PDF_d;
    if (dPDF_d > dPDF_md) dPDF_md = dPDF_d;

  }

  /* close file handle */
  fclose(fp);

  /* print info on screen */
  printf("\tmaximal difference:\n");
  printf("\t\tCDF  = %g\n",CDF_md);
  printf("\t\tPDF  = %g\n",PDF_md);
  printf("\t\tdPDF = %g\n",dPDF_md);

  return 0;

#undef BUFSIZE
#undef MAX_FPARAMS
} 

/*---------------------------------------------------------------------------*/

