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

#define DEBUG 1

/*---------------------------------------------------------------------------*/

#include "test_with_Mathematica.h"

/* maximal difference allowed (relative to maximal value) */
#define MAX_REL_DIFF 1e-10

/*---------------------------------------------------------------------------*/

int test_cdf_pdf( FILE *LOG, UNUR_DISTR *distr, char *datafile )
     /*----------------------------------------------------------------------*/
     /* test CDF, PDF and derivative of PDF by comparing to data             */
     /* read from file created by Mathematica.                               */
     /*                                                                      */
     /* the ordering in this file is:                                        */
     /*   n   par[1] par[2] ... par[n]   x   CDF   PDF   dPDF                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   LOG      ... file handle for log file                              */
     /*   distr    ... pointer to distribution object                        */
     /*   datafile ... name of file where data are stored                    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if test was successful                                       */
     /*   0 ... failed (difference too large)                                */
     /*  -1 ... could not run test                                           */
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
  double x;                     /* double argument */
  int k = 0;                    /* integer argument */
  double CDF_e, PDF_e, dPDF_e;  /* expected values for CDF, PDF and derivative of PDF at x */
  double CDF_o, PDF_o, dPDF_o;  /* observed values for CDF, PDF and derivative of PDF at x */
  double CDF_d, PDF_d, dPDF_d;  /* differences between observed and expected values */

  double CDF_md, PDF_md, dPDF_md;  /* maximal difference (absolute) */
  double CDF_me, PDF_me, dPDF_me;  /* maximal expected absolute values */

  int have_CDF, have_PDF, have_dPDF, have_upd_pdfarea;

  int is_DISCR;                 /* 1 if discrete distribution, 0 otherwise */
  int n_failed = 0;             /* number of failed tests */

  const char *dname;            /* name of distribution */

  int i;

  /* discrete distribution ? */
  is_DISCR = unur_distr_is_discr( distr );

  /* name of distribution */
  dname = unur_distr_get_name(distr);

  /* initialize */
  CDF_md = PDF_md = dPDF_md = 0.;
  CDF_me = PDF_me = dPDF_me = 0.;

  /* check existence of CDF, PDF and dPDF */
  if (is_DISCR) {
    have_CDF  = (unur_distr_discr_get_cdf(distr)) ? 1 : 0;
    have_PDF  = (unur_distr_discr_get_pmf(distr)) ? 1 : 0;
    have_dPDF = 0;
  }
  else { /* is_CONT */
    have_CDF  = (unur_distr_cont_get_cdf(distr))  ? 1 : 0;
    have_PDF  = (unur_distr_cont_get_pdf(distr))  ? 1 : 0;
    have_dPDF = (unur_distr_cont_get_dpdf(distr)) ? 1 : 0;
  }

  /* check whether unur_distr_cont_upd_pdfarea() works */
  if (is_DISCR)
    have_upd_pdfarea = unur_distr_discr_upd_pmfsum(distr); 
  else  /* is_CONT */
    have_upd_pdfarea = unur_distr_cont_upd_pdfarea(distr); 

  if (!have_upd_pdfarea) {
    /* if we cannot update the area below the PDF, then the
       given PDF is not a "real" PDF, i.e. not normalized */
    have_PDF = 0;
    have_dPDF = 0;
  }

  /* print info into log file */
  fprintf(LOG,"%s: %s ...\n", dname, datafile);

  if (!have_CDF)
    fprintf(LOG,"%s: no CDF!\n", dname);
  if (!have_PDF)
    fprintf(LOG,"%s: no PDF!\n", dname);
  if (!have_dPDF)
    fprintf(LOG,"%s: no dPDF!\n", dname);

  /* open file for reading */
  if ( (fp = fopen(datafile,"r")) == NULL ) {
    printf("%s: ERROR: could not open file %s \n", dname,datafile);
    fprintf(LOG,"%s: ERROR: could not open file %s \n", dname,datafile);
    return -1;
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
      printf("%s: ERROR: invalid number of parameters for distribution: %d \n", dname,n_fparams);
      fprintf(LOG,"%s: ERROR: invalid number of parameters for distribution: %d \n", dname,n_fparams);
      return -1;
    }

    /* read parameters */
    for (i=0; i<n_fparams; i++) {
      fparams[i] = strtod( ptr_buffer, &ptr_buffer );
    }

    /* read argument */
    x = strtod( ptr_buffer, &ptr_buffer );
    if (is_DISCR) 
      /* integer required, round */
      k = (int)(x+0.5);

    /* read CDF */
    CDF_e = strtod( ptr_buffer, &ptr_buffer );

    /* read PDF */
    PDF_e = strtod( ptr_buffer, &ptr_buffer );

    /* read dPDF */
    dPDF_e = strtod( ptr_buffer, &ptr_buffer );

#ifdef DEBUG
    /* print data */
    fprintf(LOG, "%s: n_fparams = %d\n", dname, n_fparams);

    /* print parameters */
    for (i=0; i<n_fparams; i++)
      fprintf(LOG, "%s: \t%d: %g\n", dname, i, fparams[i]);

    /* print argument x (or k) */
    if (is_DISCR)
      fprintf(LOG, "%s: expected k = %d:\t", dname, k);
    else  /* is_CONT */
      fprintf(LOG, "%s: expected x = %g:\t", dname, x);

    /* print CDF, PDF and derivative of PDF at x */
    fprintf(LOG, "%g, %g, %g\n", CDF_e, PDF_e, dPDF_e);
#endif

    if (is_DISCR) {
      /* set parameters for distribution */
      unur_distr_discr_set_pmfparams(distr,fparams,n_fparams);
      if (have_upd_pdfarea) unur_distr_discr_upd_pmfsum(distr); 
      
      /* compute CDF, PDF and derivative of PDF at x */
      CDF_o = (have_CDF) ? unur_distr_discr_eval_cdf(x, distr) : 0.;
      PDF_o = (have_PDF) ? unur_distr_discr_eval_pmf(x, distr) : 0.;
      dPDF_o = 0.;
    }

    else { /* is_CONT */
      /* set parameters for distribution */
      unur_distr_cont_set_pdfparams(distr,fparams,n_fparams);
      if (have_upd_pdfarea) unur_distr_cont_upd_pdfarea(distr); 
      
      /* compute CDF, PDF and derivative of PDF at x */
      CDF_o = (have_CDF) ? unur_distr_cont_eval_cdf(x, distr) : 0.;
      PDF_o = (have_PDF) ? unur_distr_cont_eval_pdf(x, distr) : 0.;
      dPDF_o = (have_dPDF) ? unur_distr_cont_eval_dpdf(x, distr) : 0.;
    }

    /* compute differnces */
    CDF_d = CDF_o - CDF_e;
    PDF_d = PDF_o - PDF_e;
    dPDF_d = dPDF_o - dPDF_e;

#ifdef DEBUG
    /* print argument x (or k) */
    if (is_DISCR)
      fprintf(LOG, "%s: observed k = %d:\t", dname, k);
    else  /* is_CONT */
      fprintf(LOG, "%s: observed x = %g:\t", dname, x);

    /* print CDF, PDF and derivative of PDF at x */
    fprintf(LOG, "%g, %g, %g\n",CDF_o,PDF_o,dPDF_o);

    /* print differences */
    fprintf(LOG, "%s: diff     x = %g:\t", dname,x);
    fprintf(LOG, "%g, %g, %g\n",CDF_d, PDF_d, dPDF_d);
    fprintf(LOG, "%s:\n", dname);
#endif
    
    /* absolute values */
    CDF_d = fabs(CDF_d);
    PDF_d = fabs(PDF_d);
    dPDF_d = fabs(dPDF_d);

    CDF_e = fabs(CDF_e);
    PDF_e = fabs(PDF_e);
    dPDF_e = fabs(dPDF_e);

    /* maximal differences */
    if (CDF_d > CDF_md)   CDF_md = CDF_d;
    if (PDF_d > PDF_md)   PDF_md = PDF_d;
    if (dPDF_d > dPDF_md) dPDF_md = dPDF_d;

    /* maximal values */
    if (CDF_e > CDF_me)   CDF_me = CDF_e;
    if (PDF_e > PDF_me)   PDF_me = PDF_e;
    if (dPDF_e > dPDF_me) dPDF_me = dPDF_e;

  }

  /* close file handle */
  fclose(fp);

  /* print info on screen */
  fprintf(LOG, "%s: \tmaximal difference:\n", dname);

  if (have_CDF) {
    fprintf(LOG,"%s: \t\tCDF  = %g ... ", dname,CDF_md);
    if (CDF_md > MAX_REL_DIFF * CDF_me) {
      fprintf(LOG, "failed!!\n");
      ++n_failed;
    }
    else
      fprintf(LOG, "ok\n");
  }

  if (have_PDF) {
    fprintf(LOG,"%s: \t\tPDF  = %g ... ", dname,PDF_md);
    if (PDF_md > MAX_REL_DIFF * PDF_me) {
      fprintf(LOG, "failed!!\n");
      ++n_failed;
    }
    else
      fprintf(LOG, "ok\n");
  }

  if (have_dPDF) {
    fprintf(LOG,"%s: \t\tdPDF = %g ... ", dname,dPDF_md);
    if (dPDF_md > MAX_REL_DIFF * dPDF_me) {
      fprintf(LOG, "failed!!\n");
      ++n_failed;
    }
    else
      fprintf(LOG, "ok\n");
  }

  fprintf(LOG, "\n----------------------------------------\n\n");

  /* print result on screen */
  printf("%-17s ... ", dname );
  if (n_failed > 0) {
    printf("failed!!\n");
    return 0;
  }
  else {
    printf("ok\n");
    return 1;
  }

#undef BUFSIZE
#undef MAX_FPARAMS
} /* end of test_cdf_pdf() */

/*---------------------------------------------------------------------------*/
