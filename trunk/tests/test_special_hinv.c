
#include <unuran.h>

/*---------------------------------------------------------------------------*/
/* global variables                                                          */

static FILE *TESTLOG;               /* test log file                         */
static FILE *UNURANLOG;             /* unuran log file                       */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Example of  a distribution with f(x)=0 at x=-0.25 and x=0.75              */
/*                                                                           */
/*          /  (1 +sin(2 Pi x))/2  if |x| <= 1                               */
/*  f(x) = <                                                                 */
/*          \  0        otherwise                                            */
/*                                                                           */
double sinpdf( double x, const UNUR_DISTR *distr )
{
  if (fabs(x) > 1.)
    return 0.;
  else
    return (0.5*(1.+sin((2.*M_PI)*x)));
} /* end of sinpdf() */

/* The derivative of the PDF of our distribution: */
double sindpdf( double x, const UNUR_DISTR *distr )
{
  if (fabs(x) > 1.)
    return 0.;
  else
    return (M_PI*cos((2.*M_PI)*x));
} /* end of sindpdf() */

/* The CDF of our distribution: */
double sincdf( double x, const UNUR_DISTR *distr )
{
  if (x < -1.)
    return 0.;
  else if(x<=1.)
    return (1.+2.*M_PI*(1+x)-cos((2.*M_PI)*x))/(4.*M_PI);
  else return 1.;
} /* end of sincdf() */

/*---------------------------------------------------------------------------*/

/* Example of  a distribution with f(x)=0 in (-0.5,0)                        */
/*                                                                           */
/*          /  Max(sin(2 Pi x)),0)Pi/2  if -1 < x <0.5                       */
/*  f(x) = <                                                                 */
/*          \  0        otherwise                                            */
/*                                                                           */
double sin0pdf( double x, const UNUR_DISTR *distr )
{
  if (x < -1.) return 0.;
  if (x <= -0.5) return sin((2.*M_PI)*x)*0.5*M_PI;
  if (x < 0.) return 0.;
  if (x <= 0.5) return sin((2.*M_PI)*x)*0.5*M_PI;
  return 0.;
} /* end of sin0pdf() */

/* The derivative of the PDF of our distribution: */
double sin0dpdf( double x, const UNUR_DISTR *distr )
{
  if (x < -1.) return 0.;
  if (x <= -0.5) return cos((2.*M_PI)*x)*M_PI*M_PI;
  if (x < 0.) return 0.;
  if (x <= 0.5) return cos((2.*M_PI)*x)*M_PI*M_PI;
  return 0.;
} /* end of sin0dpdf() */

/* The CDF of our distribution: */
double sin0cdf( double x, const UNUR_DISTR *distr )
{
  if (x < -1.) return 0.;
  if(x<= -0.5) return 0.25*(1-cos((2.*M_PI)*x));
  if (x < 0.) return 0.5;
  if (x <= 0.5) return 0.75-0.25*cos((2.*M_PI)*x);
  return 1.;
} /* end of sin0cdf() */

/*****************************************************************************/

int
hinv_error_gen_experiment( UNUR_GEN *gen,
			   int samplesize,    /* n for error experiment      */
			   int output,        /* 0-no,1- little, 2-more      */
			   int order,
			   double eu )
     /* returns 0 if MAE > eu, 1 if maxerror>eu  otherwise 2                 */
{   
  double maxerror, MAE;
  int i, nfpar;
  const double *fpar;     
  const UNUR_DISTR *distr = unur_get_distr(gen);

  unur_hinv_estimate_error(gen, samplesize, &maxerror, &MAE);

  if((output==1 && MAE > eu)||output==2){
    fprintf(TESTLOG,"\n eu %.2e order %d ;%s Distribution: ",
	   eu, order, unur_distr_get_name(distr));
    nfpar = unur_distr_cont_get_pdfparams(distr,&fpar);
    if(nfpar){ 
      fprintf(TESTLOG,"With parameters  ");
      for(i=0;i<nfpar;i++) fprintf(TESTLOG," | %f",fpar[i]);
      fprintf(TESTLOG," |");
    }
    fprintf(TESTLOG,"\n");
    fprintf(TESTLOG,"Total number of points: %d\n",
	 unur_hinv_get_n_intervals(gen));
    if(maxerror > eu)
      fprintf(TESTLOG,"UNUTEST Warning: Precision Problem, maxerror bigger than u_resolution\n");
    fprintf(TESTLOG,"\n samplesize %d maxerror %e MAE %e\n",samplesize,maxerror,MAE);
  } 

  if(MAE > eu){
    printf("(!+)");
      return 0;
  }else if(maxerror > eu){
    printf("(+?)");
      return 1;
  }else{
    printf("+");
      return 2;
  }
} /* end of hinv_error_gen_experiment() */

/*****************************************************************************/

int
hinv_error_experiment( UNUR_DISTR *distr,
		       double *cpoints,   /* pointer to c-point array        */
		       int ncp,           /* number of c-points              */
		       double eu,         /* epsilon-u                       */
		       int order,         /* order of polynomial             */
		       int samplesize,    /* n for error experiment          */
		       int output)        /* 0-no,1- little, 2-more          */
     /* returns 0 if MAE > eu, 1 if maxerror>eu  otherwise 2 */

{ 
  int i, nfpar;
  const double *fpar;     
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  nfpar = unur_distr_cont_get_pdfparams(distr,&fpar);
  par = unur_hinv_new(distr);
  unur_hinv_set_order(par,order);
  if(cpoints != 0 && ncp>0) unur_hinv_set_cpoints(par, cpoints, ncp);
  unur_hinv_set_u_resolution(par,eu);
  unur_set_debug(par,1u);
  gen = unur_init(par);

  i = hinv_error_gen_experiment(gen, samplesize, output, order, eu);

  unur_free(gen);

  return i;
} /* end of hinv_error_experiment() */

/*****************************************************************************/

int
check_domain( UNUR_GEN   *gen,
	      int samplesize,
	      double left,
	      double right,
	      int output )
     /* returns 0 if samples outside of left right were generated, 
	1 if samples were too far away from border,
	otherwise 2
     */
{ 
  int i, retvalue;
  double x, min, max, cdfleft,cdfright,cdfmin,cdfmax;
  const UNUR_DISTR *distr = unur_get_distr(gen);

  min = 1.e99;
  max = -1.e99;
  for(i=0; i<samplesize; i++){
    x = unur_sample_cont(gen);
    if(x < min) min = x;
    if(x > max) max = x;
  }
  cdfmin = unur_distr_cont_eval_cdf(min,distr);
  cdfleft = unur_distr_cont_eval_cdf(left,distr); 
  cdfmax = unur_distr_cont_eval_cdf(max,distr);
  cdfright = unur_distr_cont_eval_cdf(right,distr); 
  
  if(min < left || max > right) retvalue = 0;
  if(cdfmin-cdfleft > 100./samplesize
     ||cdfright - cdfmax > 100./samplesize )
    retvalue = 1;
  else retvalue = 2;
  
  if(output==2 || (output==1 && retvalue<=1)){   
    if(retvalue<=1) fprintf(TESTLOG,"ERROR with domain:\n");
    fprintf(TESTLOG,"left %e min %e max %e right %e\n",left,min,max,right);
    fprintf(TESTLOG,"cdfleft %e cdfmin %e cdfmax %e cdfright %e\n",cdfleft,cdfmin,cdfmax,cdfright);
  }
  
  if(retvalue == 2) printf("+");
  if(retvalue == 1) printf("(+?)");
  if(retvalue == 0) printf("+");
  
  return retvalue;
} /* check_domain() */

/*****************************************************************************/

int
chg_domain_experiment( UNUR_GEN   *gen, 
		       double eu,
		       int order,
		       int samplesize,
		       int output,
		       double left,   /*border of userdefined original domain*/
		       double right ) /*border of userdefined original domain*/
     /* returns 0 if MAE > eu or a domainerror occured,
	1 if maxerror>eu occured,
	2 otherwise 
     */
{ 
  int i,k,mini=9999;
  double newleft, newright;

  for(k=1; k<=10;k++){
    newleft = left + k*0.05*(right-left);
    newright = right - k*0.04*(right-left);
    unur_hinv_chg_truncated(gen, newleft, newright);
    i = hinv_error_gen_experiment(gen, samplesize/10, output, order, eu);
    if(i < mini) mini = i;
    i = check_domain(gen,samplesize,newleft,newright,output);
    if(i < mini) mini = i;
    unur_hinv_chg_truncated(gen, left, right);
    i = hinv_error_gen_experiment(gen, samplesize/10, output, order, eu);
    if(i < mini) mini = i;
    i = check_domain(gen,samplesize,left,right,output);
    if(i < mini) mini = i;
  }
  
  return mini;
} /* end of chg_domain_experiment() */

/*****************************************************************************/

int main()
{
  int    samplesize=10000,order=3, errorsum = 0;
  double fpar[4], eu=1.e-8,cpoints[10];     

  UNUR_DISTR *distr;
  UNUR_PAR   *par;      /* parameter object                      */
  UNUR_GEN   *gen;      /* generator object                      */

  /* open log file for unuran and set output stream for unuran messages */
  UNURANLOG = fopen( "test_special_hinv_unuran.log","w" );
  unur_set_stream( UNURANLOG );

  /* open log file for testing */
  TESTLOG = fopen( "test_special_hinv_testlog.log","w" );

  /* write header into log file */
  fprintf(TESTLOG,"\nUNURAN - Universal Non-Uniform RANdom number generator\n\n");
  fprintf(TESTLOG,"\n====================================================\n\n");
   
  for(eu=1.e-6;eu>1.e-13;eu*=0.01)
    for(order=(eu<1.e-9? 3: 1); order<=5; order+=2)
      {   
	distr = unur_distr_normal(NULL,0);
	hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2);
	unur_distr_free(distr);
	
	distr = unur_distr_cauchy(NULL,0);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	distr = unur_distr_exponential(NULL,0);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	fpar[0]= 5.; 
	distr = unur_distr_gamma(fpar,1);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	fpar[0]= 0.5; 
	distr = unur_distr_gamma(fpar,1);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	fpar[0]=2.; fpar[1]=2.;fpar[2]= 0.; fpar[3]=1.;
	distr = unur_distr_beta(fpar,4);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	fpar[0]=0.3; fpar[1]=3.;fpar[2]= 0.; fpar[3]=1.;
	distr = unur_distr_beta(fpar,4);
	if(!hinv_error_experiment(distr,NULL,0,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);

	distr = unur_distr_cont_new();
	unur_distr_set_name(distr,"sin-example");
	unur_distr_cont_set_cdf( distr, sincdf );
	unur_distr_cont_set_pdf( distr, sinpdf );
	unur_distr_cont_set_dpdf( distr, sindpdf );
	unur_distr_cont_set_domain( distr, -1., 1. );
	cpoints[0]= -0.75;
	cpoints[1]= -0.25;
	cpoints[2]=  0.25;
	cpoints[3]=  0.75;
	if(!hinv_error_experiment(distr,cpoints,4,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
	
	distr = unur_distr_cont_new();
	unur_distr_set_name(distr,"sin-0-example");
	unur_distr_cont_set_cdf( distr, sin0cdf );
	unur_distr_cont_set_pdf( distr, sin0pdf );
	unur_distr_cont_set_dpdf( distr, sin0dpdf );
	unur_distr_cont_set_domain( distr, -1., 0.5 );
	cpoints[0]= -0.75;
	cpoints[1]= -0.5;
	cpoints[2]=  0.;
	cpoints[3]=  0.25;
	if(!hinv_error_experiment(distr,cpoints,4,eu,order,samplesize,2)) errorsum++;
	unur_distr_free(distr);
      }
 
  /* extra test for changeing the domain of the generator object */

  printf("test for checking unur_chg_domain(): normal"); 
  
  distr = unur_distr_normal(NULL,0);
  par = unur_hinv_new(distr);
  
  eu=1.e-10; order=3; 
  
  unur_hinv_set_order(par,order);
  unur_hinv_set_u_resolution(par,eu);
  gen = unur_init(par);
  
  if(!chg_domain_experiment(gen,eu,order,samplesize,2, -6., 6.)) errorsum++;

  unur_free(gen);
  unur_distr_free(distr);

  /****/

  printf("sinus1");
  distr = unur_distr_cont_new();
  unur_distr_set_name(distr,"sin-example");
  unur_distr_cont_set_cdf( distr, sincdf );
  unur_distr_cont_set_pdf( distr, sinpdf );
  unur_distr_cont_set_dpdf( distr, sindpdf );
  unur_distr_cont_set_domain( distr, -1., 1. );
  cpoints[0]= -0.75;
  cpoints[1]= -0.25;
  cpoints[2]=  0.25;
  cpoints[3]=  0.75;
  par = unur_hinv_new(distr);
  unur_hinv_set_order(par,order);
  unur_hinv_set_cpoints(par, cpoints, 4);
  unur_hinv_set_u_resolution(par,eu);
  gen = unur_init(par);

  if(!chg_domain_experiment(gen,eu,order,samplesize,2, -1., 1.)) errorsum++;

  unur_distr_free(distr);
  unur_free(gen);

  printf("sinus2");
  distr = unur_distr_cont_new();
  unur_distr_set_name(distr,"sin-0-example");
  unur_distr_cont_set_cdf( distr, sin0cdf );
  unur_distr_cont_set_pdf( distr, sin0pdf );
  unur_distr_cont_set_dpdf( distr, sin0dpdf );
  unur_distr_cont_set_domain( distr, -1., 0.5 );
  cpoints[0]= -0.75;
  cpoints[1]= -0.5;
  cpoints[2]=  0.;
  cpoints[3]=  0.25;
  par = unur_hinv_new(distr);
  unur_hinv_set_order(par,order);
  unur_hinv_set_cpoints(par, cpoints, 4);
  unur_hinv_set_u_resolution(par,eu);
  gen = unur_init(par);

  if(!chg_domain_experiment(gen,eu,order,samplesize,2, -1., 0.5)) errorsum++;

  unur_distr_free(distr);
  unur_free(gen);

  /* test finished */
  printf("\n");  fflush(stdout);

  /* close log files */
  fprintf(TESTLOG,"\n====================================================\n\n");
  if (errorsum == 0)
    fprintf(TESTLOG,"All tests PASSED.\n");
  else
    fprintf(TESTLOG,"Test(s) FAILED.\n");
  
  fclose(UNURANLOG);
  fclose(TESTLOG);
  
  exit ((errorsum == 0) ? EXIT_SUCCESS : EXIT_FAILURE);

} /* end of main() */

/* ------------------------------------------------------------- */
