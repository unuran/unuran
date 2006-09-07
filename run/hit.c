#include <unuran.h>
#include <unur_struct.h>
#include <unur_source.h>
#include <unistd.h>
#include "../experiments/expect_source.h"
#include <experimental/hitrou.h>

double _unur_test_chi2test( double *prob, int *observed, int len, int classmin,
           int verbose, FILE *out );

/* global structure to hold all parameters given */
struct {
  int dim;
  long nsamples;
  long nloops;
  long skip;
  int shape;
  int type;
  int melt;
  double rel_param;
  int nhist_uv;
  int nhist_p;
  int u_planes;
  int verbose;
  int adaptive;
  double rho;
  double mixture;
  UNUR_GEN *gen;
  UNUR_DISTR *distr_normal;
} p;

#include "../experiments/expect.c"


void hit_set_dim(int *x)          {p.dim=x[0];}
void hit_set_nsamples(long *x)    {p.nsamples=x[0];}
void hit_set_nloops(long *x)      {p.nloops=x[0];}
void hit_set_skip(long *x)        {p.skip=x[0];}
void hit_set_shape(int *x)        {p.shape=x[0];}
void hit_set_type(int *x)         {p.type=x[0];}
void hit_set_melt(int *x)         {p.melt=x[0];}
void hit_set_rel_param(double *x) {p.rel_param=x[0];}
void hit_set_nhist_uv(int *x)     {p.nhist_uv=x[0];}
void hit_set_nhist_p(int *x)      {p.nhist_p=x[0];}
void hit_get_dim(int *x)          {x[0]=p.dim;}

double f_norm(double *x, void *par)
{
  long i;
  double f=0;
  
  for (i=0; i<p.dim; i++) {
    f += x[i]*x[i];
  }
  
  return sqrt(f);
}



double pdf_mix(const double *x, UNUR_DISTR *distr) {
  
  double f=0.;
  double *xmix;
  double *mean;
  int i;
  
  mean = _unur_vector_new(p.dim);
  
  xmix=_unur_vector_new(p.dim);
  memcpy(xmix, x, p.dim*sizeof(double));

    
  /* construct pdf-value as "PDF(x-p)+PDF(x+p)/2" */
  for (i=0; i<p.dim; i++) xmix[i]=x[i]-p.mixture;
  f+=unur_distr_cvec_eval_pdf(xmix, p.distr_normal);
  for (i=0; i<p.dim; i++) xmix[i]=x[i]+p.mixture;
  f+=unur_distr_cvec_eval_pdf(xmix, p.distr_normal);
  
  f/=2.;
  
  free(xmix);
  free(mean);
  
  return f;
}




void hit_separator() {
  printf("----------------------------------------------------------------\n");
}

void hit_showparameters() {
  printf("dim      = %d\n", p.dim );
  printf("nsamples = %ld\n",p.nsamples );
  printf("nloops   = %ld\n",p.nloops );
  printf("skip     = %ld\n",p.skip );
  if (p.shape==0) {
    printf("shape    = RoU shape\n" );
    printf("mix par  = %f\n", p.mixture );  
  }
  if (p.shape==1) {
    printf("shape    = testshape (rectangle) \n" );
    if (p.type==0) printf("type     = cube \n" );
    if (p.type==1) printf("type     = pancake \n" );
    if (p.type==2) printf("type     = cigar \n" );
    if (p.type==3) printf("type     = half \n" );
  }
  if (p.shape==2)  printf("shape    = testshape (simplex) \n" );
  if (p.shape==3)  printf("shape    = testshape (two simplices along the v direction) \n" );
  if (p.shape==4)  printf("shape    = testshape (two simplices along the u[0] direction) \n" );
  if (p.shape==5)  {
    printf("shape    = testshape (ellipsoid) \n" );
    if (p.type==0) printf("type     = sphere \n" );
    if (p.type==1) printf("type     = pancake \n" );
    if (p.type==2) printf("type     = cigar \n" );
    if (p.type==3) printf("type     = half \n" );  
  }
  printf("relparam = %f\n", p.rel_param );
  printf("nhist_uv = %d\n", p.nhist_uv );
  printf("nhist_p  = %d\n", p.nhist_p );
  if (p.melt==0) printf("melt     = no (normal run) \n" );
  if (p.melt==1) printf("melt     = yes (from center) \n" );
  if (p.melt==2) printf("melt     = yes (from edge) \n" );
  if (p.u_planes==0) printf("u_planes = no (infinite strip) \n" );
  if (p.u_planes==1) printf("u_planes = yes (full bounding rectangle) \n" );
  if (p.adaptive==0) printf("adaptive = no \n" );
  if (p.adaptive==1) printf("adaptive = yes \n" );
  printf("cov_rho = %f\n", p.rho );
  
  hit_separator();
}

/***************************************************************************************/

void hit_init() {
  /* setting default values */
  p.dim = 2 ;
  p.nsamples = 1;
  p.nloops = 1; /* number of sampling loops for 2. level test */
  p.skip = 0;
  p.shape = 0; /* RoU shape */
  p.type = 1; /* pancake */
  p.melt = 0; /* melting process or normal run */
  p.rel_param = 1.;
  p.nhist_uv=100; /* number of histogram bins (marginal, in each dimension) */
  p.nhist_p =10;  /* number of histogram bins for 2. level test of p-values */
  p.gen = NULL;
  p.distr_normal=NULL;
  p.u_planes = 1;
  p.verbose = 0;
  p.adaptive = 1;
  p.mixture=0.;
  p.rho=0.;
}

/***************************************************************************************/

void hit_sample_uv(double *uv) {
  double *x;

  x = _unur_vector_new(p.dim)  ; /* sample point */
  unur_sample_vec(p.gen, x);
  _unur_hitrou_get_point_current(p.gen, uv);

  free(x);
}

/***************************************************************************************/

void hit_sample_uv_array(long *nsamples, double *uv) {
  double *x;
  long i;

  x = _unur_vector_new(p.dim)  ; /* sample point */
  for (i=0; i<nsamples[0]; i++) {
    unur_sample_vec(p.gen, x);
    _unur_hitrou_get_point_current(p.gen, &uv[i*(p.dim+1)]);
  }

  free(x);
}


/***************************************************************************************/


void hit_run() {

//  UNUR_DISTR *covar_distr;
//  UNUR_PAR   *covar_par;
//  UNUR_GEN   *covar_gen;
  UNUR_DISTR *distr;
  UNUR_PAR   *par;

  double *mean, *covar, *mode;
  double *eigenvalues, *eigenvectors;
  double *x, *uv, *uv0, *uv1, *uv2;
  double *umin, *umax; /* bounding box boundaries [0,1]x[0,1]... */
  double vmax=1.; /* vmin=0. */
  double *relative_size; /* relative size of testrectangle  */

  int iuv;
  int ip;
  int ret;
  
  int *hist_uv, *hist_p;
  double *expected_uv;
  double pchi2;

  long loop;
  long pdfcount, sum_pdfcount=0;

  int i, j, d;

  double exact_integral ;
  double r0, r1;
  double integral;
  double integral_mean, integral_variance, integral_rmse;

  struct unur_funct_vgeneric faux;
  
  /* reserving storage for our arrays and matrices */
  mean = _unur_vector_new(p.dim);
  covar = _unur_vector_new(p.dim*p.dim);
  eigenvalues  = _unur_vector_new(p.dim);
  eigenvectors = _unur_vector_new(p.dim*p.dim);
  mode = _unur_vector_new(p.dim);
  
  x   = _unur_vector_new(p.dim)  ; /* sample point */
  uv  = _unur_vector_new(p.dim+1); /* (dim+1) for testrectangle !!! */
  uv0 = _unur_vector_new(p.dim+1); /* (dim+1) for testrectangle !!! */
  uv1 = _unur_vector_new(p.dim+1); /* (dim+1) for testrectangle !!! */
  uv2 = _unur_vector_new(p.dim+1); /* (dim+1) for testrectangle !!! */
  umin = _unur_vector_new(p.dim);
  umax = _unur_vector_new(p.dim);
  relative_size = _unur_vector_new(p.dim+1); /* relative size of testrectangle */

  hist_uv=malloc((p.dim+1)*p.nhist_uv*sizeof(int)); /* histograms in the (u,v)-space  */
  hist_p=malloc((p.dim+1)*p.nhist_p*sizeof(int));  /* histograms of the p-values     */
  expected_uv=_unur_vector_new(p.nhist_uv);

  // tridiagonal matrix
  for (i=0; i<p.dim; i++) {
  for (j=0; j<p.dim; j++) {
    covar[i+p.dim*j]=0;
    if (i==j) covar[i+p.dim*j]=1;
    if ((i-j)==1 || (j-i)==1) covar[i+p.dim*j]=p.rho;
//    if ((i==0 && j==(p.dim-1)) || (j==0 && i==(p.dim-1))) covar[i+p.dim*j]=p.rho;
  }}
  
  /* setting the eigenvalues of the correlation matrix */
  //for (i=0; i<p.dim; i++) eigenvalues[i]=1.; 
  /* obtaining a random correlation matrix with given eigenvalues */
  //covar_distr=unur_distr_correlation(p.dim);
  //covar_par=unur_mcorr_new(covar_distr);
  //unur_mcorr_set_eigenvalues(covar_par, eigenvalues);
  //covar_gen=unur_init(covar_par);
  //unur_sample_matr(covar_gen, covar);

  /* multinormal distribution */
  //distr = unur_distr_multinormal( p.dim, mean, covar );
  
  p.distr_normal = unur_distr_multinormal( p.dim, NULL, NULL );
  
  distr = unur_distr_cvec_new(p.dim);
  ret=unur_distr_cvec_set_mean(distr,NULL);
  ret=unur_distr_cvec_set_center(distr,NULL);
//  ret=unur_distr_cvec_set_mode(distr,mode);
  unur_distr_cvec_set_pdf(distr, pdf_mix);

  par=unur_hitrou_new(distr);
  
  unur_hitrou_set_skip(par, p.skip);
  
  if (p.shape>0) {
    /* bounding box [0,1]^(dim+1) */
    unur_hitrou_set_v(par, vmax);
    for (d=0; d<p.dim; d++) {
      umin[d]=0.; umax[d]=1.;
    }
    unur_hitrou_set_u(par, umin, umax);

    /* relative size of testrectangle */
    /* initial setting of all values to 1 */
    for (d=0; d<=p.dim; d++) {
      relative_size[d]=1.;
    }

    /* cube */
    if (p.type==0) {
      for (d=0; d<=p.dim; d++) {
        relative_size[d]=p.rel_param;
      }
    }

    /* pancake (shortest in the first (0'th) dimension) */
    if (p.type==1) relative_size[0]=p.rel_param;

    /* cigar (longest in the last dimension) */
    if (p.type==2) {
      for (d=0; d<p.dim; d++) {
        relative_size[d]=p.rel_param;
      }
    }

    /* half (shortest in the first "half" number of dimensions) */
    if (p.type==3) {
      for (d=0; d<=p.dim/2; d++) {
        relative_size[d]=p.rel_param;
      }
    }
  
    /* choosing starting point inside testrectangle */
    for (d=0; d<=p.dim; d++) {
      uv0[d]=relative_size[d]/2.; /* starting point in middle */
      uv1[d]=relative_size[d]/2.; /* starting point in middle */
      uv2[d]=1e-20; /* starting point in edge  */
    }
  
    if (p.shape==5) {
      /* ellipsoid */
      unur_hitrou_set_v(par, relative_size[p.dim]);
      for (d=0; d<p.dim; d++) {
        umin[d]=0.; umax[d]=relative_size[d];
      }
      unur_hitrou_set_u(par, umin, umax);
    }
  
  }
   
  unur_hitrou_use_bounding_rectangle(par, p.u_planes);
  unur_hitrou_set_adaptive_points(par, p.adaptive);
  
  if (p.gen) unur_free(p.gen);

  //double timing_setup, timing_sample;
  //unur_set_default_debug(0u);
  //p.gen = unur_test_timing( par, 5, &timing_setup, &timing_sample, TRUE, stdout);  
  //printf("setup=%g sample=%g\n", timing_setup, timing_sample);
  //  exit(0);
  
  p.gen=unur_init(par);
  if (p.gen==NULL) printf("p.gen=NULL\n");
  
  _unur_hitrou_set_shape(p.gen, p.shape);
  //_unur_hitrou_set_mixture_parameter(p.gen, p.mixture);

  hit_showparameters();

  printf("pdf(0)=%g\n", unur_distr_cvec_eval_pdf(x, unur_get_distr(p.gen)));

  /* resetting histogram for the 2. level test */
  memset(hist_p, 0, (p.dim+1)*p.nhist_p*sizeof(int));

  _unur_hitrou_get_point_current(p.gen, uv);
  
  for (loop=1; loop<=p.nloops; loop++) {

    /* resetting & initializing values */
    memset(hist_uv, 0, (p.dim+1)*p.nhist_uv*sizeof(int));

    _unur_hitrou_reset_pdfcount(p.gen);
    _unur_hitrou_reset_simplex_jumps(p.gen);
         
    if (p.shape==1) {
      _unur_hitrou_set_testrectangle(p.gen, relative_size);
      _unur_hitrou_set_point_current(p.gen, uv0);
    }

    if (p.shape>1) {
      _unur_hitrou_set_point_current(p.gen, uv2);
    }
    
    if (p.shape==5) {
      _unur_hitrou_set_point_current(p.gen, uv0);
    }
    
   
    /* mathematica output of points*/
    FILE *fm;
    fm=fopen("math.txt","w");
    fprintf(fm,"points={");

    
    /* sampling ... */
    for (i=1; i<=p.nsamples; i++) {

      if (p.melt==1) _unur_hitrou_set_point_current(p.gen, uv1);
      if (p.melt==2) _unur_hitrou_set_point_current(p.gen, uv2);

      unur_sample_vec(p.gen, x);
      
      /* mathematica output of points*/
      fprintf(fm,"{");
      for (j=0; j<p.dim; j++) {fprintf(fm,"%4.2f",x[j]); if (j<(p.dim-1)) fprintf(fm,",");}
      fprintf(fm,"}");
      if (i<p.nsamples ) fprintf(fm,",");    
      
      /* increase appropriate histogram bins of the uv histograms */
      if (p.shape>0) {
        for (d=0; d<=p.dim; d++) {
          iuv = (int) (p.nhist_uv*(uv[d])/(relative_size[d])) ;
          hist_uv[d*p.nhist_uv+iuv]+=1;
        }
      }
	
    } /* next sample */
    
    /* mathematica output of points*/
    fprintf(fm,"};\n");
    fclose(fm);
    

    pdfcount = _unur_hitrou_get_pdfcount(p.gen);
    sum_pdfcount += pdfcount;
    hit_separator();
    printf("loop=%ld\n", loop);
    printf("pdfcount=%ld\n", pdfcount);
    
    if (p.shape==1) {
      /* testrectangle */
      for (d=0; d<=p.dim; d++) {
        for (iuv=0; iuv<p.nhist_uv; iuv++) {
          /* expected values in each histogram bin */
          expected_uv[iuv]=p.nsamples/(double)p.nhist_uv;
        }

        pchi2= _unur_test_chi2test( expected_uv, &hist_uv[d*p.nhist_uv], p.nhist_uv, 5, p.verbose, stdout );
        printf("coord=%2d p=%8.5f \n", d, pchi2);

        /* increase appropriate bin of the p-histogram */
        if (pchi2>=0 && pchi2<=1) {
          ip = (int) (p.nhist_p*pchi2) ;
          hist_p[d*p.nhist_p+ip]+=1;
        }

      }
      hit_separator();
    }


    if (p.shape==2) {
      /* simplex */
      for (d=0; d<=p.dim; d++) {
        for (iuv=0; iuv<p.nhist_uv; iuv++) {
          r0=iuv/(double)p.nhist_uv;
          r1=(iuv+1)/(double)p.nhist_uv;
          /* expected values in each histogram bin */
          expected_uv[iuv]=p.nsamples*(pow(1-r0, p.dim+1)-pow(1-r1, p.dim+1));
        }

        pchi2= _unur_test_chi2test( expected_uv, &hist_uv[d*p.nhist_uv], p.nhist_uv, 5, p.verbose, stdout );
        printf("coord=%2d p=%8.5f \n", d, pchi2);

        /* increase appropriate bin of the p-histogram */
        if (pchi2>=0 && pchi2<=1) {
          ip = (int) (p.nhist_p*pchi2) ;
          hist_p[d*p.nhist_p+ip]+=1;
        }

      }
      hit_separator();
    }


    if (p.shape==3) {
      /* two simplex stacked along the v-direction */
      for (d=0; d<=p.dim; d++) {
        for (iuv=0; iuv<p.nhist_uv; iuv++) {
          r0=iuv/(double)p.nhist_uv;
          r1=(iuv+1)/(double)p.nhist_uv;
          /* expected values in each histogram bin */
          expected_uv[iuv]=p.nsamples*(pow(1-r0, p.dim+1)-pow(1-r1, p.dim+1));
          if (d==p.dim) {
            expected_uv[iuv]=(r1<=.5) ?
               p.nsamples*(pow(1-2*r0, p.dim+1)-pow(1-2*r1, p.dim+1)):
               p.nsamples*(pow(1-2*(r0-.5), p.dim+1)-pow(1-2*(r1-.5), p.dim+1));
          }
        }

        pchi2= _unur_test_chi2test( expected_uv, &hist_uv[d*p.nhist_uv], p.nhist_uv, 5, p.verbose, stdout );
        printf("coord=%2d p=%8.5f \n", d, pchi2);

        /* increase appropriate bin of the p-histogram */
        if (pchi2>=0 && pchi2<=1) {
          ip = (int) (p.nhist_p*pchi2) ;
          hist_p[d*p.nhist_p+ip]+=1;
        }

      }
      hit_separator();
    }


    if (p.shape==4) {
      /* two simplex stacked along the u[0]-direction */
      for (d=0; d<=p.dim; d++) {
        for (iuv=0; iuv<p.nhist_uv; iuv++) {
          r0=iuv/(double)p.nhist_uv;
          r1=(iuv+1)/(double)p.nhist_uv;
          /* expected values in each histogram bin */
          expected_uv[iuv]=p.nsamples*(pow(1-r0, p.dim+1)-pow(1-r1, p.dim+1));
          if (d==0) {
            expected_uv[iuv]=(r1<=.5) ?
               p.nsamples*(pow(1-2*r0, p.dim+1)-pow(1-2*r1, p.dim+1)):
               p.nsamples*(pow(1-2*(r0-.5), p.dim+1)-pow(1-2*(r1-.5), p.dim+1));
          }
        }

        pchi2= _unur_test_chi2test( expected_uv, &hist_uv[d*p.nhist_uv], p.nhist_uv, 5, p.verbose, stdout );
        printf("coord=%2d p=%8.5f \n", d, pchi2);

        /* increase appropriate bin of the p-histogram */
        if (pchi2>=0 && pchi2<=1) {
          ip = (int) (p.nhist_p*pchi2) ;
          hist_p[d*p.nhist_p+ip]+=1;
        }

      }
      hit_separator();
    }
    
    printf("simplex_jumps=%ld\n", _unur_hitrou_get_simplex_jumps(p.gen));
    
  } /* next loop*/

  hit_separator();
  printf("sum_pdfcount=%ld\n", sum_pdfcount);
  printf("average_pdfcount=%f\n", sum_pdfcount/((double)p.nsamples*p.nloops));
  hit_separator();

  unur_urng_reset(NULL);
  _unur_hitrou_set_point_current(p.gen, uv);
  faux.f = f_norm;
#if 1 
  for (loop=1; loop<=p.nloops; loop++) {
    _unur_expect_evaluate_single(p.dim, &faux, p.gen, p.nsamples, &integral);
//    printf("Integral = %g\n", integral);
  }      

  unur_urng_reset(NULL);
  _unur_hitrou_set_point_current(p.gen, uv);
  
  exact_integral=sqrt(2./M_PI);
  for (d=1; d<p.dim; d++) exact_integral=d/exact_integral;
  
  _unur_expect_evaluate(p.dim, &faux, p.gen, p.nsamples, p.nloops, &exact_integral, 
  &integral_mean, &integral_variance, &integral_rmse);
  printf("Mean     = %g\n", integral_mean);
  printf("Variance = %g\n", integral_variance);
  printf("RMSE     = %g\n", integral_rmse);
#endif  
   
  if (p.nloops>1 && p.shape>0 && p.melt==0) {
    printf("2.level test : \n");
    for (d=0; d<=p.dim; d++) {
      pchi2= _unur_test_chi2test( NULL, &hist_p[d*p.nhist_p], p.nhist_p, 5, p.verbose, stdout );
      printf("coord=%2d p=%8.5f \n", d, pchi2);
    }
  }

  unur_distr_free(p.distr_normal);
  unur_distr_free(distr);
  //unur_distr_free(covar_distr);
  //unur_free(covar_gen);

  free(eigenvalues);
  free(eigenvectors);
  free(mean); free(covar); free(mode);
  free(x); free(uv); free(uv0); free(uv1); free(uv2);
  free(umin), free(umax);
  free(relative_size);
  free(hist_uv); free(hist_p);
  free(expected_uv);

}

/***************************************************************************************/

int main(int argc, char *argv[])
{

  hit_init();

  char c;

  /* read options */
  while ((c = getopt(argc, argv, "d:s:n:f:r:t:l:u:p:m:a:b:v:x:c:h")) != -1) {
    switch (c) {
    case 'd':     /* dimension */
      p.dim=atoi(optarg);
      break;
    case 's':     /* skip size */
      p.skip=atoi(optarg);
      break;
    case 'f':     /* shape flag  */
      p.shape=atoi(optarg);
      break;
    case 't':     /* type of testrectangle  */
      p.type=atoi(optarg);
      break;
    case 'r':     /* relative size parameter */
      p.rel_param=atof(optarg);
      break;
    case 'x':     /* mixture parameter */
      p.mixture=atof(optarg);
      break;
    case 'c':     /* covariance rho */
      p.rho=atof(optarg);
      break;
    case 'n':     /* number of samples */
      p.nsamples=atoi(optarg);
      break;
    case 'l':     /* number of loops  */
      p.nloops=atoi(optarg);
      break;
    case 'u':     /* number of uv histogram bins */
      p.nhist_uv=atoi(optarg);
      break;
    case 'p':     /* number of 2. level histogram bins */
      p.nhist_p=atoi(optarg);
      break;
    case 'm':     /* melting */
      p.melt=atoi(optarg);
      break;
    case 'b':     /* all bounding planes */
      p.u_planes=atoi(optarg);
      break;
    case 'a':     /* line segment adaptive */
      p.adaptive=atoi(optarg);
      break;
    case 'v':     /* verbosity */
      p.verbose=atoi(optarg);
      break;
    case 'h':     /* help */
      printf("options\n" );
      printf(" -d dim        : dimension (default=2) \n" );
      printf(" -s skip       : skip parameter (default=0) \n" );
      printf(" -n samples    : number of samples (default=1) \n" );
      printf(" -f flag_shape : 0=RoU, 1=rectangle, 2=simplex, 3=2s_v, 4=2s_u0, 5=ellipsoid (default=0)\n" );
      printf(" -t type       : 0=cube/sphere, 1=pancake, 2=cigar, 3=half (default=1) \n" );
      printf(" -r rel_param  : relative size parameter (default=1.) \n" );
      printf(" -l loops      : number of loops for 2. level test, if > 1 (default=1) \n" );
      printf(" -u hist_uv    : number of histogram bins for uv-values (default=100) \n" );
      printf(" -p hist_p     : number of histogram bins for 2. level test (default=10) \n" );
      printf(" -m melt       : melting flag 0:no, 1:center, 2:edge (default=0)\n");
      printf(" -b u_planes   : calculate all bounding planes 0:no, 1:yes (default=1)\n");
      printf(" -a adaptive   : use line segment adaptive 0:no, 1:yes (default=1)\n");
      printf(" -x parameter  : mixture parameter (default=0)\n");
      printf(" -c cov_rho    : covariance parameter (default=0)\n");
    
      
      return 0;
      break;
    default:
      break;
    }
  }

  hit_run();

  if (p.gen) unur_free(p.gen);

  return 0;
}
