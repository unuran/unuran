/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Examples                                                                 *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unuran.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

/*---------------------------------------------------------------------------*/
/* uniform pseudo- and quasi-random number generators                        */

static int n_urng;
/* counter for urng */

double uniform_prng(void);
/* uniform PRNG (with counter) */

void seed_uniform_prng(int seed); 
/* seed uniform PRNG */

/*---------------------------------------------------------------------------*/
/* Run MC computations                                                       */

typedef UNUR_GEN *MAKE_GEN(UNUR_DISTR *distr);
typedef double MC_COMP(MAKE_GEN *make_gen);
/* type for function */

int run_MC (MC_COMP *MC, MAKE_GEN *make_gen, char *title);
/* master routine for running MC simulations */

double MC_naive(MAKE_GEN *make_gen);
/* naive MC, uniformly distributed points in unit cube */

double MC_importance(MAKE_GEN *make_gen);
UNUR_GEN *use_inversion(UNUR_DISTR *distr);
UNUR_GEN *use_tdr_bad(UNUR_DISTR *distr);
UNUR_GEN *use_tdr_better(UNUR_DISTR *distr);
UNUR_GEN *use_tdr_good(UNUR_DISTR *distr);
UNUR_GEN *use_tdr_excellent(UNUR_DISTR *distr);
/* use importance sampling, use exponential distribution */

double MC_hat_importance(MAKE_GEN *make_gen);
/* use importance sampling, use tdr hat for exponential distribution */

double MC_xhat_importance(MAKE_GEN *make_gen);
/* use importance sampling, use tdr, dim+1 urng, use exponential distribution */

double MC_xsmoothed_importance(MAKE_GEN *make_gen);
/* use importance sampling, use smoothed tdr, dim+1 urng, use exponential distribution */

double MC_const_importance(MAKE_GEN *make_gen);
/* use importance sampling, rejection from constant hat for exponential distribution */

double MC_wus(MAKE_GEN *make_gen);
/* use weighted uniform sampling, use exponential distribution */

/*---------------------------------------------------------------------------*/
/* timings                                                                   */

static struct timeval tv;
#define get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )

/*---------------------------------------------------------------------------*/
/* printing                                                                  */

void print_result( char *title, double xsum, double xsumsqu, double time );
/* print result */


/*---------------------------------------------------------------------------*/

static int sample_size = 10000;    /* sample size for MC integration         */
static int n_samples = 10;         /* total number of samples                */

/*---------------------------------------------------------------------------*/
/* Test functions                                                            */

/* dimension */
static int dim = 5;

/* domain for function to be integrated [0,domain)^dim */
static double domain = 3.;

static UNUR_DISTR *distr_normal = NULL;
  
double funct( double *x )
{
  int i;
  double tmp = 1.;
  
  if (distr_normal == NULL) distr_normal = unur_distr_normal(NULL,0);
  
  for (i=0; i<dim; i++)
    tmp *= unur_distr_cont_eval_pdf(x[i],distr_normal);
  return tmp;
} /* end of funct() */

double funct_integral_exact( void )
{
  static double result = -1.;

  if (result < 0.) { 
    if (distr_normal == NULL) distr_normal = unur_distr_normal(NULL,0);
    result = pow(unur_distr_cont_eval_cdf(domain,distr_normal)-0.5,dim);
  }

  return result;
} /* end of funct_integral_exact() */

double funct_domain_area( void )
{
  static double area = -1;
  int i;
  if (area < 0.) {
    area = 1.;
    for (i=0; i<dim; i++) area *= domain;
  }
  return area;

} /* end of funct_domain_area() */

/*---------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
  char c;
  int seed = 123;

  /* read options                                                            */
  while ((c = getopt(argc, argv, "n:m:d:b:s:")) != -1) {
    switch (c) {
    case 'n':     /* sample size */
      sample_size = atoi(optarg);
      break;
    case 'm':     /* number of samples */
      n_samples = atoi(optarg);
      break;
    case 'd':     /* dimension */
      dim = atoi(optarg);
      if (dim<=0) exit (EXIT_FAILURE);
      break;
    case 'b':     /* right boundary of domain  */
      domain = atof(optarg);
      if (domain<=0.) exit (EXIT_FAILURE);
      break;
    case 's':     /* seed */
      seed = atoi(optarg);
      if (seed<=0) seed *= -1;
      break;
    default:
      exit (EXIT_FAILURE);
    }
  }

  /* set and seed uniform PRNG */
  seed_uniform_prng(seed);
  unur_set_default_urng (uniform_prng);

  run_MC(MC_naive, NULL, "Naive MC");

/*    run_MC(MC_importance, use_inversion, "Importance MC, inversion"); */

/*    run_MC(MC_importance, use_tdr_bad, "Importance MC, TDR, bad hat"); */
/*    run_MC(MC_importance, use_tdr_better, "Importance MC, TDR, better hat"); */
/*    run_MC(MC_importance, use_tdr_good, "Importance MC, TDR, good hat"); */
/*    run_MC(MC_importance, use_tdr_excellent, "Importance MC, TDR, excellent hat"); */

/*    run_MC(MC_hat_importance, use_tdr_bad, "hat importance MC, TDR, bad hat"); */
/*    run_MC(MC_hat_importance, use_tdr_better, "hat importance MC, TDR, better hat"); */
/*    run_MC(MC_hat_importance, use_tdr_good, "hat importance MC, TDR, good hat"); */
/*    run_MC(MC_hat_importance, use_tdr_excellent, "hat importance MC, TDR, excellent hat"); */

/*    run_MC(MC_xhat_importance, use_tdr_bad, "single importance MC, TDR, bad hat"); */
  run_MC(MC_xhat_importance, use_tdr_better, "single importance MC, TDR, better hat");
/*    run_MC(MC_xhat_importance, use_tdr_good, "single importance MC, TDR, good hat"); */
/*    run_MC(MC_xhat_importance, use_tdr_excellent, "single importance MC, TDR, excellent hat"); */

/*    run_MC(MC_xsmoothed_importance, use_tdr_better, "single importance MC, smoothed TDR, better hat"); */

/*    run_MC(MC_const_importance, NULL, "Importance MC, rejection from constant hat"); */

/*    run_MC(MC_wus, NULL, "weighted uniform sampling MC"); */

  /* free working space and exit */
  exit (EXIT_SUCCESS);

  return EXIT_SUCCESS;
} /* end of main() */

/*---------------------------------------------------------------------------*/

double uniform_prng(void)
     /* uniform PRNG (with counter) */
{
  ++n_urng;
  return unur_urng_MRG31k3p();
} /* uniform_prng() */

void seed_uniform_prng( int seed )
     /* seed uniform PRNG */
{
  unur_urng_MRG31k3p_seed( seed );
} /* seed_uniform_prng() */
  
/*---------------------------------------------------------------------------*/

int run_MC (MC_COMP *MC, MAKE_GEN *make_gen, char *title)
{
  double integral;
  double sum, sumsqu;
  double start, stop;
  int m;

  /* reset counter for urng */
  n_urng = 0;

  /* init */
  sum = sumsqu = 0.;
  start = get_time();

  for (m = 0; m < n_samples; ++m) {
    integral = MC(make_gen);
    sum += integral;
    sumsqu += integral*integral;
  }

  stop = get_time();

  print_result( title, sum, sumsqu, stop-start );

  return 1;
} /* end of run_MC() */

/*---------------------------------------------------------------------------*/

double MC_naive( MAKE_GEN *dummy )
{
  double *x;         /* aux array for storing vector */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double area = funct_domain_area();

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  for (n = 0; n < sample_size; ++n) {
    /* make point uniformly distributed in [0,domain)^dim */
    for (i = 0; i<dim; i++) x[i] = domain * uniform_prng(); 
    /* compute function value */
    sum += funct(x);
  }

  /* free working space and exit */
  free(x);

  /* return result */
  return (area * sum/sample_size);

} /* end of MC_naive() */

/*---------------------------------------------------------------------------*/

double MC_importance( MAKE_GEN *make_gen )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double norm;
  UNUR_DISTR *distr;
  UNUR_GEN *gen;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  gen = make_gen(distr);
  
  /* normalization constant for truncated density */
  norm = pow(unur_distr_cont_eval_cdf(domain,distr),-dim);

  for (n = 0; n < sample_size; ++n) {
    /* make point with importance distribution */
    fx = 1.;
    for (i = 0; i<dim; i++) {
      x[i] = unur_sample_cont(gen);
      fx *= unur_distr_cont_eval_pdf(x[i],distr);
    }
    /* compute function value */
    sum += funct(x) / (norm * fx);
  }

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);
  unur_free(gen);

  /* return result */
  return (sum / sample_size);

} /* end of MC_importance() */

UNUR_GEN *use_inversion(UNUR_DISTR *distr)
{
  UNUR_PAR *par = unur_cstd_new(distr);
  unur_cstd_set_variant(par,UNUR_STDGEN_INVERSION);
  return unur_init(par);
}

UNUR_GEN *use_tdr_bad(UNUR_DISTR *distr)
{
  UNUR_PAR *par = unur_tdr_new(distr);
  unur_tdr_set_cpoints(par,1,NULL);
  unur_tdr_set_usedars(par,0);
  unur_tdr_set_max_sqhratio(par,0.5);
  unur_tdr_set_variant_ps(par);
  unur_tdr_set_c(par,-0.5);
  return unur_init(par);
}

UNUR_GEN *use_tdr_better(UNUR_DISTR *distr)
{
  UNUR_PAR *par = unur_tdr_new(distr);
  unur_tdr_set_cpoints(par,3,NULL);
  unur_tdr_set_usedars(par,1);
  unur_tdr_set_max_sqhratio(par,0.9);
  unur_tdr_set_variant_ps(par);
  unur_tdr_set_c(par,-0.5);
  return unur_init(par);
}

UNUR_GEN *use_tdr_good(UNUR_DISTR *distr)
{
  UNUR_PAR *par = unur_tdr_new(distr);
  unur_tdr_set_cpoints(par,3,NULL);
  unur_tdr_set_usedars(par,1);
  unur_tdr_set_max_sqhratio(par,0.99);
  unur_tdr_set_variant_ps(par);
  unur_tdr_set_c(par,-0.5);
  return unur_init(par);
}

UNUR_GEN *use_tdr_excellent(UNUR_DISTR *distr)
{
  UNUR_PAR *par = unur_tdr_new(distr);
  unur_tdr_set_cpoints(par,3,NULL);
  unur_tdr_set_usedars(par,1);
  unur_tdr_set_max_sqhratio(par,0.999);
  unur_tdr_set_variant_ps(par);
  unur_tdr_set_c(par,-0.5);
  return unur_init(par);
}

/*---------------------------------------------------------------------------*/

double MC_hat_importance( MAKE_GEN *make_gen )
{
  double *x;         /* aux array for storing vector */
  double hx;         /* store value of hat of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double norm;
  UNUR_DISTR *distr;
  UNUR_GEN *gen;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  gen = make_gen(distr);
  
  /* normalization constant for hat of truncated density */
  norm = pow(unur_tdr_get_hatarea(gen),-dim);

  for (n = 0; n < sample_size; ++n) {
    /* make point with importance distribution */
    hx = 1.;
    for (i = 0; i<dim; i++) {
      double fxi,hxi,sqxi;
      x[i] = unur_tdr_eval_invcdfhat(gen,uniform_prng(),&hxi,&fxi,&sqxi);
      hx *= hxi;
    }
    /* compute function value */
    sum += funct(x) / (norm * hx);
  }

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);
  unur_free(gen);

  /* return result */
  return (sum / sample_size);

} /* end of MC_hat_importance() */

/*---------------------------------------------------------------------------*/

double MC_xhat_importance( MAKE_GEN *make_gen )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  double hx;         /* store value of hat of density for importance sampling */
  double sqx;        /* store value of squeeze of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double norm;
  double V;
  UNUR_DISTR *distr;
  UNUR_GEN *gen;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  gen = make_gen(distr);
  
  /* normalization constant for hat of truncated density */
  norm = pow(unur_distr_cont_eval_cdf(domain,distr),-dim);

  for (n = 0; n < sample_size; ++n) {
    /* make point with importance distribution */
    while (1) {
      hx = fx = sqx = 1.;
      for (i = 0; i<dim; i++) {
	double fxi,hxi,sqxi;
	x[i] = unur_tdr_eval_invcdfhat(gen,uniform_prng(),&hxi,&fxi,&sqxi);
	hx *= hxi;
	fx *= fxi;
	sqx *= sqxi;
      }
      V = uniform_prng() * hx;
      if (V <= sqx || V <= fx) break;
    }
    /* compute function value */
    sum += funct(x) / (norm * fx);
  }

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);
  unur_free(gen);

  /* return result */
  return (sum / sample_size);

} /* end of MC_xhat_importance() */

/*---------------------------------------------------------------------------*/

double MC_xsmoothed_importance( MAKE_GEN *make_gen )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  double hx;         /* store value of hat of density for importance sampling */
  double sqx;        /* store value of squeeze of density for importance sampling */
  int i;             /* counting dimensions     */ 
  double sum = 0.;   /* store sum of function values */
  double w;          /* weight of generated point */
  double wsum = 0.;  /* sum of weights */
  double norm;
  double V;
  UNUR_DISTR *distr;
  UNUR_GEN *gen;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  gen = make_gen(distr);
  
  /* normalization constant for hat of truncated density */
  norm = pow(unur_distr_cont_eval_cdf(domain,distr),-dim);

  /* make point with importance distribution */
  do {
    hx = fx = sqx = 1.;
    for (i = 0; i<dim; i++) {
      double fxi,hxi,sqxi;
      x[i] = unur_tdr_eval_invcdfhat(gen,uniform_prng(),&hxi,&fxi,&sqxi);
      hx *= hxi;
      fx *= fxi;
      sqx *= sqxi;
    }
    V = uniform_prng() * hx;
    /* compute weights */
    if (V <= sqx || V <= fx) 
      w = 1.;
    else
      w = 0.;
    wsum += w;
    /* compute function value */
    sum += w * funct(x) / (norm * fx);
  } while (wsum < sample_size-0.5);

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);
  unur_free(gen);

  /* return result */
  return (sum / sample_size);

} /* end of MC_xsmoothed_importance() */

/*---------------------------------------------------------------------------*/

double MC_const_importance( MAKE_GEN *dummy )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double M;
  UNUR_DISTR *distr;
  double norm;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  
  /* normalization constant for hat of truncated density */
  norm = pow(unur_distr_cont_eval_cdf(domain,distr),-dim);
  M = unur_distr_cont_eval_pdf(0.,distr);

  for (n = 0; n < sample_size; ++n) {
    do {
      /* make point uniformly distributed in [0,domain)^dim */
      fx = 1.;
      for (i = 0; i<dim; i++) {
	x[i] = domain * uniform_prng(); 
	fx *= unur_distr_cont_eval_pdf(x[i],distr);
      } 
    } while (uniform_prng() * M >= fx);
    /* compute function value */
    sum += funct(x) / (norm * fx);
  }

  /* free working space and exit */
  free(x);

  /* return result */
  return (sum/sample_size);

} /* end of MC_const_importance() */

/*---------------------------------------------------------------------------*/

double MC_wus( MAKE_GEN *dummy )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n;             /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double fsum = 0.;  /* store sum of density values */
  double norm;
  UNUR_DISTR *distr;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* importance distribution */
  distr = unur_distr_exponential(NULL,0);
  unur_distr_cont_set_domain(distr,0.,domain);
  
  /* normalization constant for truncated density */
  norm = pow(unur_distr_cont_eval_cdf(domain,distr),-dim);

  for (n = 0; n < sample_size; ++n) {
    fx = 1.;
    for (i = 0; i<dim; i++) {
      x[i] = domain * uniform_prng(); 
      fx *= unur_distr_cont_eval_pdf(x[i],distr);
    }
    /* compute function value */
    sum += funct(x);
    fsum += fx;
  }

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);

  /* return result */
  return (sum / (norm * fsum));

} /* end of MC_wus() */

/*---------------------------------------------------------------------------*/

void print_result( char *title, double xsum, double xsumsqu, double time )
     /* print result */
{
  double mean, stddev, rmse;
  double exact = funct_integral_exact();

  static double stddev_unit = -1.;
  static double rmse_unit = -1.;

  mean = xsum / n_samples;
  stddev = sqrt((xsumsqu - n_samples*mean*mean)/(n_samples-1));
  rmse = sqrt((xsumsqu - 2*exact*xsum + n_samples*exact*exact)/n_samples);
  mean /= exact;
  stddev /= exact;
  rmse /= exact;

  if (stddev_unit < 0.) stddev_unit = stddev;
  if (rmse_unit < 0.) rmse_unit = rmse;

  printf ("%s\n", title);
  printf ("\tdim     = %d\n", dim);
  printf ("\tdomain  = [0,%g)\n", domain);
  printf ("\tsamples = %d \t (sample size = %d)\n", n_samples, sample_size);
  printf ("\t#urn    = %g\n", ((double)n_urng)/(n_samples*sample_size));
  printf ("\texact   = %g\n", 1.);
  printf ("\tmean    = %g\n", mean);
  printf ("\terror   = %g\n", fabs(mean-1.));
  printf ("\tstddev  = %g\t(Reff = %g)\n", stddev,stddev_unit/stddev);
  printf ("\trmse    = %g\t(Reff = %g)\n", rmse,rmse_unit/rmse);
  printf ("\ttime    = %d ms\n", (int)((time/n_samples)/1.e3));
  printf ("\n");
} /* end of print_result() */

/*---------------------------------------------------------------------------*/
