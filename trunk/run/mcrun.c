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

/*---------------------------------------------------------------------------*/

double funct_integral_exact (int dim );

static struct timeval tv;
#define get_time() ( gettimeofday(&tv, NULL), ((tv).tv_sec * 1.e6 + (tv).tv_usec) )

/*---------------------------------------------------------------------------*/

const int dim = 50;             /* dimensions where function lives on (<=100) */
const int sample_size = 100000;    /* sample size for MC integration         */
const int n_samples = 10;         /* total number of samples                */

int n_urng;                    /* counter for urng */

typedef double MC_COMP(int dim);

/*---------------------------------------------------------------------------*/

double funct( int dim, double *x )
     /* function to be integrated          */
     /* f(x) = \sum_{i=1}^n exp( -i*x_i^2) */
     /* domain = [0,1]^n                   */
     /*                                    */
     /* dim ... dimension (n)              */
     /* x   ... pointer to argment         */
{
  int i;
  double tmp = 0.;
  
  for (i=0; i<dim; i++)
    tmp -= (i+1.)*x[i]*x[i];

  return exp(tmp);

} /* end of funct() */

/*---------------------------------------------------------------------------*/

double uniform(void)
     /* uniform PRNG (with counter) */
{
  ++n_urng;
  return unur_urng_MRG31k3p();
} /* uniform */

void seed_uniform( int seed )
     /* seed uniform PRNG */
{
  unur_urng_MRG31k3p_seed( seed );
} /* seed_uniform() */
  
/*---------------------------------------------------------------------------*/

void print_result( char *title, int dim, int n_samples, 
		   double mean, double stddev, double time )
     /* print result */
{
  double exact = funct_integral_exact(dim);

  printf ("%s\n", title);
  printf ("\tdim     = %d\n", dim);
  printf ("\tsamples = %d \t (sample size = %d)\n", n_samples, sample_size);
  printf ("\t#urn    = %g\n", ((double)n_urng)/(n_samples*sample_size));
  printf ("\tmean    = %g\n", mean);
  printf ("\tstddev  = %g (%g%%)\n", stddev, 100.*stddev/mean);
  printf ("\texact   = %g\n", exact);
  printf ("\terror   = %g (%g%%)\n", mean-exact, 100.*(mean-exact)/exact);
  printf ("\ttime    = %d ms\n", (int)((time/n_samples)/1.e3));
  printf ("\n");
} /* end of print_result() */

/*---------------------------------------------------------------------------*/

double MC_naive( int dim )
{
  double *x;         /* aux array for storing vector */
  int i;             /* counting dimensions     */ 
  int n_sample;      /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double result;     /* result of MC integration */

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  for (n_sample = 0; n_sample < sample_size; ++n_sample) {
    /* make point uniformly distributed in [0,1]^dim */
    for (i = 0; i<dim; i++) x[i] = uniform(); 
    /* compute function value */
    sum += funct(dim,x);
  }

  /* free working space and exit */
  free(x);

  /*result */
  result = sum / n_sample;

  /* return result */
  return result;

} /* end of MC_naive() */

/*---------------------------------------------------------------------------*/

double MC_importance( int dim )
{
  double *x;         /* aux array for storing vector */
  double fx;         /* store value of density for importance sampling */
  int i;             /* counting dimensions     */ 
  int n_sample;      /* counter for sample size */
  double sum = 0.;   /* store sum of function values */
  double result;     /* result of MC integration */
  double norm;
  UNUR_DISTR *distr;
  UNUR_GEN *gen;

  /* allocate working space */
  x = malloc (dim * sizeof(double));

  /* get generator for marginal distributions */
  distr = unur_str2distr("exponential");
  gen = unur_str2gen("exponential & method=tdr; cpoints=(1); usedars=on; max_sqhratio=0.5");
  
  /* normalization constant for truncated density */
  norm = 1.;
  for (i = 0; i<dim; i++)
    norm /= unur_distr_cont_eval_cdf(sqrt(i+1.),distr) / sqrt(i+1.);

  for (n_sample = 0; n_sample < sample_size; ++n_sample) {

    /* make point with importance distribution in [0,1]^dim */
    fx = 1.;
    for (i = 0; i<dim; i++) {
      unur_tdr_chg_truncated(gen,0,sqrt(i+1.));
      x[i] = unur_sample_cont(gen);
      fx *= unur_distr_cont_eval_pdf(x[i],distr);
      x[i] /= sqrt(i+1.);
    }
    /* compute function value */
    sum += funct(dim,x) / (norm * fx);
  }

  /* free working space and exit */
  free(x);
  unur_distr_free(distr);
  unur_free(gen);

  /* result */
  result = sum / sample_size;

  /* return result */
  return result;

} /* end of MC_importance() */

/*---------------------------------------------------------------------------*/

int run_MC (MC_COMP *MC, char *title, int dim, int n_samples)
{
  double integral;
  double sum, sumsqu;
  double mean, stddev;
  double start, stop;
  int sample;

  /* reset counter for urng */
  n_urng = 0;

  /* init */
  sum = sumsqu = 0.;
  start = get_time();

  for (sample = 0; sample < n_samples; ++sample) {
    integral = MC( dim );
    sum += integral;
    sumsqu += integral*integral;
  }

  stop = get_time();

  mean = sum / n_samples;
  stddev = sqrt(sumsqu / n_samples - mean*mean);

  print_result( title, dim, n_samples, mean, stddev, stop-start );

  return 1;

} /* end of run_MC() */

/*---------------------------------------------------------------------------*/

int main()
{
  /* set and seed uniform PRNG */
  seed_uniform(123);
  unur_set_default_urng (uniform);


  run_MC(MC_naive, "Naive MC", dim, n_samples);
  run_MC(MC_importance, "Importance MC", dim, n_samples);

  /* free working space and exit */
  exit (EXIT_SUCCESS);

}

/*---------------------------------------------------------------------------*/

double funct_integral_exact (int dim )
{
  double data[100] = 
    { 0.74682413281242703, 0.59814400666130410, 0.50434356023143881, 0.44104069538121084, \
      0.39571230961051354, 0.36160814735365850, 0.33490105817655928, 0.31330868732130717, \
      0.29540244941984041, 0.28024739050664274, 0.26720674335186228, 0.25583143052938306, \
      0.24579504080564116, 0.23685407997867067, 0.22882279832973735, 0.22155672794739224, \
      0.21494160010412462, 0.20888568914041529, 0.20331440033476209, 0.19816636482997365, \
      0.19339056992497953, 0.18894421535395370, 0.18479108808019950, 0.18090031363879588, \
      0.17724538509027910, 0.17380339947509397, 0.17055445132438055, 0.16748114642265756, \
      0.16456820862355452, 0.16180215937964007, 0.15917105461020218, 0.15666426716443734, \
      0.15427230582764102, 0.15198666383034363, 0.14979969134027405, 0.14770448757545967, \
      0.14569480906717378, 0.14376499129129276, 0.14190988142510872, 0.14012478040994822, \
      0.13840539283493047, 0.13674778342396566, 0.13514833912180541, 0.13360373594713928, \
      0.13211090992020037, 0.13066703148589484, 0.12926948294637178, 0.12791583849331106, \
      0.12660384649325114, 0.12533141373155003, 0.12409659136408728, 0.12289756236218159, \
      0.12173263026670363, 0.12060020909304461, 0.11949881425029427, 0.11842705435636683, \
      0.11738362384644484, 0.11636729628544092, 0.11537691830657834, 0.11441140410797112, \
      0.11346973044749433, 0.11255093208348863, 0.11165409761511313, 0.11077836568159475, \
      0.10992292148434307, 0.10908699360000995, 0.10826985105616018, 0.10747080064435617, \
      0.10668918444820866, 0.10592437756635953, 0.10517578601248646, 0.10444284477629169, \
      0.10372501603109048, 0.10302178747507787, 0.10233267079464885, 0.10165720023929806, \
      0.10099493129864889, 0.10034543947307326, 0.099708319130176738, 0.099083182440150275, \
      0.098469658383639779, 0.097867391826367331, 0.097276042655260457, 0.096695284971315475, \
      0.096124806334843435, 0.095564307059127779, 0.095013499548866183, 0.094472107680079087, \
      0.093939866218447792, 0.093416520273298812, 0.092901824784681209, 0.092395544041192357, \
      0.091897451226397168, 0.091407327991858247, 0.090924964054951337, 0.090450156819783457, \
      0.089982711019661553, 0.089522438379678589, 0.089069157298092852, 0.088622692545275801 };

  double result = 1.;
  int i;

  for (i=0; i<dim; i++) 
    result *= data[i];

  return result;
  
} /* end of funct_integral_exact() */

/*---------------------------------------------------------------------------*/













