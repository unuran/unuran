/* UNURAN program example_cont.c                                   */

#include <unuran.h>

/* --------------------------------------------------------------- */

/* PDF -- piecewise linear                                         */
double mypdf( double x, UNUR_DISTR *dummy )
{
  double fx;

  if ( x <= 0.0 )
    fx = 0.0;
  else if ( x <= 1.0 )
    fx = 0.5 * x; 
  else if ( x <= 2.0 )
    fx = 0.5;
  else if ( x <= 3.0 )
    fx = -0.5 * x + 1.5;
  else
    fx = 0.0;

  return fx;
} /* end of mypdf                                                  */

/* DPDF -- the derivative of the PDF                               */
double mydpdf ( double x, UNUR_DISTR *dummy)
{
  double dfx;

  if ( x <= 0.0 )
    dfx = 0.0;
  else if ( x <= 1.0 )
    dfx = 0.5; 
  else if ( x <= 2.0 )
    dfx = 0.0;
  else if ( x <= 3.0 )
    dfx = -0.5;
  else
    dfx = 0.0;

  return dfx;
} /* end of mydpdf                                                 */

/* --------------------------------------------------------------- */

int main()
{
  int    i;

  UNUR_DISTR *distr;  /* distribution                              */
  UNUR_PAR   *par;    /* parameter                                 */
  UNUR_GEN   *gen;    /* generator                                 */
  UNUR_URNG  *ug;     /* uniform generator object                  */

  /* Create an empty distribution object                           */
  distr = unur_distr_cont_new();

  /* "Fill" the distribution object -- the provided information
     must fulfill the requirements of the method choosen below.    */
  unur_distr_cont_set_pdf(distr,  mypdf);    /* set pdf ...        */
  unur_distr_cont_set_dpdf(distr, mydpdf);   /* ... and derivative */
  unur_distr_cont_set_mode(distr, 1.5);      /* set mode           */

  /* Choose method "TDR" and create parameter object               */
  par = unur_tdr_new(distr);

  /* Set some parameters of the method TDR                         */
  unur_tdr_set_variant_ia(par);
  unur_tdr_set_max_sqhratio(par, 0.99);
  unur_tdr_set_c(par, -0.5);
  unur_tdr_set_max_intervals(par, 1000);
  unur_tdr_set_cpoints(par, 10, NULL);

  fprintf(stderr,"junk 5\n");
  
  /* creating generator object -- now sampling is possible         */
  gen = unur_init(par);
  /* Always perform this check                                     */
  if ( gen == NULL ){
     fprintf(stderr, "Error creating generation object\n");
     return (1);
  }

  /* destroy distribution object                                   */
  unur_distr_free(distr);

  /* sample some random numbers                                    */
  for (i=0;i<100;i++)
     printf( "%f\n", unur_sample_cont(gen) );

  /* destroy generator object                                      */
  unur_free(gen);

  return (0);
}


