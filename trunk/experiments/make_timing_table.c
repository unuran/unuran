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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include <unuran.h>
#include <unuran_tests.h>

/* ------------------------------------------------------------------------- */

/* Program name                                                              */
static char *progname;

/*---------------------------------------------------------------------------*/

extern char *_unur_parser_prepare_string( const char *str );

/*---------------------------------------------------------------------------*/

/* get unit for relative timing */
double get_timing_unit(void);

/* read config file */
int read_config_file ( const char *filename, 
		       struct unur_slist *distr_str_list, struct unur_slist *meth_str_list ); 

/* make distribution objects */
struct unur_slist *make_distr_list ( struct unur_slist *distr_str_list );

/* print legend for distributions and methods */
int print_legend ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list ); 

/* print label for distribution and methods with index n*/
int print_label ( int n, char ltype );

#define LABEL_DISTR 'A'  /* label used for distributions */
#define LABEL_METH  '1'  /* label used for methods */

/* compute timings table */
double *compute_timings ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list,
			  int samplesize, double duration ); 

/* print timings table */
int print_timings ( double *timings,
		    struct unur_slist *distr_str_list, struct unur_slist *meth_str_list,
		    int samplesize, int rowentry );

#define ROW_DISTRIBUTION   1  /* print distributions on row, methods in columns */
#define ROW_METHOD         2  /* print distributions on columns, methods in rows*/

/* print usage */
void print_usage(void);

/*---------------------------------------------------------------------------*/

int main (int argc, char *argv[])
{
  /* defaults */
  int samplesize = 1000;    /* size of samples */
  double duration = 0.1;    /* duration in seconds for timing generation of a sample */

  struct unur_slist *distr_str_list; /* list of strings for distributions */
  struct unur_slist *meth_str_list;  /* list of strings for methods */
  char *conffile;                    /* name of configuration file */

  double *time_0;


  /* read parameters */
  progname = argv[0];
  if (argc<2) {
    print_usage();
    exit (EXIT_FAILURE);
  }
  /* name of configuration file */
  conffile = argv[1];
  
  /* sample size */
  if (argc >= 3)
    samplesize = atoi(argv[2]);

  /* switch off all debugging and logging information*/
  unur_set_default_debug(0u);

  /* create lists for distributions and methods */
  distr_str_list = _unur_slist_new();
  meth_str_list  = _unur_slist_new();

  /* read config file */
  read_config_file(conffile, distr_str_list, meth_str_list);

  /* print legend */
  print_legend(distr_str_list,meth_str_list);
  
  /* make timings */
  time_0 = compute_timings(distr_str_list,meth_str_list,samplesize,duration);

  /* print timings */
  print_timings(time_0,distr_str_list,meth_str_list,samplesize,ROW_DISTRIBUTION);
  print_timings(time_0,distr_str_list,meth_str_list,samplesize,ROW_METHOD);

  /* free memory */
  _unur_slist_free(distr_str_list);
  _unur_slist_free(meth_str_list);
  free(time_0);

  exit (EXIT_SUCCESS);

} /* end of main() */

/*---------------------------------------------------------------------------*/

double
get_timing_unit(void)
     /* get unit for relative timings                                        */
     /* (use generation of exponential random variate via inversion)         */
{
  UNUR_DISTR *distr;    /* pointer to working distribution object */
  UNUR_PAR *par;        /* pointer to working parameter object */
  double timing_unit;   /* timing result for basis of relative timings */

  distr = unur_distr_exponential(NULL,0);
  par = unur_cstd_new(distr);
  timing_unit = unur_test_timing_exponential(par, 5);
  free(par);

  return timing_unit;
} /* end of get_timing_unit() */

/*---------------------------------------------------------------------------*/

int
read_config_file ( const char *filename, 
		   struct unur_slist *distr_str_list, struct unur_slist *meth_str_list )
     /* read config file */
{
#define LINELENGTH  1024      /* max length of lines allowed    */

  char line[LINELENGTH];      /* input buffer */
  char *str;                  /* pointer to working string */

  FILE *fh;                   /* file handle for input stream */
  int line_no;                /* counter for lines */

  /* open the file with the data read-only */
  fh = fopen(filename, "r");
  if (fh == NULL) {
    fprintf(stderr,"error: cannot open config file `%s'.\n",filename);
    exit(EXIT_FAILURE);
  }

  /* read lines until eof */
  line_no = 0;
  while( fgets(line, LINELENGTH, fh) ) {
    ++line_no;

    /* make a working copy of the line;                           */
    /* remove all white spaces and convert to lower case letters. */
    str = _unur_parser_prepare_string( line );

    /* ignore all lines that do not start with a letter */
    if ( ! isalpha(str[0]) || str[0] == '#') {
      free(str);
      continue;
    }
    
    /* store distribution object */
    if ( !strncmp( str, "distr", 5) ) {
      _unur_slist_append(distr_str_list,str);
      continue;
    }

    /* store method (parameter object) */
    if ( !strncmp( str, "method=", 7) ) {
      _unur_slist_append(meth_str_list,str);
      continue;
    }

    fprintf(stderr,"syntax error in line %d: %s\n",line_no,str);
    exit (EXIT_FAILURE);
  }

  /* close input stream */
  fclose(fh);

  /* check number of distributions and methods */
  if ( _unur_slist_length(distr_str_list) == 0 ||
       _unur_slist_length(meth_str_list)  == 0 ) {
    fprintf(stderr,"error: no distributions or methods given.\n");
    exit (EXIT_FAILURE);
  }

  /* o.k. */
  return 1;

#undef LINELENGTH

} /* end of read_config_file() */

/*---------------------------------------------------------------------------*/

struct unur_slist *
make_distr_list ( struct unur_slist *distr_str_list )
     /* make distribution objects */
{
  struct unur_slist *distr_list; /* list of distributions */
  int n_distr;        /* number of distributions */
  char *str;          /* pointer to working string */
  UNUR_DISTR *distr;  /* pointer to distribution object */
  int i;

  /* get number of distributions */
  n_distr = _unur_slist_length(distr_str_list);

  /* create lists for distributions and methods */
  distr_list = _unur_slist_new();

  /* get all distribution objects */
  for (i=0; i<n_distr; i++) {
    str = _unur_slist_get(distr_str_list,i);
    distr = unur_str2distr(str);
    if (distr == NULL) {
      fprintf(stderr,"syntax error: %s\n",str);
      exit (EXIT_FAILURE);
    }
    else
      _unur_slist_append(distr_list,distr);
  }

  /* o.k. */
  return distr_list;
} /* end of make_distr_list() */

/*---------------------------------------------------------------------------*/

int 
print_legend ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list )
     /* print legend for distributions and methods */
{
  int n_distr;        /* number of distributions */
  int n_meth;         /* number of methods */
  char *str;          /* pointer to working string */
  int i,k;

  /* get number of distributions and methods */
  n_distr = _unur_slist_length(distr_str_list);
  n_meth  = _unur_slist_length(meth_str_list);

  /* print legend */
  printf("%d distributions:\n",n_distr);
  for (i=0; i<n_distr; i++) {
    str = _unur_slist_get(distr_str_list, i);
    print_label(i, LABEL_DISTR);
    printf(" ... %s\n", str);
  }
  printf("\n");

  printf("%d methods:\n",n_meth);
  for (k=0; k<n_meth; k++) {
    str = _unur_slist_get(meth_str_list, k);
    print_label(k, LABEL_METH);
    printf(" ... %s\n", str);
  }
  printf("\n");

  /* o.k. */
  return 1;
} /* end of print_legend() */

/*---------------------------------------------------------------------------*/

int
print_label ( int n, char ltype )
     /* print label for distribution and methods with index n*/
{
  switch(ltype) {
  case 'a':
    printf(" [%c]", 'a'+n);
    break;
  case 'A':
    printf(" [%c]", 'A'+n);
    break;
  case '1':
  default:
    printf("[%02d]", 1+n);
    break;
  }

  return 1;
} /* end of print_label() */

/*---------------------------------------------------------------------------*/

double *
compute_timings ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list,
		  int samplesize, double duration )
     /* compute timings table */
{
  struct unur_slist *distr_list; /* list of distributions */
  int n_distr;        /* number of distributions */
  int n_meth;         /* number of methods */
  char *str;          /* pointer to working string */
  UNUR_DISTR *distr;  /* pointer to working distribution object */
  UNUR_PAR *par;      /* pointer to working parameter object */
  int i,k;
  struct unur_slist *mlist;  /* list of allocated memory blocks in running _unur_str2par() */
  double *timing;     /* timing results */
  double timing_unit; /* unit for timing result */

  /* get unit for relative timings */
  timing_unit = get_timing_unit();

  /* get all distribution objects */
  distr_list = make_distr_list(distr_str_list);

  /* get number of distributions and methods */
  n_distr = _unur_slist_length(distr_str_list);
  n_meth  = _unur_slist_length(meth_str_list);

  /* allocate array for timings */
  timing = malloc(n_distr * n_meth * sizeof(double));

  /* make timings */
  for (i=0; i<n_distr; i++) {
    distr = _unur_slist_get(distr_list, i);
    for (k=0; k<n_meth; k++) {
      str = _unur_slist_get(meth_str_list, k);
      par = _unur_str2par(distr, str, &mlist);
      if (par) {
	timing[i*n_meth+k] = unur_test_timing_total(par, samplesize, duration );
	timing[i*n_meth+k] /= samplesize * timing_unit;
	free(par);
	_unur_slist_free(mlist);
      }
      else {
	timing[i*n_distr+k] = -1.;  /* no timing result */
      }
    }  
  }

  /* free memory */
  _unur_slist_free(distr_list);

  return timing;
} /* end of compute_timings() */

/*---------------------------------------------------------------------------*/

int 
print_timings ( double *timings,
		struct unur_slist *distr_str_list, struct unur_slist *meth_str_list,
		int samplesize, int rowentry )
     /* print timings table */
{
  int n_distr, n_meth;          /* number of distributions and methods */
  int n_row, n_col;             /* number of elements in rows and columns */
  char rltype, cltype;          /* label types for rows and columns */
  int row, col;                 /* indices for row and columns */
  int idx_time;                 /* position in timing table */

  /* get number of distributions and methods */
  n_distr = _unur_slist_length(distr_str_list);
  n_meth  = _unur_slist_length(meth_str_list);

  /* where to put distributions */
  switch (rowentry) {
  case ROW_DISTRIBUTION:   
    /* print distributions on row, methods in columns */
    n_row = n_distr;
    rltype = LABEL_DISTR;
    n_col = n_meth;
    cltype = LABEL_METH;
    break;
  case ROW_METHOD:
  default:
    /* print distributions on columns, methods in rows*/
    n_row = n_meth;
    rltype = LABEL_METH;
    n_col = n_distr;
    cltype = LABEL_DISTR;
    break;
  }

  printf("Average generation times (including setup) for sample of size %d.\n",samplesize);
  printf("Timings are relative to generation of expontential random variate\n");
  printf("using inversion within UNU.RAN environment\n");
  printf("(timing unit = %g microseconds)\n\n",get_timing_unit());

  /* print table header */
  printf("    ");
  for (row=0; row<n_row; row++) {
    printf("      ");
    print_label(row,rltype);
  }
  printf("\n");

  /* print table rows */ 
  for (col=0; col<n_col; col++) {
    print_label(col,cltype);
    for (row=0; row<n_row; row++) {
      switch (rowentry) {
      case ROW_DISTRIBUTION:   
	/* print distributions on row, methods in columns */
	idx_time = row*n_col+col;
	break;
      case ROW_METHOD:
      default:
	/* print distributions on columns, methods in rows*/
	idx_time = col*n_row+row;
	break;
      }
      if (timings[idx_time] > 0.) {
	printf("%10.4g",timings[idx_time]);
      }
      else {
	printf("%10s","--");
      }
    }
    printf("\n");
  }
  printf("\n");

  /* o.k. */
  return 1;
} /* end of print_timings() */

/*---------------------------------------------------------------------------*/

void
print_usage(void)
     /* print usage */
{
  fprintf(stderr,"\n%s conffile [samplesize]\n",progname);
  fprintf(stderr,"\n");
  fprintf(stderr,"Compute average generation time (including setup) for a sampling.\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"Arguments:\n");
  fprintf(stderr,"\tconffile   ... file with list of distributions and methods\n");
  fprintf(stderr,"\tsamplesize ... size of sample\n");
  fprintf(stderr,"\n");

} /* end of print_usage() */

/*---------------------------------------------------------------------------*/
