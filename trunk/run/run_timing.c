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

#include <ctype.h>
#include <malloc.h>
#include <string.h>

#include <unuran.h>
#include <unuran_tests.h>

/*---------------------------------------------------------------------------*/

extern char *_unur_parser_prepare_string( const char *str );

/*---------------------------------------------------------------------------*/

/* read config file */
int read_config_file ( const char *filename, 
		       struct unur_slist *distr_str_list, struct unur_slist *meth_str_list ); 

/* make distribution objects */
struct unur_slist *make_distr_list ( struct unur_slist *distr_str_list );

/* print legend for distributions and methods */
int print_legend ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list ); 

/* print name of distribution with index n*/
int print_name_distr ( struct unur_slist *distr_str_list, int n );

/* print name of method with index n*/
int print_name_meth ( struct unur_slist *meth_str_list, int n );

/* compute timings table */
double *compute_timings ( struct unur_slist *distr_str_list, struct unur_slist *meth_str_list,
			  int samplesize, double duration ); 

/*---------------------------------------------------------------------------*/

int main()
{
  struct unur_slist *distr_str_list; /* list of strings for distributions */
  struct unur_slist *meth_str_list;  /* list of strings for methods */

  double *time_0;

  char fileconf[] = "t.conf";

  /* create lists for distributions and methods */
  distr_str_list = _unur_slist_new();
  meth_str_list  = _unur_slist_new();

  /* read config file */
  read_config_file(fileconf, distr_str_list, meth_str_list);

  /* print legend */
  print_legend(distr_str_list,meth_str_list);
  
  /* make timings */
  time_0 = compute_timings(distr_str_list,meth_str_list,10,0.1);

  /* free memory */
  _unur_slist_free(distr_str_list);
  _unur_slist_free(meth_str_list);

  exit (EXIT_SUCCESS);

} /* end of main() */

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
    print_name_distr(distr_str_list,i);
    printf(" ... %s\n", str);
  }
  printf("\n");

  printf("%d methods:\n",n_meth);
  for (k=0; k<n_meth; k++) {
    str = _unur_slist_get(meth_str_list, k);
    print_name_meth(meth_str_list,k);
    printf(" ... %s\n", str);
  }
  printf("\n");

  /* o.k. */
  return 1;
} /* end of print_legend() */

/*---------------------------------------------------------------------------*/

int
print_name_distr ( struct unur_slist *distr_str_list, int n )
     /* print name of distribution with index n */
{
  printf("[%c]", 'A'+n);
  return 1;
} /* end of print_name_distr() */

/*---------------------------------------------------------------------------*/

int
print_name_meth ( struct unur_slist *meth_str_list, int n )
     /* print name of method with index n*/
{
  printf("[%d]", n);
  return 1;
} /* end of print_name_meth() */

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
	timing[i*n_distr+k] = unur_test_timing_total(par, samplesize, duration );
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
