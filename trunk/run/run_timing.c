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

char *_unur_parser_prepare_string( const char *str );

/*---------------------------------------------------------------------------*/

int main()
{

#define LINELENGTH  1024      /* max length of lines allowed    */

  char line[LINELENGTH];      /* input buffer */

  char *str;

  struct unur_slist *distr_str_list; /* list of strings for distributions */
  struct unur_slist *meth_str_list;  /* list of strings for methods */

  struct unur_slist *distr_list; /* list of distributions */
  struct unur_slist *meth_list;   /* list of methods */

  int n_distr;        /* number of distributions */
  int n_meth;         /* number of methods */

  UNUR_DISTR *distr;
  UNUR_PAR *par;
  UNUR_GEN *gen;
  int i;

  FILE *fhconf;   /* file handler for config file */

  char fileconf[] = "t.conf";

  /* make lists for distributions and methods */
  distr_str_list = _unur_slist_new();
  meth_str_list  = _unur_slist_new();
  distr_list = _unur_slist_new();
  meth_list  = _unur_slist_new();

  /* open the file with the data */
  fhconf = fopen(fileconf, "r");
  if (fhconf == NULL) {
    fprintf(stderr,"error\n");
    exit(EXIT_FAILURE);
  }


  /* read lines until eof */
  while( fgets(line, LINELENGTH, fhconf) ) {

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

    fprintf(stderr,"syntax error: %s\n",str);
  }

  /* close input stream */
  fclose(fhconf);

  /* get number of distributions and methods */
  n_distr = _unur_slist_length(distr_str_list);
  n_meth  = _unur_slist_length(meth_str_list);

  fprintf(stderr,"n_distr = %d, n_meth = %d\n",n_distr,n_meth);
  
  /* get all distribution objects */
  for (i=0; i<n_distr; i++) {
    str = _unur_slist_get(distr_str_list,i);
    distr = unur_str2distr(str);
    if (distr == NULL)
      fprintf(stderr,"syntax error: %s\n",str);
    else
      _unur_slist_append(distr_list,distr);
  }

  /* update number of distributions */
  n_distr = _unur_slist_length(distr_list);

  fprintf(stderr,"n_distr = %d, n_meth = %d\n",n_distr,n_meth);

  /* free memory */
  _unur_slist_free(distr_str_list);
  _unur_slist_free(meth_str_list);
  _unur_slist_free(distr_list);
  _unur_slist_free(meth_list);

  exit (EXIT_SUCCESS);

#undef LINELENGTH

} /* end of main() */

/*---------------------------------------------------------------------------*/













