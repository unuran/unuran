#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unuran.h>

/* Maximal length of provides string allowed */
#define MAX_STRLENGTH (1024)
/* Maximal entries of any list allowed */
#define MAX_ELEM (1024)

#define UNKNOWN (-2)
#define UNDEF (-1)

/* types of distributions */
#define CONT  (0)
#define CVEC  (1)
#define DISCR (2)
#define CEMP  (3)
#define CVEMP (4)

/* methods for generating random numbers*/
#define AROU (0)
#define CSTD (1)
#define NINV (2)
#define SROU (3)
#define SSR  (4)
#define TABL (5)
#define TDR  (6)
#define UTDR (7)
#define DARI (8)
#define DAU  (9)
#define DGT  (10)
#define DSTD (11)


/* function prototypes */

void         subst_whitespace(char *, char);
char       * elim_whitespace(char *);
int          parselist(char *, double *);
UNUR_DISTR * make_distr_obj(char *);
UNUR_PAR   * make_par_obj(UNUR_DISTR *, char *);
UNUR_GEN   * make_gen_obj (char *);






/**********************************************************************/
/*                                                                    */
/* function: parselist()                                              */
/*                                                                    */
/* called by: make_distr_obj()                                        */
/*            make_par_obj()                                          */
/*                                                                    */
/* gets a string containing a list of comma separated doubles         */
/* which is terminited via a ')' and a pointer to a double array      */
/* and extracts the doubles from the list to the array                */
/*                                                                    */
/* !!!ATTENTION!!!                                                    */
/* there is no white space allowed within the list                    */
/*                                                                    */
/**********************************************************************/
int parselist( char *liststr, double *ptr_to_list ){

  int no_of_elem = 0; /* number of elements in list */

  /* extract doubles from string and write them to array of doubles */
  while ( *liststr != ')' ){   /* end of list is indicated by right bracket */
    /* there can only be ',' or '(' in front of the numbers */
    while ( *liststr == ',' || *liststr == '('){
      liststr++;               /* next char */
    }
    
    /* extract double and write to array */
    ptr_to_list[no_of_elem] = strtod(liststr, &liststr);
    no_of_elem++; /* goto next entry of array */

  } /* end while -- all elements of list read*/

  return (no_of_elem);

} /* end of parselist() */




/**********************************************************************/
/*                                                                    */
/* function: make_distr_obj()                                         */
/*                                                                    */
/* called by: make_gen_obj()                                          */
/*                                                                    */
/* gets a string with information about the distribution,             */
/* generates the corresponding distribution object and returns it     */
/*                                                                    */
/* The string consists of key=value entries separeted by ';'          */
/*                                                                    */
/* three possibilities:                                               */
/* key=value                                                          */
/* key=(komma separated list of numbers)                              */
/* key=value(komma separated list of numbers)                         */
/*                                                                    */
/**********************************************************************/
UNUR_DISTR *make_distr_obj(char *str){

  UNUR_DISTR *distr;

  double list[MAX_ELEM]; /* value contains a list of numbers */
  int no_of_elem;        /* size of that list */

  int type = UNDEF;  /* type of distribution (now not defined)   */
  char * tmpstr;     /* temporary pointer to string */

  char * token;      /* char pointer to tokenize the string */
  char *key, *value; /* the key and its value */

  /* tokenize the string -- split at ';' */
  for ( token  = strtok(str, ";"); 
        token != NULL;     
        token  = strtok(NULL, ";") ){

    /* determine key and value */
    key = token;
    /* one ore more '=' seperates key from value */
    value = strchr(key, '=');
    for ( *value='\0', value++; *value=='='; value++ ); 


    /* split value into value and list (one of them must exist) */
    no_of_elem = 0;
    /* is a list beginning with '(' in the value included? */
    if (NULL != strchr(value, '(') ){

      tmpstr = strchr(value, '(');
      *tmpstr = '\0'; /* terminate value */
      tmpstr++;       /* points to begin of list (if available) */

      /* extract list of doubles from string  */
      no_of_elem = parselist(tmpstr ,list);

    }

    /* Now: key, value and list are determined               */


    /* ------------------------------------------- */
    /*                                             */
    /* distribution                                */
    /*                                             */
    /* ------------------------------------------- */
    if ( !strcmp( key , "distr") ){

      /****************************/
      /* continuous distributions */
      /****************************/
      if ( !strcmp(value, "beta") ){
	type = CONT;
	distr = unur_distr_beta(list, no_of_elem);
      }
      else if ( !strcmp(value, "cauchy") ){
	type = CONT;
	distr = unur_distr_cauchy(list, no_of_elem);
      }
      else if ( !strcmp(value, "chi") ){
	type = CONT;
	distr = unur_distr_chi(list, no_of_elem);
      }
      else if ( !strcmp(value, "chisquare") ){
	type = CONT;
	distr = unur_distr_chisquare(list, no_of_elem);
      }
      else if ( !strcmp(value, "exponential") ){
	type = CONT;
	distr = unur_distr_exponential(list, no_of_elem);
      }
      else if ( !strcmp(value, "extremeI") ){
	type = CONT;
	distr = unur_distr_extremeI(list, no_of_elem);
      }
      else if ( !strcmp(value, "extremeII") ){
	type = CONT;
	distr = unur_distr_extremeII(list, no_of_elem);
      }
      else if ( !strcmp(value, "gamma") ){
	type = CONT;
	distr = unur_distr_gamma(list, no_of_elem);
      }
      else if ( !strcmp(value, "laplace") ){
	type = CONT;
	distr = unur_distr_laplace(list, no_of_elem);
      }
      else if ( !strcmp(value, "logistic") ){
	type = CONT;
	distr = unur_distr_logistic(list, no_of_elem);
      }
      else if ( !strcmp(value, "lomax") ){
	type = CONT;
	distr = unur_distr_lomax(list, no_of_elem);
      }
      else if ( !strcmp(value, "normal") ){
	type = CONT;
	distr = unur_distr_normal(list, no_of_elem);
      }
      else if ( !strcmp(value, "pareto") ){
	type = CONT;
	distr = unur_distr_pareto(list, no_of_elem);
      }
      else if ( !strcmp(value, "powerexponential") ){
	type = CONT;
	distr = unur_distr_powerexponential(list, no_of_elem);
      }
      else if ( !strcmp(value, "rayleigh") ){
	type = CONT;
	distr = unur_distr_rayleigh(list, no_of_elem);
      }
      else if ( !strcmp(value, "triangular") ){
	type = CONT;
	distr = unur_distr_triangular(list, no_of_elem);
      }
      else if ( !strcmp(value, "uniform") ){
	type = CONT;
	distr = unur_distr_uniform(list, no_of_elem);
      }
      else if ( !strcmp(value, "weibull") ){
	type = CONT;
	distr = unur_distr_weibull(list, no_of_elem);
      }
      /**************************/
      /* discrete distributions */
      /**************************/
      else if ( !strcmp(value, "binomial") ){
	type = DISCR;
	distr = unur_distr_binomial(list, no_of_elem);
      }
      else if ( !strcmp(value, "geometric") ){
	type = DISCR;
	distr = unur_distr_geometric(list, no_of_elem);
      }
      else if ( !strcmp(value, "hypergeometric") ){
	type = DISCR;
	distr = unur_distr_hypergeometric(list, no_of_elem);
      }
      else if ( !strcmp(value, "logarithmic") ){
	type = DISCR;
	distr = unur_distr_logarithmic(list, no_of_elem);
      }
      else if ( !strcmp(value, "negativebinomial") ){
	type = DISCR;
	distr = unur_distr_negativebinomial(list, no_of_elem);
      }
      else if ( !strcmp(value, "poisson") ){
	type = DISCR;
	distr = unur_distr_poisson(list, no_of_elem);
      }
      else{
	printf("Unknown distribution: %s\n", value);
	break;
      }
    }

    /* ------------------------------------------- */
    /*                                             */
    /* domain                                      */
    /*                                             */
    /* ------------------------------------------- */
    else if ( !strcmp( key , "domain") ){
      if (type == CONT){
	unur_distr_cont_set_domain( distr, list[0], list[1]);
      }
      else if (type == DISCR){
	unur_distr_discr_set_domain( distr, list[0], list[1]);
      }
      else{
	fprintf(stderr, "Unkown type of distribution while parsing domain!\n");
	break;
      }
    }
    else {
      fprintf(stderr, "Unknown key: %s\n", key);
    }


  } /* end while -- all tokens handled */


  return ( distr );

}
/* end of make_distr_obj()                                            */



/**********************************************************************/
/*                                                                    */
/* function: make_par_obj()                                           */
/*                                                                    */
/* called by: make_gen_obj()                                          */
/*                                                                    */
/* gets a pointer to a distribution object and a string with          */
/* information about the distribution,                                */
/* generates the corresponding parameter object and returns it        */
/*                                                                    */
/* The string consists of key=value entries separeted by ';'          */
/*                                                                    */
/* three possibilities:                                               */
/* key=value                                                          */
/* key=(komma separated list of numbers)                              */
/* key=value(komma separated list of numbers)                         */
/*                                                                    */
/**********************************************************************/
UNUR_PAR *make_par_obj(UNUR_DISTR *distr, char *methstr){

  UNUR_PAR *par;


  int method = UNDEF; /* method to generate random numbers */

  char * token;  /* char pointer to tokenize the string */
  char *key, *value; /* the key and its value */
  char *tmpstr;  /* temporary pointer to char */

  double dblvalue; /* double representation of value */

  double list[MAX_ELEM]; /* value contains a list of double numbers */
  int no_of_elem;        /* size of that list */


  /* tokenize the string */
  for ( token = strtok( methstr, ";");
       token != NULL;
       token = strtok(NULL, ";") ){


    /* determine key and value */
    key = token;
    /* one ore more '=' seperates key from value */
    value = strchr(key, '=');
    for ( *value='\0', value++; *value=='='; value++ ); 


    /* split key into key and list (if possible) */
    no_of_elem = 0;
    if (NULL != strchr(value, '(') ){

      tmpstr = strchr(value, '(');
      *tmpstr = '\0'; /* terminate value */
      tmpstr++;

      /* extract list of doubles from string  */
      no_of_elem = parselist(tmpstr ,list);

    }

    /* extract double value from string value */
    dblvalue = atof(value);

    /* Now: key, value, dblvalue and list are set           */


    /* ************************** */
    /* determine choosen method   */
    /* ************************** */
    if ( !strcmp( key , "method") ){
      /* continuous distribution */
      if ( !strcmp( value , "arou") ){
	method = AROU;
	par = unur_arou_new(distr);
      }
      else if ( !strcmp( value , "cstd") ){
	method = CSTD;
	par = unur_cstd_new(distr);
      }
      else if ( !strcmp( value , "ninv") ){
	method = NINV;
	par = unur_ninv_new(distr);
      }
      else if ( !strcmp( value , "srou") ){
	method = SROU;
	par = unur_srou_new(distr);
      }
      else if ( !strcmp( value , "ssr") ){
	method = SSR;
	par = unur_ssr_new(distr);
      }
      else if ( !strcmp( value , "tabl") ){
	method = TABL;
	par = unur_tabl_new(distr);
      }
      else if ( !strcmp( value , "utdr") ){
	method = UTDR;
	par = unur_utdr_new(distr);
      }
      /* discrete distribution */
      else if ( !strcmp( value , "dari") ){
	method = DARI;
	par = unur_dari_new(distr);
      }
      else if ( !strcmp( value , "dau") ){
	method = DAU;
	par = unur_dau_new(distr);
      }
      else if ( !strcmp( value , "dgt") ){
	method = DGT;
	par = unur_dgt_new(distr);
      }
      else if ( !strcmp( value , "dstd") ){
	method = DSTD;
	par = unur_dstd_new(distr);
      }
      else {
	method = UNKNOWN;
	fprintf(stderr, "Please report: Unknown method!\n");
	break;
      }
    }
    /* ****************************************** */
    /*                                            */
    /* set parameters depending on choosen method */
    /*                                            */
    /* ****************************************** */

    /* ****************************************** */
    /* AROU                                       */
    /* ****************************************** */
    if ( method == AROU && strcmp(key, "method") ){
      if ( !strcmp(key, "max_sqhratio") ){
	unur_arou_set_max_sqhratio(par, dblvalue);
      }
      else if ( !strcmp(key, "max_segments") ){
	unur_arou_set_max_segments(par, (int) dblvalue);
      }
      else if ( !strcmp(key, "cpoints") ){
	unur_arou_set_cpoints(par,  no_of_elem, list);
      }
      else if ( !strcmp(key, "center") ){
	unur_arou_set_center(par, dblvalue);
      }
      else if ( !strcmp(key, "usecenter") ){
	unur_arou_set_usecenter(par, (int) dblvalue);
      }
      else if ( !strcmp(key, "guidefactor") ){
	unur_arou_set_guidefactor(par, dblvalue);
      }
      else if ( !strcmp(key, "verify") ){
	unur_arou_set_verify(par, (int) dblvalue);
      }
      else if ( !strcmp(key, "pedantic") ){
	unur_arou_set_pedantic(par, (int) dblvalue);
      }
      else {
	fprintf (stderr, "Unknown option for method AROU: %s\n", key);
      }
    }
    /* ****************************************** */
    /* CSTD                                       */
    /* ****************************************** */
    else if ( method == CSTD && strcmp(key, "method") ){
      if ( !strcmp(key, "variant") ){
	unur_cstd_set_variant(par, (unsigned) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method CSTD: %s\n", key);
      }
    }
    /* ****************************************** */
    /* NINV                                       */
    /* ****************************************** */
    else if ( method == NINV && strcmp(key, "method") ){
      if ( !strcmp(key, "usenewton") ){
	unur_ninv_set_usenewton(par);
      }
      else if ( !strcmp(key, "useregula") ){
	unur_ninv_set_useregula(par);
      }
      else if ( !strcmp(key, "max_iter") ){
	unur_ninv_set_max_iter(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "x_resolution") ){
	unur_ninv_set_x_resolution(par, dblvalue );
      }
      else if ( !strcmp(key, "start") ){
	unur_ninv_set_start(par, list[0], list[1] );
      }
      else if ( !strcmp(key, "table") ){
	unur_ninv_set_table(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method NINV: %s\n", key);
      }
    }
    /* ****************************************** */
    /* SROU                                       */
    /* ****************************************** */
    else if ( method == SROU && strcmp(key, "method") ){
      if ( !strcmp(key, "cdfatmode") ){
	unur_srou_set_cdfatmode(par, dblvalue );
      }
      else if ( !strcmp(key, "pdfatmode") ){
	unur_srou_set_pdfatmode(par, dblvalue );
      }
      else if ( !strcmp(key, "usesqueeze") ){
	unur_srou_set_usesqueeze(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "usemirror") ){
	unur_srou_set_usemirror(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "verify") ){
	unur_srou_set_verify(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method SROU: %s\n", key);
      }
    }
    /* ****************************************** */
    /* SSR                                       */
    /* ****************************************** */
    else if ( method == SSR && strcmp(key, "method") ){
      if ( !strcmp(key, "cdfatmode") ){
	unur_ssr_set_cdfatmode(par, dblvalue );
      }
      else if ( !strcmp(key, "pdfatmode") ){
	unur_ssr_set_pdfatmode(par, dblvalue  );
      }
      else if ( !strcmp(key, "usesqueeze") ){
	unur_ssr_set_usesqueeze(par,  (int) dblvalue );
      }
      else if ( !strcmp(key, "verify") ){
	unur_ssr_set_verify(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method SSR: %s\n", key);
      }
    }
    /* ****************************************** */
    /* TABL                                       */
    /* ****************************************** */
    else if ( method == TABL && strcmp(key, "method") ){
      if ( !strcmp(key, "variant_setup") ){
	unur_tabl_set_variant_setup(par, (unsigned) dblvalue  );
      }
      else if ( !strcmp(key, "variant_splitmode") ){
	unur_tabl_set_variant_splitmode(par, (unsigned) dblvalue  );
      }
      else if ( !strcmp(key, "max_sqhratio") ){
	unur_tabl_set_max_sqhratio(par, dblvalue  );
      }
      else if ( !strcmp(key, "max_intervals") ){
	unur_tabl_set_max_intervals(par, (int) dblvalue  );
      }
      else if ( !strcmp(key, "areafraction") ){
	unur_tabl_set_areafraction(par, dblvalue );
      }
      else if ( !strcmp(key, "nstp") ){
	unur_tabl_set_nstp(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "slopes") ){
	fprintf(stderr, "No slopes support from string interface!");
      }
      else if ( !strcmp(key, "guidefactor") ){
	unur_tabl_set_guidefactor(par, dblvalue  );
      }
      else if ( !strcmp(key, "boundary") ){
	unur_tabl_set_boundary(par, list[0], list[1]  );
      }
      else if ( !strcmp(key, "verify") ){
	unur_tabl_set_verify(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method TABL: %s\n", key);
      }
    }
    /* ****************************************** */
    /* TDR                                        */
    /* ****************************************** */
    else if ( method == TDR && strcmp(key, "method") ){
      if ( !strcmp(key, "c") ){
	unur_tdr_set_c(par, dblvalue );
      }
      else if ( !strcmp(key, "variant_gw") ){
	unur_tdr_set_variant_gw(par);
      }
      else if ( !strcmp(key, "variant_ps") ){
	unur_tdr_set_variant_ps(par);
      }
      else if ( !strcmp(key, "variant_ia") ){
	unur_tdr_set_variant_ia(par);
      }
      else if ( !strcmp(key, "max_sqhratio") ){
	unur_tdr_set_max_sqhratio(par,  dblvalue);
      }
      else if ( !strcmp(key, "max_intervals") ){
	unur_tdr_set_max_intervals(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "cpoints") ){
	unur_tdr_set_cpoints(par, no_of_elem, list);
      }
      else if ( !strcmp(key, "center") ){
	unur_tdr_set_center(par, dblvalue );
      }
      else if ( !strcmp(key, "usecenter") ){
	unur_tdr_set_usecenter(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "usemode") ){
	unur_tdr_set_usemode(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "guidefactor") ){
	unur_tdr_set_guidefactor(par, dblvalue );
      }
      else if ( !strcmp(key, "verify") ){
	unur_tdr_set_verify(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "pedantic") ){
	unur_tdr_set_pedantic(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method TDR: %s\n", key);
      }
    }
    /* ****************************************** */
    /* UTDR                                       */
    /* ****************************************** */
    else if ( method == UTDR && strcmp(key, "method") ){
      if ( !strcmp(key, "pdfatmode") ){
	unur_utdr_set_pdfatmode(par, dblvalue );
      }
      else if ( !strcmp(key, "cpfactor") ){
	unur_utdr_set_cpfactor(par, dblvalue );
      }
      else if ( !strcmp(key, "deltafactor") ){
	unur_utdr_set_deltafactor(par, dblvalue );
      }
      else if ( !strcmp(key, "verify") ){
	unur_utdr_set_verify(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method UTDR: %s\n", key);
      }
    }
    /* ****************************************** */
    /* DARI                                       */
    /* ****************************************** */
    else if ( method == DARI && strcmp(key, "method") ){
      if ( !strcmp(key, "squeeze") ){
	unur_dari_set_squeeze(par, (char) dblvalue );
      }
      else if ( !strcmp(key, "tablesize") ){
	unur_dari_set_tablesize(par, (int) dblvalue );
      }
      else if ( !strcmp(key, "cpfactor") ){
	unur_dari_set_cpfactor(par, dblvalue );
      }
      else if ( !strcmp(key, "verify") ){
	unur_dari_set_verify(par, (int) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method DARI: %s\n", key);
      }
    }
    /* ****************************************** */
    /* DAU                                        */
    /* ****************************************** */
    else if ( method == DAU && strcmp(key, "method") ){
      if ( !strcmp(key, "urnfactor") ){
	unur_dau_set_urnfactor(par, dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method DAU: %s\n", key);
      }
    }
    /* ****************************************** */
    /* DGT                                        */
    /* ****************************************** */
    else if ( method == DGT && strcmp(key, "method") ){
      if ( !strcmp(key, "guidefactor") ){
	unur_dgt_set_guidefactor(par, dblvalue );
      }
      else if ( !strcmp(key, "variant") ){
	unur_dgt_set_variant(par,  (unsigned) dblvalue );
      }
      else {
	fprintf (stderr, "Unknown option for method DGT: %s\n",key);
      }
    }
    /* ****************************************** */
    /* DSTD                                       */
    /* ****************************************** */
    else if ( method == DSTD && strcmp(key, "method") ){
      if ( !strcmp(key, "variant") ){
	unur_dstd_set_variant(par, (unsigned) dblvalue  );
      }
      else {
	fprintf (stderr, "Unknown option for method DSTD: %s\n", key);
      }
    }
    /* no method defined -> error */
    else if ( !strcmp(key, "method") && method == UNDEF){
      fprintf(stderr, "Error: No method defined!\n");
    }

  } /* end while -- all tokens handled */


  /* return the parameter object */
  return ( par );

}
/* end of make_par_obj()                                              */


/**********************************************************************/
/*                                                                    */
/* function: subst_whitespce()                                        */
/*                                                                    */
/* finds all occurences of sep_char in str and replaces the           */
/* surrounding whitespace with the separation character sep_char      */
/*                                                                    */
/**********************************************************************/
void subst_whitespace(char *str, char sep_char){

  char *tmpstr;

  for (tmpstr = str; tmpstr != NULL; tmpstr =strchr(tmpstr, sep_char) ){

    for ( tmpstr--; isspace(*tmpstr); tmpstr-- )
      *tmpstr = sep_char;

    for ( tmpstr++; isspace(*tmpstr) || (*tmpstr == sep_char); tmpstr++ )
      *tmpstr = sep_char;
  }

}
/* end of subst_whitespace()                                          */


/**********************************************************************/
/*                                                                    */
/* function: elim_whitespce()                                         */
/*                                                                    */
/* gets a string, cuts off leading and terminating white space and    */
/* replaces all other white space within the string with the nearby   */
/* separation charakter                                               */
/*                                                                    */
/* returns pointer to the string                                      */
/*                                                                    */
/**********************************************************************/
char *elim_whitespace(char *str){

  char * tmpstr;

  /* remove initial white space */
  while (*str == ' ') str++;

  /* remove terminating white space */
  tmpstr = strchr(str, '\0');
  tmpstr--;
  while ( isspace(*tmpstr) ){
    *tmpstr = '\0';
    tmpstr--;
  }


  /* substitute whitespace whith nearby characters*/
  subst_whitespace(str, ':');
  subst_whitespace(str, ';');
  subst_whitespace(str, ',');
  subst_whitespace(str, '=');
  subst_whitespace(str, ')');
  subst_whitespace(str, '(');

  return ( str );
}
/* end of elim_whitespace()                                           */





/**********************************************************************/
/*                                                                    */
/* function: make_gen_obj                                             */
/*                                                                    */
/* gets a string with information about the distribution and          */
/* the desired method and                                             */
/* generates the corresponding generator object invoking              */
/* the funktions                                                      */
/*   make_dist_obj() and                                              */
/*   make_par_obj()                                                   */
/*                                                                    */
/* The string consists of two substrings (separated by ':')           */
/* (distributio info, method info )                                   */
/* and each of this strings consists of key=value                     */
/*  entries separeted by ';'                                          */
/*                                                                    */
/* three possibilities:                                               */
/* key=value                                                          */
/* key=(komma separated list of numbers)                              */
/* key=value(komma separated list of numbers)                         */
/*                                                                    */
/**********************************************************************/
UNUR_GEN *make_gen_obj (char *str){

  UNUR_DISTR *distr;       /* distribution object */
  UNUR_PAR *par;           /* parameter object */
  UNUR_GEN *gen;           /* generator object    */

  char *diststr; /* string describing distribution */
  char *methstr; /* string describing method       */
  char *urngstr; /* string holding info for uniform generator */
  char *tmpchar; /* temporary pointer to char      */


  /* convert string to lowercase */
  tmpchar = str;
  while ( *tmpchar != '\0' ){
    *tmpchar = tolower(*tmpchar);
    tmpchar++;
  }

  /* remove initial and terminating white space and
     substitute other whitespace within string with nearby
     separation characters                                 */
  tmpchar = elim_whitespace( str );

  /* split info about distribution, method und uniform generator */
  diststr = strtok(tmpchar, ":");
  methstr = strtok(NULL, ":");
  urngstr = strtok(NULL, ":");


  /* generate distribution object */
  distr = make_distr_obj(diststr);

  /* generate parameter object */
  if ( methstr != NULL ){ /* method info is provided */
    par = make_par_obj(distr, methstr);
  }
  else{ /* no info about method provided -> standard method */
    par = unur_cstd_new(distr);
  }

  /* generate generator object */
  gen = unur_init(par);

  /* return generator object to calling routine */
  return (gen);

}





/* main -- test program */
int main(){

  int i;

  UNUR_GEN *gen;
  char str[] = "    distr  =normal   ((((1,  1  )  ; domain  =  (  -1,  1) : method= arou   (1,2)  :   weiter =(2,2)  ";
  gen = make_gen_obj(str);

  for ( i=0; i<15; i++){
    printf("random number: %f\n", unur_sample_cont(gen));
  }

  return (0);
}
