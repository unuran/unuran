
#include <source_unuran.h>

char            *readln  (char *s);
void  show_symb_tab      (void);


#define MAXLENGTH 256

int main()
{
  struct ftreenode *parsetree;
  struct ftreenode *dev_tree;
/*    char            input_string[MAXLENGTH]; */
/*    int             errcode, errpos; */
/*    double          argument; */
  double          value;

#if 0
  char fstr[MAXLENGTH] = "f(x) = -exp(x-3.) +  ln(2*3-5) + pi + x + 3 * x + pi + 1 + 3.4 * x^2 - (5 * x^4 + 365.4 + ln(x) + not(x>9) + log(2,exp(x)) + 1.5e-1-1.098612) + (x < 10^2)    ";

  parsetree = _unur_fstr2tree_DefFunct(fstr);

#else

  char fstr[MAXLENGTH] = "-exp(x-3.) +  ln(2*3-5) + pi + x + 3 * x + pi + 1 + 3.4 * x^2 - (5 * x^4 + 365.4 + ln(x) + not(x>9) + log(2,exp(x)) + 1.5e-1-1.098612) + (x < 10^2)    ";

  parsetree = _unur_fstr2tree(fstr);

#endif
    
  if  (parsetree == NULL)  {  
    printf("Fehler!\n");
        exit (EXIT_FAILURE);
  };


  printf("%s\n",_unur_fstr_tree2string(parsetree,"y","F"));
  


/*      _unur_fstr_debug_tree(parsetree); */
  fflush(stdout);

  /*-----------------------------------------------------------------*/
  /* Funktionsauswertung */

/*    printf("\n Argument: "); */
  /*       readln(input_string); */
  /*       argument=atof(input_string); */
  /*       value=_unur_fstr_eval_tree(parsetree,atof(input_string)); */

  value = _unur_fstr_eval_tree(parsetree,3.);
  printf("\n Wert: %f \n",value);


  dev_tree = _unur_fstr_make_derivative(parsetree);

   
  /* Speicher fuer tree freigeben */
  _unur_fstr_free(parsetree);
  _unur_fstr_free(dev_tree);

  
  exit(EXIT_SUCCESS);
}  



































