
#include <source_unuran.h>

char            *readln  (char *s);
void  show_symb_tab      (void);


#define MAXLENGTH 256

int main()
{
  struct treenode *parsetree;
/*    struct treenode *dev_tree; */
/*    char            input_string[MAXLENGTH]; */
/*    int             errcode, errpos; */
/*    double          argument; */
  double          value;

  char fstr[MAXLENGTH] = "f(x) = -exp(x-3.) +  ln(2*3-5) + pi + x + 3 * x + pi + 1 + 3.4 * x^2 - 5 * x^4 + 365.4 + ln(x) + not(x>9) + log(2,exp(x)) + 1.5e-1-1.098612 + (x < 10^2)    ";




/*    _unur_fstr_init(); */
  
  /* Einlesen einer Funktion als string */
/*    printf("\n\nFunktion eingeben:\n"); */
  /*       readln(input_string); */
  /*       if (strcmp(input_string,"") == 0) break; */
  /*       parsetree = _unur_fstr2tree(input_string,&errcode,&errpos); */

  parsetree = _unur_fstr2tree(fstr);
    
  if  (parsetree == NULL)  {  
    printf("Fehler!\n");
/*          exit (EXIT_FAILURE); */
  };

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

  /* if (value==nan) {
     printf("Fehler\n");
     break; }
  */     

  /*-----------------------------------------------------------------*/
  /* Stringausgabe    */

/*    Ntree2string(parsetree,input_string); */
/*    printf("\nParse-Baum als String:\n%s\n",input_string);  */

  /*-----------------------------------------------------------------*/
  /*  Ableitung */

  /*       readln(input_string); */
    
  /*     show_symb_tab(); */

/*    printf("\n Wert: %f \n",value); */

   
  /* Speicher fuer tree freigeben */
/*    _unur_fstr_free(parsetree); */

  
  exit(EXIT_SUCCESS);
}  



































