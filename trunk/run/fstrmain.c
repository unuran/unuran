#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "pars.h"

char            *readln  (char *s);


#define MAXLENGTH 50

int main()
{
  struct treenode *parsetree,*dev_tree;
  char            *input_string;
  int             errcode, errpos;
  double          argument;

   _unur_fstr_init();
  

do {
  /* Einlesen einer Funktion als string */
     printf("\n\nFunktion eingeben:\n");
     input_string= (char *)malloc(MAXLENGTH*sizeof(char));
     readln(input_string);
     if (strcmp(input_string,"") == 0) break;
     parsetree = _unur_fstr2tree(input_string,&errcode,&errpos);
    

     if  (errcode>0)  {  
              printf("Fehler\n");
              break;
            };
  /*-----------------------------------------------------------------*/
  /* Funktionsauswertung */

   do {
     printf("\n Argument: ");
     readln(input_string);
     argument=atof(input_string); 
     printf("\n Wert: %f \n", _unur_fstr_eval_tree(parsetree,atof(input_string)));
    } while (0); 

/*-----------------------------------------------------------------*/
 /* Stringausgabe    */

 /*    Ntree2string(parsetree,input_string);
       printf("\nParse-Baum als String:\n%s\n",input_string); */

 /*-----------------------------------------------------------------*/
 /*  Ableitung */

   do {
     printf("\nArgument fuer Ableitung:\n");readln(input_string);
     dev_tree=_unur_fstr_make_derivative(parsetree);
    
     printf("\n Wert: %f \n", _unur_fstr_dev_eval_tree(dev_tree,atof(input_string)));
    } while (0); 
  
   
   show_tree(parsetree);
   
   /* Speicher fuer tree freigeben */
     _unur_fstr_free(parsetree);
     _unur_fstr_free(dev_tree);
     free (input_string);
 
     
  
  

}    while (0);
 return 0;


}  



































