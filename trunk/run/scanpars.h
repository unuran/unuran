
#include <stdio.h> 
#include <stdlib.h>
#include <string.h> 


#define SYMBLENGTH 20 
#define MAXLENGTH  256 
#define NOSYMB  0 
#define ADD_OP  1 
#define REL_OP  2 
#define MUL_OP  3 
#define HPR_OP  4 
#define OTHERS  5 
#define SCONST  6 
#define UCONST  7 
#define UIDENT  8 
#define UFUNCS  9 
#define SFUNCS 10 

#ifndef FALSE
#define FALSE  (0)                   /* boolean false                   */
#define TRUE   (!FALSE)              /* boolean true                    */
#endif

struct treenode { 
  char            symb[SYMBLENGTH];  /* zeigt auf Symbol aus Symboltab. */ 
  int             token;             /* Token des Symbols               */ 
  int             symbkind;          /* Art des Symbols (REL_OP etc.)   */ 
  float           val;               /* aktueller arithmetischer Wert   */ 
  struct treenode *left;             /* Zeiger auf linken Sohn          */ 
  struct treenode *right;            /* Zeiger auf rechten Sohn         */ 
}; 

struct symbols { 
  char            name[SYMBLENGTH];  /* Name des Symbols (z. B. "SIN")  */ 
  int             info;              /* Prioritaet bzw. Argumentanzahl  */ 
  double           val;               /* Konstanten: numerischer Wert    */ 
  double           (*vcalc)(int t, double l, double r);        
                                     /* Zeiger auf Berechnungsfunktion  */ 
  char            *(*dcalc)(char *par,struct treenode *w,
                            char *l, char *r, char *dl, char *dr,char *s);
                                     /* Zeiger auf Ableitungsfunktion   */ 
  struct treenode *tree;             /* Bei UFUNCS: Zeiger auf Baum     */ 
};



















