/************************************************************************
 * MODUL pars.c                                                         *
 *                                                                      *
 *                                                                      * 
 ************************************************************************/

#include "pars.h" 

#include <math.h>
#include <time.h>
#include <stdlib.h>

/*long _STKSIZ = 500000; grosser Stack f. Rekursion in derive_expression*/

/************************************************************************/

/* --- Importe aus scanner.c und parser.c: --- */

extern struct symbols symbol[];
extern int            ros, roe, aos, aoe, mos, moe, hos, hoe, oss, ose; 
extern int            scs, sce, ucs, uce, uis, uie, ufs, ufe, sfs, sfe; 
extern char           *errorstrings[]; 
extern struct treenode *string2tree(char *function, int *errcodep, 
                                                    int *errposp);
extern int  nxt_symbol  (char function[], int *scanposp,  char symb[],
                             int *tokenp, int *symbkindp, int *errposp);
extern void init_scanner(void), show_symb_tab(void); 
extern void pt_error    (char *function, int errcode, int errpos);
extern void show_tree   (struct treenode *root);
extern void _unur_fstr_free (struct treenode *root);                                
extern char *readln     (char *s);

/************************************************************************/ 

/* --- Prototypen: --- */

char *tree2string(struct treenode *tree_root, char *ret_str);
char *tree2Cstring(struct treenode *tree_root, char *ret_str);
char *strcpy3    (char *s1, char *s2, char *s3, char *s4); 
char *strcpy5    (char *s1, char *s2, char *s3, char *s4, char *s5, char *s6); 
double tree2float (struct treenode *E_root); 
double v_dummy  (int t, double l, double r);
double v_less   (int t, double l, double r);
double v_equal  (int t, double l, double r);
double v_greater(int t, double l, double r);
double v_less_or(int t, double l, double r);
double v_unequal(int t, double l, double r);
double v_grtr_or(int t, double l, double r);
double v_or     (int t, double l, double r);
double v_xor    (int t, double l, double r);
double v_plus   (int t, double l, double r); 
double v_minus  (int t, double l, double r); 
double v_mul    (int t, double l, double r);
double v_and    (int t, double l, double r);
double v_div    (int t, double l, double r);
double v_mod    (int t, double l, double r);
double v_power  (int t, double l, double r);
double v_not    (int t, double l, double r);
double v_sconst (int t, double l, double r);
double v_uconst (int t, double l, double r);
double v_uident (int t, double l, double r);
double v_ufuncs (int t, double l, double r);
double v_exp    (int t, double l, double r);
double v_ln     (int t, double l, double r);
double v_log    (int t, double l, double r);
double v_sin    (int t, double l, double r); 
double v_cos    (int t, double l, double r);
double v_tan    (int t, double l, double r);
double v_sec    (int t, double l, double r); 
double v_sqr    (int t, double l, double r);
double v_abs    (int t, double l, double r);
void gradient(struct treenode *root);
void nxt_part_derivation(struct treenode *DP_root,struct treenode *root);
struct treenode *part_deriv(struct treenode *root,char *par,char *ret_str); 
char *derive_expression(struct treenode *E_root,char *par,char *ret_str);
char *strcpy6(char *s1, char *s2, char *s3,
	      char *s4, char *s5, char *s6, char *s7);
void str_upr(char *s);


/************************************************************************
 * Umwandlung eines Parsebaumes in einen String:                        *
 ************************************************************************/

char *tree2string(struct treenode *tree_root, char *ret_str)

/*  Wandelt einen vom Parser erzeugten Baum zurueck in einen String.
 *  Zunaechst werden der linke und der rechte Zweig des Knotens in einen
 *  String verwandelt; abhaengig vom Operator der Wurzel und den Opera-
 *  toren der Zweige muessen die Zweige ggf. eingeklammert werden, dann
 *  werden die beiden Zweige ueber den Operator der Wurzel miteinander
 *  verknuepft.
 */

{
  struct treenode *right, *left;
  char            symb_m[SYMBLENGTH]; /* Symbol an der Wurzel ("Mitte") */
  char            symb_l[SYMBLENGTH]; /* Symbol des linken Zweiges      */
  char            symb_r[SYMBLENGTH]; /* Symbol des rechten Zweiges     */
  char            str_l[MAXLENGTH];   /* kompl. linker Zweig als String */
  char            str_r[MAXLENGTH];   /* kompl. rechter Zweig           */
  char            temp_str[MAXLENGTH];/* temporaerer String             */
  int             t_m, t_l, t_r;      /* Token Mitte, links, rechts     */
  int             p_m, p_l, p_r;      /* Priorit. Mitte, links, rechts  */
  int             sk_m, sk_l, sk_r;   /* symbkind Mitte,links,rechts    */

  str_l[0] = str_r[0] = symb_l[0] = symb_r[0] = '\0';
  temp_str[0] = ret_str[0] = '\0';

  /* --- Werte der Wurzel: --- */
  strcpy(symb_m,tree_root->symb);
  sk_m    = tree_root->symbkind;
  t_m     = tree_root->token;
  p_m     = symbol[t_m].info;
  left    = tree_root->left;
  right   = tree_root->right;

  /* --- Zuerst die beiden Zweige in Strings verwandeln: --- */
  if( left != NULL ){
     /* --- Werte des linken Sohnes: --- */
     strcpy(symb_l,left->symb);
     sk_l    = left->symbkind;
     t_l     = left->token;
     p_l     = symbol[t_l].info;
     tree2string(left,str_l);
  }
  if( right != NULL ){
     /* --- Werte des rechten Sohnes: --- */
     strcpy(symb_r,right->symb);
     sk_r    = right->symbkind;
     t_r     = right->token;
     p_r     = symbol[t_r].info;
     tree2string(right,str_r);
  }

  /* --- Beide Strings ueber die Wurzel verknuepfen, ggf. Klammern: --- */

  if( t_m <= hos+1 ){
     /* --- Alle Operatoren mit zwei Eingaengen (von '<' bis '^'): --- */

     /* --- Funktionen links bzgl. Prioritaet  wie '*'-Operatoren: --- */
     if (sk_l == SFUNCS || sk_l == UFUNCS) sk_l = MUL_OP;

     /* --- Dummy-Null bei neg. Vorzeichen entfernen: --- */
     if( strcmp(symb_m,"-") == 0 && strcmp(symb_l,"0") == 0 ){
        if( p_r <= p_m && t_r < hoe 
           ){ strcpy3(ret_str,"-(",str_r,")"); 
           }else{ strcpy(ret_str,"-"); strcat(ret_str,str_r); 
	   } 
        return ret_str; 
     } 

     /* --- Linken String evtl. einklammern: --- */ 
     if( t_l < hoe ){ /* alle Operatoren incl. 'NOT' */ 
	if( strcmp(symb_m,"^") == 0
	    ){ /* -- bei '^' auch bei gleicher Prioritaet klammern: --*/
	  if( p_l <= p_m ){
		    strcpy3(temp_str,"(",str_l,")");
		    strcpy(str_l,temp_str);
	  }
	   } else{ /* --- sonst nur bei niedrigerer Prioritaet klammern,
		    aber nicht bei '+' mit '-' [z.B. f(x)=(1+2)-3]: --- */
		 if( (p_l < p_m) &&
		     !(strcmp(symb_l,"+")==0 && strcmp(symb_m,"-")==0) ){
		    strcpy3(temp_str,"(",str_l,")");
		    strcpy(str_l,temp_str);
		 }
	   }
     }

     /* --- Rechten String evtl. einklammern: --- */
     if ( t_r < hoe ){                /*alle Operatoren incl. 'NOT' */
	if( strcmp(symb_m,"/") == 0 || strcmp(symb_m,"MOD") == 0 ||
	    strcmp(symb_m,"-") == 0 || strcmp(symb_m,"^") == 0
	    ){ /* --- klammern auch bei gleicher Prioritaet: --- */
		 if( symbol[t_r].info <= symbol[t_m].info ){
		    strcpy3(temp_str,"(",str_r,")");
		    strcpy(str_r,temp_str);
		 }
	    }else{ if( symbol[t_r].info < symbol[t_m].info ){
		    strcpy3(temp_str,"(",str_r,")");
		    strcpy(str_r,temp_str);
		 }
	}
     }

     /* --- Bei einigen Operatoren Spaces drumherum setzen: --- */
     if( strcmp(symb_m,"XOR") == 0 || strcmp(symb_m,"AND") == 0 ||
	strcmp(symb_m,"OR" ) == 0 || strcmp(symb_m,"MOD") == 0 ){
	 strcpy3(temp_str," ",symb_m," ");
	 strcpy(symb_m,temp_str);
     }

     /* --- Rueckgabestring zusammensetzen aus linkem String,
			  Verknuepfungsoperator und rechtem String: --- */
     return strcpy3(ret_str,str_l,symb_m,str_r);
  }

  /*-- bei "NOT" klammern, wenn Operator geringerer Prioritaet folgt: --*/
  if( strcmp(symb_m,"NOT") == 0 ){
     if( sk_r==REL_OP || sk_r==ADD_OP||sk_r==MUL_OP||strcmp(symb_r,"^")==0
	){ strcpy3(ret_str,"NOT(",str_r,")");
	     return ret_str;
	}else{ /* --- sonst nicht klammern: --- */
	     return strcpy3(ret_str,"NOT ",str_r,"");
     }
  }

  /* --- Parameterlisten von Funktionen: --- */
  if (strcmp(symb_m,",") == 0) return strcpy3(ret_str,str_l,",",str_r);

  /* --- Funktionen: Parameterlisten in Klammern setzen: --- */
  if( sk_m == SFUNCS || sk_m == UFUNCS ){
     strcpy3(ret_str,symb_m,"(",str_r);
     return strcat(ret_str,")");
  }

  /* --- Alle Endsymbole direkt zurueck: --- */
  return strcpy(ret_str,symb_m);
}




/************************************************************************
 * Umwandlung eines Parsebaumes in einen C_code:                        *
 ************************************************************************/

char *tree2Cstring(struct treenode *tree_root, char *ret_str)

/*  Wandelt einen vom Parser erzeugten Baum zurueck in einen String.
 *  Zunaechst werden der linke und der rechte Zweig des Knotens in einen
 *  String verwandelt; abhaengig vom Operator der Wurzel und den Opera-
 *  toren der Zweige muessen die Zweige ggf. eingeklammert werden, dann
 *  werden die beiden Zweige ueber den Operator der Wurzel miteinander
 *  verknuepft.
 */

{
  struct treenode *right, *left;
  char            symb_m[SYMBLENGTH]; /* Symbol an der Wurzel ("Mitte") */
  char            symb_l[SYMBLENGTH]; /* Symbol des linken Zweiges      */
  char            symb_r[SYMBLENGTH]; /* Symbol des rechten Zweiges     */
  char            str_l[MAXLENGTH];   /* kompl. linker Zweig als String */
  char            str_r[MAXLENGTH];   /* kompl. rechter Zweig           */
  char            temp_str[MAXLENGTH];/* temporaerer String             */
  int             t_m, t_l, t_r;      /* Token Mitte, links, rechts     */
  int             p_m, p_l, p_r;      /* Priorit. Mitte, links, rechts  */
  int             sk_m, sk_l, sk_r;   /* symbkind Mitte,links,rechts    */

  str_l[0] = str_r[0] = symb_l[0] = symb_r[0] = '\0';
  temp_str[0] = ret_str[0] = '\0';

  /* --- Werte der Wurzel: --- */
  strcpy(symb_m,tree_root->symb);
  sk_m    = tree_root->symbkind;
  t_m     = tree_root->token;
  p_m     = symbol[t_m].info;
  left    = tree_root->left;
  right   = tree_root->right;

  /* --- Zuerst die beiden Zweige in Strings verwandeln: --- */
  if( left != NULL ){
     /* --- Werte des linken Sohnes: --- */
     strcpy(symb_l,left->symb);
     sk_l    = left->symbkind;
     t_l     = left->token;
     p_l     = symbol[t_l].info;
     tree2string(left,str_l);
  }
  if( right != NULL ){
     /* --- Werte des rechten Sohnes: --- */
     strcpy(symb_r,right->symb);
     sk_r    = right->symbkind;
     t_r     = right->token;
     p_r     = symbol[t_r].info;
     tree2string(right,str_r);
  }

  /* --- Beide Strings ueber die Wurzel verknuepfen, ggf. Klammern: --- */

  if( t_m <= hos+1 ){
     /* --- Alle Operatoren mit zwei Eingaengen (von '<' bis '^'): --- */

     /* --- Funktionen links bzgl. Prioritaet  wie '*'-Operatoren: --- */
     if( sk_l == SFUNCS || sk_l == UFUNCS) sk_l = MUL_OP;

     /* --- Dummy-Null bei neg. Vorzeichen entfernen: --- */
     if( strcmp(symb_m,"-") == 0 && strcmp(symb_l,"0") == 0 ){
        if( p_r <= p_m && t_r < hoe 
           ){ strcpy3(ret_str,"-(",str_r,")"); 
           }else{ strcpy(ret_str,"-"); strcat(ret_str,str_r); 
        } 
        return ret_str; 
     } 

     /* --- Linken String evtl. einklammern: --- */ 
     if( t_l < hoe ){ /* alle Operatoren incl. 'NOT' */ 
	if( strcmp(symb_m,"^") == 0
	    ){ /* -- bei '^' auch bei gleicher Prioritaet klammern: --*/
		 if( p_l <= p_m ){
		    strcpy3(temp_str,"(",str_l,")");
		    strcpy(str_l,temp_str);
		 }
	    }else{ /* --- sonst nur bei niedrigerer Prioritaet klammern,
		    aber nicht bei '+' mit '-' [z.B. f(x)=(1+2)-3]: --- */
		 if( (p_l < p_m) &&
		    !(strcmp(symb_l,"+")==0 && strcmp(symb_m,"-")==0) ){
		    strcpy3(temp_str,"(",str_l,")");
		    strcpy(str_l,temp_str);
		 }
	}
     }

     /* --- Rechten String evtl. einklammern: --- */
     if( t_r < hoe ){                 /*alle Operatoren incl. 'NOT' */
	if( strcmp(symb_m,"/") == 0 || strcmp(symb_m,"MOD") == 0 ||
	    strcmp(symb_m,"-") == 0 || strcmp(symb_m,"^") == 0
	    ){ /* --- klammern auch bei gleicher Prioritaet: --- */
		 if( symbol[t_r].info <= symbol[t_m].info ){
		    strcpy3(temp_str,"(",str_r,")");
		    strcpy(str_r,temp_str);
		 }
	    }else{ if( symbol[t_r].info < symbol[t_m].info ){
		    strcpy3(temp_str,"(",str_r,")");
		    strcpy(str_r,temp_str);
		 }
	}
     }

     /* --- Bei einigen Operatoren Spaces drumherum setzen: --- */
     if( strcmp(symb_m,"XOR") == 0 || strcmp(symb_m,"AND") == 0 ||
	strcmp(symb_m,"OR" ) == 0 || strcmp(symb_m,"MOD") == 0 ){
	 strcpy3(temp_str," ",symb_m," ");
	 strcpy(symb_m,temp_str);
     }

     /* --- Rueckgabestring zusammensetzen aus linkem String,
			  Verknuepfungsoperator und rechtem String: --- */
/*TIRLER C-Ausgabe */
       if (strcmp(symb_m,"^") == 0)
	 { return strcpy5(ret_str, "pow(" , str_l, "," , str_r, ")" ); }
 


     return strcpy3(ret_str,str_l,symb_m,str_r);
  }

  /*-- bei "NOT" klammern, wenn Operator geringerer Prioritaet folgt: --*/
  if( strcmp(symb_m,"NOT") == 0 ){
     if( sk_r==REL_OP || sk_r==ADD_OP||sk_r==MUL_OP||strcmp(symb_r,"^")==0
	){ strcpy3(ret_str,"NOT(",str_r,")");
	     return ret_str;
	}else{ /* --- sonst nicht klammern: --- */
	     return strcpy3(ret_str,"NOT ",str_r,"");
     }
  }

  /* --- Parameterlisten von Funktionen: --- */
  if (strcmp(symb_m,",") == 0) return strcpy3(ret_str,str_l,",",str_r);

  /* --- Funktionen: Parameterlisten in Klammern setzen: --- */
  if( sk_m == SFUNCS || sk_m == UFUNCS ){
     strcpy3(ret_str,symb_m,"(",str_r);
     return strcat(ret_str,")");
  }

  /* --- Alle Endsymbole direkt zurueck: --- */
  return strcpy(ret_str,symb_m);
}
/*********************************** ************************************/

char *strcpy3(char *s1, char *s2, char *s3, char *s4)

/*  Verkettet die Strings s2, s3 und s4 und kopiert sie in s1. */

{
  *s1 = '\0';
  strcat(s1,s2); strcat(s1,s3); strcat(s1,s4);
  return s1;
}


/*********************************** ************************************/
/* TIRLER */
char *strcpy5(char *s1, char *s2, char *s3, char *s4, char *s5, char *s6)

/*  Verkettet die Strings s2, bis s6  und kopiert sie in s1. */

{
  *s1 = '\0';
  strcat(s1,s2); strcat(s1,s3); strcat(s1,s4); strcat(s1,s5); strcat(s1,s6);
  return s1;
}

/************************************************************************
 * Prozeduren zur numerischen Bewertung eines Parse-Baumes:             *
 ************************************************************************/

double tree2float(struct treenode *E_root)

/*  Erwartet in E_root den Zeiger auf den Parsebaum einer Expression und
 *  liefert als Ergebnis den numerischen Wert des Baumes.
 */

{
  double val_l, val_r;

  /* --- Werte der Wurzel: --- */
  if (E_root->left  != NULL)  val_l = tree2float(E_root->left);
  if (E_root->right != NULL)  val_r = tree2float(E_root->right);
  //  printf("ergebnis: %f - %f - %s - %f \n",val_l, val_r,
  //      symbol[E_root->token].name , symbol[E_root->token].val      );
  return (*symbol[E_root->token].vcalc)(E_root->token,val_l,val_r);
  
}
/*********************************** ************************************/

/* #pragma warn -par */    /* Compiler-Warnung "Parameter nicht benutzt" aus,*/
		      /* denn Parameter t wird hier selten gebraucht.   */
double v_dummy  (int t, double l, double r) { return 0.0; }
double v_less   (int t, double l, double r) { return (l <  r); }
double v_equal  (int t, double l, double r) { return (l == r); }
double v_greater(int t, double l, double r) { return (l  > r); }
double v_less_or(int t, double l, double r) { return (l <= r); }
double v_unequal(int t, double l, double r) { return (l != r); }
double v_grtr_or(int t, double l, double r) { return (l >= r); }
double v_or   (int t, double l, double r) { return   l || r;  }
double v_xor  (int t, double l, double r) { return !(l == r); }
double v_plus (int t, double l, double r) { return   l +  r;  }
double v_minus(int t, double l, double r) { return   l -  r;  }
double v_mul(int t, double l, double r) { return l *  r; }
double v_and(int t, double l, double r) { return l && r; }
double v_div(int t, double l, double r) { return l /  r; }
double v_mod(int t, double l, double r) { return (long)l % (long)r; }
double v_power(int t, double l, double r) { return pow(l,r); }
double v_not  (int t, double l, double r) { return !r   ; }
double v_sconst (int t, double l, double r) { return symbol[t].val; }
double v_uconst (int t, double l, double r) { return symbol[t].val; }
double v_uident (int t, double l, double r) { return symbol[t].val; }
/* verschachtelte userdefinierte Funktionen noch nicht implementiert: */
double v_ufuncs (int t, double l, double r) { return 0; }
double v_exp (int t, double l, double r) { return exp (  r); }
double v_ln  (int t, double l, double r) { return log (  r); }
double v_log (int t, double l, double r) { return log (r)/log(l); }
double v_sin (int t, double l, double r) { return sin (  r); }
double v_cos (int t, double l, double r) { return cos (  r); }
double v_tan (int t, double l, double r) { return tan (  r); }
double v_sec (int t, double l, double r) { return 1/cos( r); }
double v_sqr (int t, double l, double r) { return sqrt(  r); }
double v_abs (int t, double l, double r) { return abs (  r); }
/* #pragma warn .par      */                /* Compiler-Warnungen wieder ein */

/************************************************************************
 * Prozeduren zur analytischen Ableitung eines Parse-Baumes:            *
 ************************************************************************/

void gradient(struct treenode *root)

/*  Erwartet als Input in root den Zeiger auf eine komplette Funktions-
 *  definition ("fname(x,y,...)=...").
 *  Alle partiellen Ableitungen werden auf dem Bildschirm ausgegeben.
 */

{
  struct treenode *DFD_root = root->left; /* zgt. auf DefFuncDesignator */
  struct treenode *DP_root  =DFD_root->right; /*  auf DefParameterlist  */

  printf("\n\nGradient:");
  nxt_part_derivation(DP_root,root);
}
/************************************************************************/

void nxt_part_derivation(struct treenode *DP_root,struct treenode *root)

/*  Hangelt sich rekursiv durch die Parameterliste DP_root der Funktions-
 *  definition, auf die root zeigt.
 *  Alle partiellen Ableitungen werden auf dem Bildschirm ausgegeben.
 */

{
  struct treenode *d_tree;
  char par[SYMBLENGTH], deriv[MAXLENGTH];

  strcpy(par,DP_root->symb);
  if( strcmp(par,",") == 0 ){
     nxt_part_derivation(DP_root->left,root);
     strcpy(par,DP_root->right->symb);
  }
  d_tree = part_deriv(root,par,deriv);
  printf("\n%s",deriv);
}
/*********************************** ************************************/

struct treenode *part_deriv(struct treenode *root,char *par,char *ret_str)

/*  Erwartet als Input in root den Zeiger auf eine komplette Funktions-
 *  definition ("fname(x,y,...)=...") und in *par den String des
 *  Parameters, nach dem abgeleitet werden soll.
 *  Der Parsebaum der partiellen Ableitung wird zurueckgeliefert, in
 *  ret_str die partielle Ableitung als String. Der Funktionsname der
 *  Ableitung wird aus dem Namen der Stammfunktion durch Anhaengen
 *  des Parameters erzeugt ("fname_x(...)=...").
 */

{
  struct treenode *E_root  =root->right;    /*zeigt auf Expression      */
  struct treenode *DFD_root=root->left;    /*zeigt auf DefFuncDesignator*/
  struct treenode *DP_root =DFD_root->right;/*zeigt auf DefParameterlist*/
  struct treenode *parsetree;
  char            temp_str[MAXLENGTH];
  int             errcode, errpos;

  strcpy(temp_str,tree2string(DP_root,ret_str)); /*Par.liste als String */
  strcpy6(ret_str,DFD_root->symb,"_",par,"(",temp_str,")=");/*neuer Name*/
  derive_expression(E_root,par,temp_str);
  strcat(ret_str,temp_str);
  parsetree = string2tree(ret_str,&errcode,&errpos);
  if( errcode
     ){ pt_error(ret_str,errcode,errpos);
     }else{ tree2string(parsetree,ret_str);
  }
  return parsetree;
}
/*********************** *********************** ************************/

char *derive_expression(struct treenode *E_root,char *par,char *ret_str)

/*  Erwartet in E_root den Zeiger auf den Parsebaum einer Expression,
 *  in par den String des Parameters, nach dem abgeleitet werden soll,
 *  und liefert in ret_str die Ableitung als String zurueck.
 *  ACHTUNG: Wegen Rekursion sehr viel Stack-Speicher-Verbrauch!!!
 *                           =====================================
 */

{
  struct treenode *right, *left;
  char            l[MAXLENGTH], dl[MAXLENGTH];/* linker Zweig/Ableitung */
  char            r[MAXLENGTH], dr[MAXLENGTH];/* rechtr Zweig/Ableitung */
  char            temp_str[MAXLENGTH];
  int             tm;

  l[0] = r[0] = dl[0] = dr[0] = temp_str[0] = ret_str[0] = '\0';
  /* --- Werte der Wurzel: --- */
  tm    = E_root->token; /* Token der Wurzel ("Mitte") */
  left  = E_root->left;
  right = E_root->right;
  if( left != NULL ){
     tree2string(left,temp_str);
     strcpy3(l,"(",temp_str,")");        /* linker Zweig, eingeklammert */
     derive_expression(left,par,temp_str);
     strcpy3(dl,"(",temp_str,")");/* Abl. d.linken Zweigs,eingeklammert */
  }
  if( right != NULL ){
     tree2string(right,temp_str);
     strcpy3(r,"(",temp_str,")");       /* rechter Zweig, eingeklammert */
     derive_expression(right,par,temp_str);
     strcpy3(dr,"(",temp_str,")");/* Abl. d.rechtn Zweigs,eingeklammert */
  }
  return (*symbol[tm].dcalc)(par,E_root,l,r,dl,dr,ret_str);
}
/**************** ****************** ****************** *****************/

char *d_dummy(char *par,struct treenode *w,char *l,char *r,char *dl,
	      char *dr,char *s)
	     { return strcpy(s,"dummy"); }

char *d_const(char *par,struct treenode *w,char *l,char *r,char *dl,
	      char *dr,char *s)  /* Ableitung von Konstanten/Parametern */
	     {
	       if( strcmp(w->symb,par) == 0
		  ){ strcpy(s,"1");        /* Abl. des Parameters = 1 */
		  }else{ strcpy(s,"0");        /* Abl. von Konstanten = 0 */
	       }
	       return s;
	     }

char *d_add(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)        /* Additionsregel: (l+r)' = l'+r'  */
	   { return strcpy3(s,dl,w->symb,dr); }

char *d_mul(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s) /* Multiplikationsregel: (lr)' = l'r+lr'  */
	   { return strcat(strcpy6(s,dl,"*",r,"+",l,"*"),dr); }

char *d_div(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s) /*Quotientenregel: (l/r)' = (l'r-lr')/r^2 */
	   {
	     strcpy6(s,"(",dl,"*",r,"-",l);
	     strcpy6(l,s,"*",dr,")/",r,"^2");
	     return strcpy(s,l);
	   }

char *d_power(char *par,struct treenode *w,char *l,char *r,char *dl,
	      char *dr,char *s) /*Potenzregel:(l^r)'=l^r*(r'LN l+r*l'/l)*/
				/*   r=const.:(l^r)'=rl^(r-1)l'         */
	     {
	       if( strstr(r,par) == 0
		  ){ /* --- im Exponenten kein unabh. Parameter: --- */
		       strcpy6(s,r,"*",l,"^(",r,"-1)*");
		       strcat (s,dl);
		  }else{ /* --- auch im Exponenten unabh. Parameter: --- */
		       strcpy6(s,l,"^",r,"*(",dr,"*LN(");
		       strcpy6(dr,s,l,")+",r,"*",dl);
		       strcpy6(s,dr,"/",l,")","","");
	       }
	       return s;
	     }

char *d_ufuncs(char *par,struct treenode *w,char *l,char *r,char *dl,
	       char *dr,char *s)    /* Abl. der benutzerdef. Funktionen */
	      {
		struct treenode *parsetree;
		parsetree = part_deriv(symbol[w->token].tree,par,s);
		return strcpy(s,parsetree->left->symb);
	      }

char *d_exp(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (EXP(r))' = r'*EXP(r)        */
	   { return strcpy3(s,dr,"*EXP",r); }

char *d_ln(char *par,struct treenode *w,char *l,char *r,char *dl,
	   char *dr,char *s)            /* (LN(r))' = r'/r              */
	  { return strcpy3(s,dr,"/",r); }

char *d_log(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (LOG(l,r))' = r'/(r*LN(l))   */
	   {
	     /*l,r,dl,dr neu berechnen,da in r die ganze Arg.liste steht*/
	     tree2string      (w->right->left ,l);
	     derive_expression(w->right->left ,par,dl);
	     tree2string      (w->right->right,r);
	     derive_expression(w->right->right,par,dr);
	     return strcpy6(s,dr,"/(",r,"*LN(",l,"))");
	   }

char *d_sin(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SIN(r))' = r'*COS(r)        */
	   { return strcpy3(s,dr,"*COS",r); }

char *d_cos(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (COS(r))' = -r'*SIN(r)       */
	   { return strcpy6(s,"-",dr,"*SIN",r,"",""); }

char *d_tan(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (TAN(r))' = r'*(SEC(r))^2    */
	   { return strcpy6(s,dr,"*(SEC",r,")^2","",""); }

char *d_sec(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SEC(r))' = r'*TAN(r)*SEC(r) */
	   { return strcpy6(s,dr,"*TAN",r,"*SEC",r,""); }

char *d_sqr(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SQR(r))' = r'/(2*SQR(r))    */
	   { return strcpy6(s,dr,"/(2*SQR(",r,"))" ,"",""); }

char *d_abs(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (ABS(r))' = "ABSerror"       */
	   { return strcpy(s,"ABSerror"); }


/**************** ****************** ****************** *****************/

char *strcpy6(char *s1, char *s2, char *s3,
	      char *s4, char *s5, char *s6, char *s7)

/*  Verkettet die Strings s2 bis s7 und kopiert sie in s1. */

{
  *s1 = '\0';
  strcat(s1,s2); strcat(s1,s3); strcat(s1,s4);
  strcat(s1,s5); strcat(s1,s6); strcat(s1,s7);
  return s1;
}

/**************** ****************** ****************** *****************/

void str_upr(char *s)

/*  Wandelt Strings in Grossbuchstaben. */

{
  do
    *s=((*s>96)&&(*s<123)) ? *s-=32 : *s;
  while(*(++s));
}




/***************************************************************************************/
struct treenode *_unur_fstr2tree(char *function, int *errcodep, int *errposp)

{ char            *fvonx;
  struct treenode *root;

 
  str_upr(function);
  fvonx= (char *)malloc(MAXLENGTH*sizeof(char));     
  strcpy(fvonx,"f(x)=");  
  strcat(fvonx,function);
  function = fvonx; 
  free(fvonx); 
  root=string2tree(function,errcodep,errposp);
  return root;

}

/***************************************************************************************/
double  _unur_fstr_eval_tree(struct treenode *E_root, double argument)
{
  double           result;
  struct treenode *froot;
  int             ftok,xtok;

  xtok=find_index("X",uis,ufe,"");
  ftok=find_index("F",uis,ufe,"");
  froot=(*symbol[ftok].tree).right;             /* Achtung Fehler in Beschreibung !!! */
  symbol[xtok].val= argument;
  result=tree2float(froot);
  return result;
  }



/***************************************************************************************/
char *Ntree2string(struct treenode *tree_root, char *ret_str)

{
  double           result;
  struct treenode  *froot;
  int              ftok;

  ftok=find_index("F",uis,ufe,"");
  froot=(*symbol[ftok].tree).right;          

  tree2Cstring(froot,ret_str); 
  return ret_str;
}

/***************************************************************************************/
struct treenode *_unur_fstr_make_derivative(struct treenode *root)

{ 
  struct treenode *E_root  =root->right;    /*zeigt auf Expression      */
  // struct treenode *DFD_root=root->left;    /*zeigt auf DefFuncDesignator*/
  // struct treenode *DP_root =DFD_root->right;/*zeigt auf DefParameterlist*/
   struct treenode *parsetreeh,*parsetree;
   char            temp_str[MAXLENGTH];
   int             errcode, errpos;
  
   char x='X';
   char            ret_str[MAXLENGTH];


   //   printf("\nAbleitungxxxx:%s\n",input_string);
   parsetreeh= part_deriv(root,&x,ret_str);
    // printf("\nAbleitungxxxx:%s\n",ret_str);
   return parsetreeh;


  


   //   strcpy(ret_str,"f(x)=");
   // printf("\nAbleitung1:%s\n",temp_str);
   // part_deriv(E_root,&x,temp_str);
   //printf("\nAbleitung2:%s\n",temp_str);
   strcat(ret_str,temp_str);
   parsetree = string2tree(ret_str,&errcode,&errpos);
   // if( (errcode>0) return NULL;
  //    ){ pt_error(ret_str,errcode,errpos);
  //   }else{ tree2string(parsetree,ret_str);
  // }

   // derive_expression(E_root,*x,temp_str);
   // parsetree = string2tree(temp_str,&errcode,&errpos);
    
   // printf("\nAbleitung:%s\n",temp_str);
   // printf("\nAbleitung:%s\n",ret_str);
    printf("\nx:%c\n",x);
    printf("-------------------------------------------");
   
   return parsetree;
}

/***************************************************************************************/












