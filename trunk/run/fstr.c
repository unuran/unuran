
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

/* --- Prototypen: --- */
void   _unur_fstr_init   (void);
void  show_symb_tab      (void);
void clear_symbol_area(int start, int end);
static char  get_ch_after_spaces(char *function, int *scanposp);
static char *get_ds             (char *function, int  *scanposp, char *ds);
static char *get_sf             (char *function, int  *scanposp, char *sf);
static int   get_uc_symbol      (char *function, int *scanposp, char *uc);
int   nxt_symbol         (char function[], int *scanposp, char symb[],
			      int *tokenp, int *symbkindp, int *errposp);
static int   get_id_symbol      (char *function, int  *scanposp, char *id);
static int   get_le_symbol      (char *function, int  *scanposp, char *le);
static int   get_ge_symbol      (char *function, int  *scanposp, char *ge);
int   find_kind          (int token);
int   find_index         (char *symb, int start, int end, int nxt_c);

/* --- Bereichs-Start- und -Endemarkierungen in der Symboltabelle: --- */
int ros, roe, aos, aoe, mos, moe, hos, hoe, oss, ose, scs, sce;
int ucs, uce, uis, uie, ufs, ufe, sfs, sfe;





/* --- forward-declarations aus Modul CALC.C; Bewertungsfunktionen: --- */
extern double v_exp(), v_dummy  (), v_or   (), v_mul(), v_sconst();
extern double v_ln (), v_less   (), v_xor  (), v_and(), v_uconst();
extern double v_log(), v_equal  (), v_plus (), v_div(), v_uident();
extern double v_sin(), v_greater(), v_minus(), v_mod();
extern double v_cos(), v_less_or();
extern double v_tan(), v_unequal(), v_power(), v_ufuncs();
extern double v_sec(), v_grtr_or(), v_not  ();
extern double v_sqr();
extern double v_abs();

/* --- forward-declarations aus Modul CALC.C; Ableitungsfunktionen: --- */
extern char *d_exp(), *d_dummy ();
extern char *d_ln (), *d_add   ();
extern char *d_log(), *d_mul   ();
extern char *d_sin(), *d_div   ();
extern char *d_cos(), *d_power ();
extern char *d_tan(), *d_const ();
extern char *d_sec(), *d_par   ();
extern char *d_sqr(), *d_ufuncs();
extern char *d_abs();


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


/************************************************************************/
/* --- SYMBOLTABELLE, mit Funktionsadressen --- */

struct symbols symbol[] = {   
  {"ros",0, .0,v_dummy,  d_dummy, 0},            /* RelationalOperators */
  {"<",  1, .0,v_less,   d_const, 0}, {"=",  1, .0,v_equal,  d_const, 0},
  {">",  1, .0,v_greater,d_const, 0}, {"<=", 1, .0,v_less_or,d_const, 0},
  {"<>", 1, .0,v_unequal,d_const, 0}, {">=", 1, .0,v_grtr_or,d_const, 0},
  {"roe",0, .0,v_dummy,  d_dummy, 0},

  {"aos",0, .0,v_dummy,  d_dummy, 0},                /* AddingOperators */
  {"OR", 2, .0,v_or,     d_const, 0}, {"XOR",2, .0,v_xor,    d_const, 0},
  {"+",  2, .0,v_plus,   d_add,   0}, {"-",  3, .0,v_minus,  d_add,   0},
  {"aoe",0, .0,v_dummy,  d_dummy, 0},

  {"mos",0, .0,v_dummy,  d_dummy, 0},           /* MultiplyingOperators */
  {"*",  4, .0,v_mul,    d_mul,   0}, {"AND",4, .0,v_and,    d_const, 0},
  {"/",  4, .0,v_div,    d_div,   0}, {"MOD",4, .0,v_mod,    d_const, 0},
  {"moe",4, .0,v_dummy,  d_dummy, 0},

  {"hos",0, .0,v_dummy,  d_dummy, 0},          /* HighPriorityOperators */
  {"^",  5, .0,v_power,  d_power, 0}, {"NOT",6, .0,v_not,    d_const, 0},
  {"hoe",0, .0,v_dummy,  d_dummy, 0},

  {"oss",0, .0,v_dummy,  d_dummy, 0},                   /* OtherSymbols */
  {"(",  0, .0,v_dummy,  d_dummy, 0}, {")",  0, .0,v_dummy,  d_dummy, 0},
  {",",  0, .0,v_dummy,  d_dummy, 0},
  {"ose",0, .0,v_dummy,  d_dummy, 0},

  {"scs",0, .0,v_dummy,  d_dummy, 0},                /* SystemConstants */
  {"0",  0,0.0,v_sconst, d_const, 0}, {"1",  0,1.0,v_sconst, d_const, 0},
  {"2",  0,2.0,v_sconst, d_const, 0}, {"3",  0,3.0,v_sconst, d_const, 0},
  {"4",  0,4.0,v_sconst, d_const, 0}, {"5",  0,5.0,v_sconst, d_const, 0},
  {"6",  0,6.0,v_sconst, d_const, 0}, {"7",  0,7.0,v_sconst, d_const, 0},
  {"8",  0,8.0,v_sconst, d_const, 0}, {"9",  0,9.0,v_sconst, d_const, 0},
  {"PI", 0,3.1415927,v_sconst, d_const, 0},
  {"E",  0,2.7182818,v_sconst, d_const, 0},
  {"sce",0, .0,v_dummy,  d_dummy, 0},

  {"ucs",0, .0,v_dummy,  d_dummy, 0},           /* UserdefinedConstants */
  {"",   0, .0,v_uconst, d_const, 0}, {"",   0, .0,v_uconst, d_const, 0},
  {"",   0, .0,v_uconst, d_const, 0}, {"",   0, .0,v_uconst, d_const, 0},
  {"",   0, .0,v_uconst, d_const, 0}, {"",   0, .0,v_uconst, d_const, 0},
  {"",   0, .0,v_uconst, d_const, 0}, {"",   0, .0,v_uconst, d_const, 0},
  {"uce",0, .0,v_dummy,  d_dummy, 0},

  {"uis",0, .0,v_dummy,  d_dummy, 0},         /* UserdefinedIdentifiers */
  {"",   0, .0,v_uident, d_const, 0}, {"",   0, .0,v_uident, d_const, 0},
  {"",   0, .0,v_uident, d_const, 0}, {"",   0, .0,v_uident, d_const, 0},
  {"",   0, .0,v_uident, d_const, 0}, {"",   0, .0,v_uident, d_const, 0},
  {"",   0, .0,v_uident, d_const, 0}, {"",   0, .0,v_uident, d_const, 0},
  {"uie",0, .0,v_dummy,  d_dummy, 0},

  {"ufs",0, .0,v_dummy,  d_dummy, 0},           /* UserdefinedFunctions */
  {"",   0, .0,v_ufuncs, d_ufuncs,0}, {"",   0, .0,v_ufuncs, d_ufuncs,0},
  {"",   0, .0,v_ufuncs, d_ufuncs,0}, {"",   0, .0,v_ufuncs, d_ufuncs,0},
  {"",   0, .0,v_ufuncs, d_ufuncs,0}, {"",   0, .0,v_ufuncs, d_ufuncs,0},
  {"",   0, .0,v_ufuncs, d_ufuncs,0}, {"",   0, .0,v_ufuncs, d_ufuncs,0},
  {"ufe",0, .0,v_dummy,  d_dummy, 0},

  {"sfs",0, .0,v_dummy,  d_dummy, 0},                /* SystemFunctions */
  {"EXP",1, .0,v_exp,    d_exp,   0}, {"LN", 1, .0,v_ln,     d_ln,    0},
  {"LOG",2, .0,v_log,    d_log,   0}, {"SIN",1, .0,v_sin,    d_sin,   0},
  {"COS",1, .0,v_cos,    d_cos,   0}, {"TAN",1, .0,v_tan,    d_tan,   0},
  {"SEC",1, .0,v_sec,    d_sec,   0}, {"SQR",1, .0,v_sqr,    d_sqr,   0},
  {"ABS",1, .0,v_abs,    d_abs,   0},
  {"sfe",0, .0,v_dummy,  d_dummy, 0},

  {"end",0, .0,v_dummy,  d_dummy, 0},
};

/************************************************************************/

void  _unur_fstr_init(void)

/*  Ermittelt die Indices der Start- und Endsymbole in der Symboltabelle
 *  und loescht alle benutzerdefinierten Konstanten und Identifier.
 */

{
  int i = -1;
  char *s;

  do {
     s = symbol[++i].name;
     if         (strcmp(s,"ros") == 0) ros = i;
	else if (strcmp(s,"roe") == 0) roe = i;
	else if (strcmp(s,"aos") == 0) aos = i;
	else if (strcmp(s,"aoe") == 0) aoe = i;
	else if (strcmp(s,"mos") == 0) mos = i;
	else if (strcmp(s,"moe") == 0) moe = i;
	else if (strcmp(s,"hos") == 0) hos = i;
	else if (strcmp(s,"hoe") == 0) hoe = i;
	else if (strcmp(s,"oss") == 0) oss = i;
	else if (strcmp(s,"ose") == 0) ose = i;
	else if (strcmp(s,"scs") == 0) scs = i;
	else if (strcmp(s,"sce") == 0) sce = i;
	else if (strcmp(s,"ucs") == 0) ucs = i;
	else if (strcmp(s,"uce") == 0) uce = i;
	else if (strcmp(s,"uis") == 0) uis = i;
	else if (strcmp(s,"uie") == 0) uie = i;
	else if (strcmp(s,"ufs") == 0) ufs = i;
	else if (strcmp(s,"ufe") == 0) ufe = i;
	else if (strcmp(s,"sfs") == 0) sfs = i;
	else if (strcmp(s,"sfe") == 0) sfe = i;
  } while (!(strcmp(s,"end") == 0));
  /* --- Benutzerdefinierte Eintraege loeschen: --- */
  clear_symbol_area(ucs,uce);
  clear_symbol_area(uis,uie);
  clear_symbol_area(ufs,ufe);
}

/*********************************** ************************************/

void clear_symbol_area(int start, int end)

/* Loescht Symbole aus der Symboltabelle von Index start bis Index end  */
{
  symbol[start].info=0; /* zeigt immer auf den letzten belegten Eintrag */
  do {
      symbol[++start].name[0] = '\0';
      symbol[  start].val     = 0.0;
  } while (!(start >= end-1));
}

/************************************************************************/

void show_symb_tab(void)         /* Gibt die aktuelle Symboltabelle aus */
{
  int i = 0;

  printf("\nSymboltabelle:\n");
  do {
       printf("\n%4d: %7s inf:%3d val:%f ",
       i, symbol[i].name,symbol[i].info,symbol[i].val);
  } while (!(strcmp(symbol[++i].name,"end") == 0));
  printf("\nros=%4d hos=%4d ucs=%4d sfs=%4d",ros,hos,ucs,sfs);
  printf("\nroe=%4d hoe=%4d uce=%4d sfe=%4d",roe,hoe,uce,sfe);
  printf("\naos=%4d oss=%4d uis=%4d"        ,aos,oss,uis    );
  printf("\naoe=%4d ose=%4d uie=%4d"        ,aoe,ose,uie    );
  printf("\nmos=%4d scs=%4d ufs=%4d"        ,mos,scs,ufs    );
  printf("\nmoe=%4d sce=%4d ufe=%4d"        ,moe,scs,ufe    );
  printf("\n");
}
/************************************************************************/

int nxt_symbol(char function[],int *scanposp, char symb[], int *tokenp,
					     int *symbkindp, int *errposp)

/*  Liefert aus dem String function[] in symb[] das naechste Symbol ab
 *  der Position *scanposp. Nach der Bearbeitung zeigt *scanposp auf das
 *  Zeichen, das dem ermittelten Symbol unmittelbar folgt. Leerzeichen
 *  werden ignoriert. In *tokenp wird der Index des Symbols in der
 *  Symboltabelle. geliefert. *symbkindp ist die Art des Symbols
 *  (UIDENT, SCONST etc.), *errposp zeigt auf das Zeichen, auf das
 *  *scanposp beim Aufruf zeigte.
 *  Der zurueckgelieferte Funktionswert ist der Fehlercode.
 */

 {
  char c;
  int  dummy = 0, nxt_c, errcode = 0;

  if (*scanposp > strlen(function)) return 1;
  c = get_ch_after_spaces(function, scanposp);
  *errposp = *scanposp;
  if ( (c >= '0' && c <= '9') || c == '.') {           /* UnsignedConstant */
     errcode = get_uc_symbol(function,scanposp,symb);
     if ((*tokenp = find_index(symb,scs,uce,dummy)) <= 0) return 2; }
  else if (c >=  'A' && c <= 'Z') {                       /* Identifier */
     errcode = get_id_symbol(function,scanposp,symb);
     nxt_c = function[*scanposp];
     if ((*tokenp = find_index(symb,aos,sfe,nxt_c)) <= 0) return 3; }
  else if (c == '<') {             /* kleiner, kleiner gleich, ungleich */
     errcode = get_le_symbol(function,scanposp,symb);
     *tokenp = find_index(symb,ros,roe,dummy);                      }
  else if (c == '>') {                     /* groesser, groesser gleich */
     errcode = get_ge_symbol(function,scanposp,symb);
     *tokenp = find_index(symb,ros,roe,dummy);                      }
  else { symb[0] = c; symb[1] = '\0';           /* alle anderen Zeichen */
     (*scanposp)++;
     *tokenp = find_index(symb,ros,sfe,dummy);
  }
  *symbkindp = find_kind(*tokenp);         /* Art des Symbols ermitteln */
  if (*symbkindp == 0) errcode = 4;        /*       unbekanntes Zeichen */
  return errcode;
 }
/*********************************** ************************************/

static char get_ch_after_spaces(char *function, int *scanposp)

/*  šberliest Spaces und liefert n„chstes Zeichen im String function[].
 *  scanpos zeigt anschliessend auf dieses erste Nicht-Space-Zeichen.
 */

{
  char c;

  while ((c = function[*scanposp]) == ' ') (*scanposp)++;
  return c;
}
/*********************************** ************************************/

static int get_uc_symbol(char *function, int *scanposp, char *uc)

/* UnsignedConstant ::= UnsignedInteger | UnsignedReal
 * UnsignedInteger  ::= DigitSequence
 * UnsignedReal ::= UnsignedInteger "." DigitSequence ["E" ScaleFactor] |
 *                  UnsignedInteger "E" ScaleFactor
 */

{
  char s[SYMBLENGTH], c;

  strcpy(uc,get_ds(function,scanposp,s));
  c = function[*scanposp];
  if( c == 'E' ){ /* wenn "E" folgt, wird noch ein ScaleFactor gelesen */
     strcat(uc,"E");
     (*scanposp)++;
     strcat(uc,get_sf(function,scanposp,s));
  }
  if( c == '.' ){       /* wenn "." folgt, koennen noch Ziffern kommen */
     strcat(uc,".");
     (*scanposp)++;
     strcat(uc,get_ds(function,scanposp,s));
     c = function[*scanposp];
     if( c == 'E' ){  /* wenn "E" folgt, wird noch ScaleFactor gelesen */
	strcat(uc,"E");
	(*scanposp)++;
	strcat(uc,get_sf(function,scanposp,s));
     }
  }
  return 0;
}
/*********************** *********************** ************************/

static char *get_ds(char *function, int *scanposp, char *ds)

/*  DigitSequence   ::= Digit { Digit }
 *  Digit           ::= "0" | "1" | "2" | ... | "8" | "9"
 */

{
  int i = 0;

  while ( (ds[i] = function[*scanposp]) >= '0' && ds[i] <= '9' ) {
     i++;
     (*scanposp)++;
  }
  ds[i] = '\0';
  return ds;
}
/*********************** *********************** ************************/

static char *get_sf(char *function, int *scanposp, char *sf)

/*  ScaleFactor ::= [Sign] DigitSequence
 *  Sign        ::= "+" | "-"
 */

{
  char ds[SYMBLENGTH], c;

  sf[0] = sf[1] = '\0';
  if( (c=function[*scanposp]) == '+' || c == '-' ){
     sf[0] = c;
     (*scanposp)++;
  }
  return strcat(sf,get_ds(function,scanposp,ds));
}
/*********************************** ************************************/

static int get_id_symbol(char *function, int *scanposp, char *id)

/*  Identifier ::= Letter { Letter | Digit }
 *  Letter     ::= "A" | "B" | ... | "Z" | "_"
 *  Digit      ::= "0" | "1" | ... | "9"
 */

{
  while ( ((*id = function[(*scanposp)]) >= 'A' && *id <= 'Z')
     || *id == '_' ||              ( *id >= '0' && *id <= '9')) {
     id++;
     (*scanposp)++;
  }
  *id = '\0';
  return 0;
}
/*********************************** ************************************/

static int get_le_symbol(char *function, int *scanposp, char *le)
						 /*  *le="<","<=","<>". */
{
  while ((*le = function[*scanposp]) == '<' || *le == '=' || *le == '>') {
    le++;
    (*scanposp)++;
  }
  *le = '\0';
  return 0;
}
/*********************************** ************************************/

static int get_ge_symbol(char *function, int *scanposp, char *ge)
						 /*  *ge= ">" oder ">=".*/
{
  while ( (*ge = function[*scanposp]) == '>' || *ge == '=') {
    ge++;
    (*scanposp)++;
  }
  *ge = '\0';
  return 0;
}
/*********************************** ************************************/

int find_kind(int token)            /*  Ermittelt Art des Tokens token. */
{
  if      (token > ros && token < roe) return REL_OP;
  else if (token > aos && token < aoe) return ADD_OP;
  else if (token > mos && token < moe) return MUL_OP;
  else if (token > hos && token < hoe) return HPR_OP;
  else if (token > oss && token < ose) return OTHERS;
  else if (token > scs && token < sce) return SCONST;
  else if (token > ucs && token < uce) return UCONST;
  else if (token > uis && token < uie) return UIDENT;
  else if (token > ufs && token < ufe) return UFUNCS;
  else if (token > sfs && token < sfe) return SFUNCS;
  else                                 return NOSYMB;
}
/*********************************** ************************************/

int find_index(char *symb, int start, int end, int nxt_c)

/*  Ermittelt den Index von symb[] in der Symboltabelle *symbol[] und
 *  traegt unbekannte Symbole ggf. in die Tabelle ein.
 *  In nxt_c wird der Character nach symb uebergeben; wenn hier '(' steht,
 *  ist der Identifier eine Funktion.
 */

{
  int i, nxt_free_place;

  for (i = start + 1; i < end; i++)
 {     if (strcmp(symb,symbol[i].name) == 0) break;
 // printf("\n %d: %s - %s -> %s",i,symb, symbol[i].name, symbol[64].name);
 }  
if (i < end ) return i;

  /*  symb ist noch nicht in Symboltabelle => eintragen, wenn Platz: */
  if( start == scs ){
     /* --- Symbol ist Userdefined Constant: --- */
     nxt_free_place = ucs + symbol[ucs].info + 1;
     if( nxt_free_place >= uce
	){ /* --- kein Platz mehr: --- */
	     return 0;
	}else{ strcpy(symbol[nxt_free_place].name,symb);
	     symbol[nxt_free_place].val = atof(symb);
	     symbol[ucs].info++;     /* hier Zeiger auf letzten Eintrag */
	     return nxt_free_place;
     }
  }
  if( start == aos ){
     /* --- Symbol ist Identifier: --- */
     if( nxt_c == '('
	){ /* --- Symbol ist Userdefined Function --- */
	     nxt_free_place = ufs + symbol[ufs].info + 1;
	     if( nxt_free_place >= ufe
		){ /* --- kein Platz mehr: --- */
		     return 0;
		}else{ strcpy(symbol[nxt_free_place].name,symb);
		     symbol[nxt_free_place].info = 0;
		     symbol[ufs].info++;
		     return nxt_free_place;
	     }
	}else{ /* --- Symbol ist Userdefined Identifier --- */
	     nxt_free_place = uis + symbol[uis].info + 1;
	     if( nxt_free_place >= uie
		){ /*** kein Platz mehr: ***/
		     return 0;
		}else{ strcpy(symbol[nxt_free_place].name,symb);
		     symbol[nxt_free_place].val = 0.0;
		     symbol[uis].info++;
		     return nxt_free_place;
	     }
     }
  }
  return 0;
}



/* --- Import aus SCANNER.C: --- */ 
extern int ros, roe, aos, aoe, mos, moe, hos, hoe, oss, ose, scs, sce; 
extern int ucs, uce, uis, uie, ufs, ufe, sfs, sfe; 

extern struct symbols symbol[]; 

extern void init_scanner(void);
extern void show_symb_tab(void);
extern int  nxt_symbol(char function[], int *scanposp,  char symb[],
                       int  *tokenp,    int *symbkindp, int *errposp);
extern int  find_kind(int token);
extern int  find_index(char *symb, int start, int end, int nxt_c);

/* --- Prototypen: --- */ 
struct treenode *string2tree(char *function, int *errcodep, int *errposp);
static struct treenode *FuncDefinition   (char *f, int *spp, int *ecp, int *epp);
static struct treenode *DefFuncDesignator(char *f, int *spp, int *ecp, int *epp, 
                                                                int *ftp);
static struct treenode *DefParameterlist (char *f, int *spp, int *ecp, int *epp, 
                                                            int *paranzp);
static struct treenode *Expression       (char *f, int *spp, int *ecp, int *epp);
static struct treenode *SimpleExpression (char *f, int *spp, int *ecp, int *epp);
static struct treenode *VTerm            (char *f, int *spp, int *ecp, int *epp);
static struct treenode *Term             (char *f, int *spp, int *ecp, int *epp);
static struct treenode *Factor           (char *f, int *spp, int *ecp, int *epp);
static struct treenode *bas_exp          (char *f, int *spp, int *ecp, int *epp);
static struct treenode *FuncDesignator   (char *f, int *spp, int *ecp, int *epp);
static struct treenode *ActualParameterlist(char *f, int *spp, int *ecp, 
                                                 int *epp, int corr_panz);
static struct treenode *create_node(char *symb, int token, int symbkind,
                double val, struct treenode *left, struct treenode *right);
static struct treenode *tnalloc (void);
static struct treenode *set_err (int err_nr, int *errcodep);
char            *readln  (char *s);
static char            *erase_back_blancs (char *s); 
static int             simplification(char *symb, int t, struct treenode *l, 
                                                      struct treenode *r);
void            check_reorg (struct treenode *root);
void            show_tree   (struct treenode *root);
void            _unur_fstr_free (struct treenode *root);                            
               
void            pt_error    (char *function, int errcode, int errpos);

char *errorstrings[] = { 
     "OK",                                                /* 0*/ 
     "scan pointer too big",                              /* 1*/ 
     "unsigned constant area full",                       /* 2*/ 
     "identifier area full",                              /* 3*/ 
     "4", "5", "6", "7", 
     "8", "9", "10", 
     "operator expected        in string2tree",           /*11*/ 
     "expected symbol: '='     in FunctionDefinition",    /*12*/ 
     "user function expected   in DefFuncDesignator",     /*13*/ 
     "expected symbol: '('     in DefFuncDesignator",     /*14*/ 
     "expected symbol: ')'     in DefFuncDesignator",     /*15*/ 
     "user identifier expected in DefParameterlist",      /*16*/ 
     "user identifier expected in DefParameterlist",      /*17*/ 
     "expected symbol: ')'     in bas_exp",               /*18*/ 
     "unknown symbol           in bas_exp",               /*19*/ 
     "function expected        in FuncDesignator",        /*20*/ 
     "expected symbol: '('     in FuncDesignator",        /*21*/ 
     "expected symbol: ')'     in FuncDesignator",        /*22*/ 
     "expected symbol: ')'     in ActualParameterlist",   /*23*/ 
     "expected symbol: ','     in ActualParameterlist",   /*24*/ 
}; 
/************************************************************************/ 

struct treenode *string2tree(char *function, int *errcodep, int *errposp) 

/*  Wandelt den String function[] in einen Parse-Baum um. function[] muss 
 *  in der Form f(a,b,...)=Expression vorliegen. In *errcodep und *errposp 
 *  wird im Fehlerfall eine Fehlernummer und die -position im String 
 *  function[] zurueckgegeben, sonst ist *errcodep Null. 
 */ 

{ 
  struct treenode *root;
  char            f[MAXLENGTH]; 
  int             scanpos = 0; 

  *errcodep = 0; 
  strcpy(f,function);     /* damit function[] unveraendert bleiben kann */ 
  erase_back_blancs(f);   /* fuer Fehler Nr. 11 noetig! */ 
  str_upr(f);
  root = FuncDefinition(f, &scanpos, errcodep, errposp); 
  if (scanpos != strlen(f)) return set_err(11,errcodep); 
  return root; 
} 
/*********************************** ************************************/ 

static char *erase_back_blancs(char *s) /* Entf. in s[] hintere Blancs. */ 

{ 
  unsigned long i = strlen(s) - 1; 

  while (s[i] == ' ') s[i--] = '\0'; 
  return s; 
 } 
/*********************************** ************************************/ 

static struct treenode *FuncDefinition(char *f,int *spp,int *ecp,int *epp) 

/*  FuncDefinition ::= DefFuncDesignator "=" Expression 
 * 
 *                    "=" 
 *                   /   \ 
 *  DefFuncDesignator     Expression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, ft; 
  double           v = 0.0; 

  l = DefFuncDesignator(f,spp,ecp,epp,&ft); if (*ecp) return NULL; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,"=") != 0) return set_err(12,ecp); 
  r = Expression(f,spp,ecp,epp); if (*ecp) return NULL; 
  root = create_node(symb,t,sk,v,l,r); 
  symbol[ft].tree = root; 
  return root; 
} 
/*********************** *********************** ************************/ 

static struct treenode *DefFuncDesignator(char *f, int *spp, int *ecp, 
                                                   int *epp, int *ftp) 

/*  DefFuncDesignator ::= Identifier "(" DefParameterlist ")" 
 * 
 *       Identifier 
 *      /          \ 
 *  NULL            DefParameterlist 
 */ 

{ 
  struct treenode *root, *r; 
  char            symb[SYMBLENGTH], fsymb[SYMBLENGTH]; 
  int             t, sk, fsk, paranz; 
  double           fv = 0.0; 

  *ecp = nxt_symbol(f,spp,fsymb,ftp,&fsk,epp); 
  if (fsk != UFUNCS) return set_err(13,ecp); 
  if( symbol[*ftp].info != 0 ){ 
     /* --- Funktion war schon vorhanden => alten Baum loeschen: --- */ 
     _unur_fstr_free(symbol[*ftp].tree); 
  } 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,"(") != 0) return set_err(14,ecp); 
  r = DefParameterlist(f,spp,ecp,epp,&paranz); if (*ecp) return NULL; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,")") != 0) return set_err(15,ecp); 
  symbol[*ftp].info = paranz; 
  root = create_node(fsymb,*ftp,fsk,fv,NULL,r); 
  return root; 
} 
/**************** ****************** ****************** *****************/ 

static struct treenode *DefParameterlist(char *f, int *spp, int *ecp, 
                                                  int *epp, int *paranzp) 

/*  DefParameterlist ::= "(" Identifier { "," Identifier } ")" 
 * 
 *       Identifier                                    "," 
 *      /          \      oder:                       /   \ 
 *  NULL            NULL         more identifiers tree     Identifier 
 *                                                        /          \ 
 *                                                    NULL            NULL 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             t, ct, sk, scanpos; 
  double           v = 0.0; 

  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (sk != UIDENT) return set_err(16,ecp); 
  *paranzp = 1; 
  root = create_node(symb,t,sk,v,NULL,NULL); 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&ct,&sk,epp); 
  while (strcmp(symb,",") == 0) {
     l = root; 
     *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
     if (sk != UIDENT) return set_err(17,ecp); 
     r = create_node(symb,t,sk,v,NULL,NULL); 
     (*paranzp)++; 
     root = create_node(",",ct,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(f,spp,symb,&ct,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 
/*********************** *********************** ************************/ 

static struct treenode *set_err (int err_nr, int *errcodep) 

/*  Setzt *errcodep gleich err_nr und gibt NULL zurueck. */ 

{ 
  printf("\n%d %s",err_nr,errorstrings[err_nr]); 
  if (*errcodep == 0) *errcodep = err_nr;   /* nur den 1. Fehler melden */ 
  return NULL; 
} 
/*********************** *********************** ************************/ 

static struct treenode *Expression(char *f, int *spp, int *ecp, int *epp) 

/* Expression ::= SimpleExpression [ RelationalOperator SimpleExpression ] 
 * 
 *                                     RelationalOperator 
 * SimpleExpression  oder:            /                  \ 
 *                    SimpleExpression                    SimpleExpression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos; 
  double           v = 0.0; 

  l = SimpleExpression(f,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if( sk == REL_OP 
     ){ r = SimpleExpression(f,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ *spp = scanpos; 
          root = l; 
  } 
  return root; 
} 
/**************** ****************** ****************** *****************/ 

static struct treenode *SimpleExpression(char *f, int *spp, 
                                                  int *ecp, int *epp) 

/*  SimpleExpression ::= VTerm { AddingOperator Term } 
 * 
 *                                AddingOperator 
 *  VTerm  oder:                 /              \ 
 *                more terms tree                Term 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos; 
  double           v = 0.0; 

  root = VTerm(f,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  while (sk == ADD_OP) {
     l = root; 
     r = Term(f,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(symb,t,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 
/************ ************** *************** ************** *************/ 

static struct treenode *VTerm(char *f, int *spp, int *ecp, int *epp) 

/*  Vterm ::= [ "+" | "-" ] Term 
 * 
 *                        "-" 
 *  Term  oder:          /   \ 
 *                    "0"     Term 
 *                   /   \ 
 *               NULL     NULL 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos = *spp; 
  double           v = 0.0; 

  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if( strcmp(symb,"-") == 0 
     ){ /* --- Term hat neg. Vorzeichen =>     --- */ 
          /* --- Vorzeichen in Operator wandeln: --- */ 
          l = create_node("0",scs+1,SCONST,0.0,NULL,NULL); 
          r = Term(f,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ /* --- Term hat pos. oder kein Vorzeichen: --- */ 
          if (strcmp(symb,"+") != 0) *spp = scanpos;  /* "+" ignorieren */ 
          root = Term(f,spp,ecp,epp); if (*ecp) return NULL; 
  } 
  return root; 
} 
/********** *********** ************ ************ *********** ***********/ 

static struct treenode *Term(char *f, int *spp, int *ecp, int *epp) 

/*  Term ::= Factor { MultiplyingOperator Factor } 
 * 
 *                                   MultiplyingOperator 
 *  Factor  oder:                   /                   \ 
 *                 more factors tree                     Factor 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos; 
  double           v = 0.0; 

  root = Factor(f,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  while (sk == MUL_OP) {
     l = root; 
     r = Factor(f,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(symb,t,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 
/******** ********* ********** *********** ********** ********* *********/ 

static struct treenode *Factor(char *f, int *spp, int *ecp, int *epp) 

/*  Factor ::= Base [ "^" Exponent ] 
 * 
 *                          "^" 
 *  bas_exp  oder:         /   \ 
 *                  bas_exp     bas_exp 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos; 
  double           v = 0.0; 

  l = bas_exp(f,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if( strcmp(symb,"^") == 0 
     ){ r = bas_exp(f,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ *spp = scanpos; 
          root = l; 
  } 
  return root; 
} 
/******* ******** ******** ********* ********* ******** ******** ********/ 

static struct treenode *bas_exp(char *f, int *spp, int *ecp, int *epp) 

/*  Base ::= Exponent ::= UnsignedConstant | Identifier | FuncDesignator | 
 *                        "NOT" Base | "(" Expression ")" 
 * 
 *       UnsignedConstant                Identifier 
 *      /                \     oder     /          \     oder 
 *  NULL                  NULL      NULL            NULL 
 * 
 *                                         "NOT" 
 *        FuncDesignator       oder       /     \        oder Expression 
 *                                    NULL       bas_exp 
 */ 

{ 
  struct treenode *root, *r; 
  struct treenode *bas_exp(); 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos = *spp; 
  double           v = 0.0; 

  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if( sk==SCONST || sk==UCONST || sk==UIDENT || strcmp(symb,"NOT")==0 
     ){ /* --- neuen Knoten bilden: --- */ 
          if( strcmp(symb,"NOT") == 0 
             ){ r = bas_exp(f,spp,ecp,epp); if (*ecp) return NULL; 
             }else{ r = NULL; 
          } 
          root = create_node(symb,t,sk,v,NULL,r); 
     }else{ /* --- neuer Knoten nicht noetig: --- */ 
          if( sk == UFUNCS || sk == SFUNCS 
             ){ *spp = scanpos; 
                  root = FuncDesignator(f,spp,ecp,epp); 
                  if (*ecp) return NULL; 
             }else{ if( strcmp(symb,"(") == 0 
                     ){ root = Expression(f,spp,ecp,epp); 
                          if (*ecp) return NULL; 
                          *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
                          if (strcmp(symb,")")!=0) return set_err(18,ecp); 
                     }else{ *spp = scanpos; 
                          return set_err(19,ecp); 
                  } 
          } 
  } 
  return root; 
} 
/******* ******* ******* ******* ******* ******* ******* ******* ********/ 

static struct treenode *FuncDesignator(char *f, int *spp, 
                                                int *ecp, int *epp) 
/*  FuncDesignator ::= FuncIdentifier "(" ActualParameterlist ")" 
 * 
 *       Identifier 
 *      /          \ 
 *  NULL            ActualParameterlist 
 */ 

{ 
  struct treenode *root, *r; 
  char            symb[SYMBLENGTH], fsymb[SYMBLENGTH]; 
  int             sk, fsk, t, ft, corr_panz; 
  double           fv = 0.0; 

  *ecp = nxt_symbol(f,spp,fsymb,&ft,&fsk,epp); 
  /*if (fsk != SFUNCS && fsk != UFUNCS) return set_err(20,ecp);*/ 
  if (fsk != SFUNCS) return set_err(20,ecp); 
  corr_panz = symbol[ft].info; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,"(") != 0) return set_err(21,ecp); 
  r = ActualParameterlist(f,spp,ecp,epp,corr_panz); if (*ecp) return NULL; 
  *ecp = nxt_symbol(f,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,")") != 0) return set_err(22,ecp); 
  root = create_node(fsymb,ft,fsk,fv,NULL,r); 
  return root; 
} 
/****** ****** ****** ****** ******* ******* ****** ****** ****** *******/ 

static struct treenode *ActualParameterlist(char *f, int *spp, int *ecp, 
                                           int *epp, int corr_panz) 
/*  ActualParameterlist ::= ActualParameter { "," ActualParameter } 
 * 
 *                                           "," 
 *  Expression  oder:                       /   \ 
 *                     more expressions tree     Expression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, ct, scanpos, panz; 
  double           v = 0.0; 

  root = Expression(f,spp,ecp,epp); if (*ecp) return NULL; 
  panz = 1; 
  scanpos = *spp; 
  *ecp = nxt_symbol(f,spp,symb,&ct,&sk,epp); 
  while (strcmp(symb,",") == 0) {
     panz++; 
     if (panz > corr_panz) return set_err(23,ecp); 
     l = root; 
     r = Expression(f,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(",",ct,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(f,spp,symb,&ct,&sk,epp); 
  }
  if (panz < corr_panz) return set_err(24,ecp); 
  *spp = scanpos; 
  return root; 
} 
/*********************** *********************** ************************/ 

static struct treenode *create_node(char *symb, int token, int symbkind,
                 double val, struct treenode *left, struct treenode *right) 

/*  Setzt im Knoten mit der Wurzel root den String, das Token und die Art 
 *  des Symbols (REL_OP, ADD_OP etc.), 
 */ 

{ 
  struct treenode *root; 
  int             sum, prod; 

  if( left != 0 && right != 0 && simplification(symb,token,left,right) 
     ){ root = left; 
          free((char*)right);   /* Speicher fuer Blatt wieder freigeben */ 
     }else{ root = tnalloc(); 
          strcpy(root->symb,symb); 
          root->token    = token; 
          root->symbkind = symbkind; 
          root->val      = val; 
          root->left     = left; 
          root->right    = right; 
  } 
  sum  = strcmp(root->symb,"+")==0 || strcmp(root->symb,"-")==0; 
  prod = strcmp(root->symb,"*")==0 || strcmp(root->symb,"/")==0; 
  if (sum || prod) check_reorg(root); 

  /* TIRLER */
  // printf("Adresse root:%p  -   Symbol: %s --- %s\n",root,symb,root->symb);


  return root; 
} 
/**************** ****************** ****************** *****************/ 

static int simplification(char *symb, int t, struct treenode *l, 
                                             struct treenode *r) 

/*  Untersucht, ob der Knoten vereinfacht werden kann; wenn ja, wird in 
 *  create_node kein neuer Knoten, sondern ein Blatt zurueckgegeben. 
 */ 

{ 
  int product  = strcmp(   symb,"*"  ) == 0; 
  int quotient = strcmp(   symb,"/"  ) == 0; 
  int power    = strcmp(   symb,"^"  ) == 0; 
  int plus     = strcmp(   symb,"+"  ) == 0; 
  int minus    = strcmp(   symb,"-"  ) == 0; 
  int and      = strcmp(   symb,"AND") == 0; 
  int mod      = strcmp(   symb,"MOD") == 0; 
  int l_0      = strcmp(l->symb,"0"  ) == 0; 
  int r_0      = strcmp(r->symb,"0"  ) == 0; 
  int l_1      = strcmp(l->symb,"1"  ) == 0; 
  int r_1      = strcmp(r->symb,"1"  ) == 0; 
  int leaves   = l->left==0 && l->right==0 && r->left==0 && r->right==0; 
  int eq_leaves= leaves && strcmp(l->symb,r->symb) == 0; 
  int l_const  = l->symbkind == SCONST || l->symbkind == UCONST; 
  int r_const  = r->symbkind == SCONST || r->symbkind == UCONST; 

  if (l_const && l_0) l_const = FALSE; /* 0-Blatt muss bleiben */ 
  /* --- Ueberpruefen, ob x/x, x^0 oder 1^x => Blatt = 1: --- */ 
  if( (quotient && eq_leaves) || (power && (r_0 || l_1)) ){ 
     strcpy(l->symb,"1"); 
     l->token    = scs + 2; 
     l->symbkind = SCONST; 
     l->val      = 1.0; 
     l->left     = NULL; 
     l->right    = NULL; 
     return TRUE; 
  } 

  /*-- Ueberpruefen, ob 0*x,x*0,0ANDx,xAND0,0/x,0^x,0MODx => Blatt=0: --*/ 
  if( ((product||and) && (l_0||r_0)) || (l_0 && (quotient||power||mod)) ){ 
     strcpy(l->symb,"0"); 
     l->token    = scs + 1; 
     l->symbkind = SCONST; 
     l->val      = 0.0; 
     l->left     = NULL; 
     l->right    = NULL; 
     return TRUE; 
  } 

  /*- Ueberpruefen, ob x+0,x-0,x*1,x/1,x MOD 1,x^1 => l. Baum zurueck:- */ 
  if( (r_0 && (plus||minus)) || (r_1 && (product||quotient||mod||power)) ){ 
     return TRUE; 
  } 

  /* --- Ueberpruefen, ob 0+x, 1*x => rechten Baum zurueck: --- */ 
  if( (l_0 && plus ) || (l_1 && product) ){ 
     strcpy(l->symb,r->symb); 
     l->token    = r->token; 
     l->symbkind = r->symbkind; 
     l->val      = r->val; 
     l->left     = r->left; 
     l->right    = r->right; 
     return TRUE; 
  } 


  /* --- Ueberpruefen, ob beide Blaetter Konstanten => ausrechnen: --- */ 
  if( l_const && r_const ){ 
     l->val = (*symbol[t].vcalc)(t,
                              symbol[l->token].val,symbol[r->token].val); 
     if( l->symbkind == SCONST && r->symbkind == SCONST && !quotient 
        ){ sprintf(l->symb,"%d",(int)l->val);/* Wert ist int         */ 
        }else{ sprintf(l->symb,"%g",l->val);     /* Wert ist evtl. double */ 
     } 
     l->token    = find_index(l->symb,scs,uce,'0'); 
     l->symbkind = find_kind(l->token); 
     l->left     = NULL; 
     l->right    = NULL; 
     return TRUE; 
  } 


  return FALSE; 
} 
/**************** ****************** ****************** *****************/ 

static struct treenode *tnalloc(void)/* Alloziert Speicher f. Baumknoten*/ 
{ 
  struct treenode *p; 

  p = ((struct treenode *) malloc(sizeof(struct treenode))); 
  return p; 
} 
/**************** ****************** ****************** *****************/ 

void check_reorg(struct treenode *root) 

/*  Ueberprueft, ob bei Addition, Subtraktion, Multiplikation oder 
 *  Division durch Umstellen der Knoten Vereinfachung moeglich ist und 
 *  fuehrt sie dann durch. 
 */ 

{ 
  struct treenode *temp; 
  int             mul     = strcmp(root->       symb,"*") == 0; 
  int             div     = strcmp(root->       symb,"/") == 0; 
  int             minus   = strcmp(root->       symb,"-") == 0; 
  int             plus    = strcmp(root->       symb,"+") == 0; 
  int             l_minus = strcmp(root->left-> symb,"-") == 0; 
  int             r_minus = strcmp(root->right->symb,"-") == 0; 
  int             r_plus  = strcmp(root->right->symb,"+") == 0; 

  if( (plus||minus) && r_minus&&strcmp(root->right->left->symb,"0")==0 ){ 
     /*  Baum1+(-Baum2) => Baum1-Baum2   | Baum1-(-Baum2) => Baum1+Baum2 
      *                                  | 
      *        +                 -       |       -                 + 
      *       / \               / \      |      / \               / \ 
      *  Baum1   -      => Baum1   Baum2 | Baum1   -      => Baum1   Baum2 
      *         / \                      |        / \ 
      *        0   Baum2                 |       0   Baum2 
      */ 
     if( plus 
        ){ strcpy(root->symb,"-"); 
        }else{ strcpy(root->symb,"+"); 
     } 
     free((char*)root->right->left);/* Speicher f."0"-Blatt wieder frei */ 
     temp        = root->right; 
     root->right = root->right->right; 
     free((char*)temp);         /* Speicher fuer "-"-Knoten wieder frei */ 
     goto cr_lbl; 
  } 

  if( plus && l_minus && strcmp(root->left->left->symb,"0") == 0 ){ 
     /* ---   -Baum1+Baum2 => Baum2-Baum1: ---
      * 
      *      +                  -
      *     / \                / \ 
      *    -   Baum2  =>  Baum2   Baum1 
      *   / \ 
      *  0   Baum1 
      */ 
     strcpy(root->symb,"-"); 
     free((char*)root->left->left);/* Speicher f."0"-Blatt  wieder frei */ 
     temp        = root->left->right; 
     free((char*)root->left); /* Speicher fuer f."-"-Knoten wieder frei */ 
     root->left  = root->right; 
     root->right = temp; 
     goto cr_lbl; 
  } 

  if( minus && (r_minus || r_plus) ){ 
     /*  B1-(B2-B3) => B1-B2+B3  |  B1-(B2+B3) => B1-B2-B3 
      *                          | 
      *     -              +     |     -               -
      *    / \            / \    |    / \             / \ 
      *  B1   -     =>   -   B3  |  B1   +     =>    -   B3 
      *      / \        / \      |      / \         / \ 
      *    B2   B3    B1   B2    |    B2   B3     B1   B2 
      */ 
     if( r_minus 
        ){ strcpy(root->symb,"+"); 
        }else{ strcpy(root->right->symb,"-"); 
             root->right->token    = find_index("-",aos,aoe,'0'); 
             root->right->symbkind = find_kind(root->right->token); 
     } 
     temp = root->right->right;              /* temp  = Baum3 */ 
     root->right->right = root->right->left; /* Baum3 = Baum2 */ 
     root->right->left  = root->left;        /* Baum2 = Baum1 */ 
     root->left         = root->right;       /* Baum1 = mBaum */ 
     root->right        = temp;              /* mBaum = Baum3 */ 
     goto cr_lbl; 
  } 

  if( (mul || div) && r_minus ){ 
     /*  B1*(-B2) => -B1*B2       |       B1/(-B2) => -B1/B2 
      * 
      *     *             -
      *    / \           / \ 
      *  B1   -     =>  0   * 
      *      / \           / \ 
      *     0   B2       B1   B2 
      */ 
     if( mul 
        ){ strcpy(root->right->symb,"*"); 
        }else{ strcpy(root->right->symb,"/"); 
     } 
     root->right->token    = find_index(root->right->symb,mos,moe,'0'); 
     root->right->symbkind = find_kind(root->right->token); 
     strcpy(root->symb,"-"); 
     temp = root->left;                      /* temp  = Baum1 */ 
     root->left        = root->right->left;  /* Baum1 = 0     */ 
     root->right->left = temp;               /*     0 = Baum1 */ 
     goto cr_lbl; 
  } 

  cr_lbl: 
  root->token    = find_index(root->symb,aos,moe,'0'); 
  root->symbkind = find_kind(root->token); 
  return; 
} 
/************************************************************************/ 

char *readln(char *s)       /* Eine Zeile Zeichen einlesen incl. Blancs */ 

{ 
  int c; 

  while ((c = getchar()) != '\n') 
       *(s++) = c; 
  *(s++) = '\0'; 
  return s; 
} 
/************************************************************************/ 

void pt_error(char *function, int errcode, int errpos) 

/*  Gibt Fehlermeldung und Fehlerposition nach Baumerzeugung aus. */ 

{ 
  int i = 0; 

  if( errcode != 0 ){ 
     printf("\n%s",function); 
     function[0] = '\0'; 
     while (i++ < errpos) strcat(function,"-"); 
     strcat(function,"^"); 
     printf("\n%s",function); 
  } 
  printf("\nErrcode %d: %-46.46s",errcode,errorstrings[errcode]); 
  if (errcode !=0) printf("at pos. %d.",errpos); 
} 
/************************************************************************/ 

void show_tree(struct treenode *root)
                                 /* Gibt rekursiv einen Parse-Baum aus. */ 
{ 
  static int count = 0; 
  int i; 

  printf("\n"); 
  count++; 
  for (i = 1; i < count; i++) printf("   "); 
  if( root != NULL 
     ){ printf("%s|root=%p|token=%d|sk=%d|val=%g|left=%p|right=%p",
          root->symb, root, root->token, root->symbkind, root->val,
          root->left, root->right); 
          show_tree(root->right); 
          show_tree(root->left); 

 /* TIRLER */
	  //  printf("Adresse root:%p  -   Symbol:    --- %s\n",root,root->symb);



     }else{ printf("NULL"); 
  } 
  count--; 
} 
/************************************************************************/ 

void _unur_fstr_free(struct treenode *root)  
                          /* Gibt Speicher fuer schon exist. Baum frei. */ 
{ 
  if( root != NULL ){
 
     _unur_fstr_free(root->right); 
     _unur_fstr_free(root->left); 
     free((char*)root); 
  } 
} 
/************************************************************************/ 

/************************************************************************
 * MODUL pars.c                                                         *
 *                                                                      *
 *                                                                      * 
 ************************************************************************/

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
  double          result;
  struct treenode *froot;
  int             ftok,xtok;

  xtok=find_index("X",uis,ufe,0);
  ftok=find_index("F",uis,ufe,0);
    froot=(*symbol[ftok].tree).right;             /* Achtung Fehler in Beschreibung !!! */
  //   froot=symbol[ftok].tree;  
  show_tree(E_root);
  symbol[xtok].val= argument;
  show_tree(E_root); 
 result=tree2float(froot);
  // show_tree(E_root);
  return result;
  }



/***************************************************************************************/
char *Ntree2string(struct treenode *tree_root, char *ret_str)

{
  struct treenode  *froot;
  int              ftok;

  ftok=find_index("F",uis,ufe,0);
  froot=(*symbol[ftok].tree).right;          

  tree2Cstring(froot,ret_str); 
  return ret_str;
}

/***************************************************************************************/
struct treenode *_unur_fstr_make_derivative(struct treenode *root)

{ 
   struct treenode *parsetreeh;
  
   char            x[MAXLENGTH];
   char            ret_str[MAXLENGTH];
  
   strcpy(x,"X"); 
   parsetreeh=part_deriv(root,x,ret_str);
   return parsetreeh;
}

/***************************************************************************************/
double  _unur_fstr_dev_eval_tree(struct treenode *E_root, double argument)
{
  double           result;
  struct treenode *froot;
  int             ftok,xtok;

  xtok=find_index("X",uis,ufe,0);
  ftok=find_index("F_X",uis,ufe,0);
  froot=(*symbol[ftok].tree).right;             /* Achtung Fehler in Beschreibung !!! */
  symbol[xtok].val= argument;
  result=tree2float(froot);
  return result;
  }























