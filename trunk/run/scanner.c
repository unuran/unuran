
#include "scanpars.h"

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
extern float v_exp(), v_dummy  (), v_or   (), v_mul(), v_sconst();
extern float v_ln (), v_less   (), v_xor  (), v_and(), v_uconst();
extern float v_log(), v_equal  (), v_plus (), v_div(), v_uident();
extern float v_sin(), v_greater(), v_minus(), v_mod();
extern float v_cos(), v_less_or();
extern float v_tan(), v_unequal(), v_power(), v_ufuncs();
extern float v_sec(), v_grtr_or(), v_not  ();
extern float v_sqr();
extern float v_abs();

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
  if (c >= '0' && c <= '9' || c == '.') {           /* UnsignedConstant */
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
  IF c == 'E' THEN /* wenn "E" folgt, wird noch ein ScaleFactor gelesen */
     strcat(uc,"E");
     (*scanposp)++;
     strcat(uc,get_sf(function,scanposp,s));
  ENDIF
  IF c == '.' THEN       /* wenn "." folgt, koennen noch Ziffern kommen */
     strcat(uc,".");
     (*scanposp)++;
     strcat(uc,get_ds(function,scanposp,s));
     c = function[*scanposp];
     IF c == 'E' THEN  /* wenn "E" folgt, wird noch ScaleFactor gelesen */
	strcat(uc,"E");
	(*scanposp)++;
	strcat(uc,get_sf(function,scanposp,s));
     ENDIF
  ENDIF
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
  IF (c=function[*scanposp]) == '+' || c == '-' THEN
     sf[0] = c;
     (*scanposp)++;
  ENDIF
  return strcat(sf,get_ds(function,scanposp,ds));
}
/*********************************** ************************************/

static int get_id_symbol(char *function, int *scanposp, char *id)

/*  Identifier ::= Letter { Letter | Digit }
 *  Letter     ::= "A" | "B" | ... | "Z" | "_"
 *  Digit      ::= "0" | "1" | ... | "9"
 */

{
  while ( (*id = function[(*scanposp)]) >= 'A' && *id <= 'Z'
     || *id == '_' ||               *id >= '0' && *id <= '9') {
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
  IF start == scs THEN
     /* --- Symbol ist Userdefined Constant: --- */
     nxt_free_place = ucs + symbol[ucs].info + 1;
     IF nxt_free_place >= uce
	THEN /* --- kein Platz mehr: --- */
	     return 0;
	ELSE strcpy(symbol[nxt_free_place].name,symb);
	     symbol[nxt_free_place].val = atof(symb);
	     symbol[ucs].info++;     /* hier Zeiger auf letzten Eintrag */
	     return nxt_free_place;
     ENDIF
  ENDIF
  IF start == aos THEN
     /* --- Symbol ist Identifier: --- */
     IF nxt_c == '('
	THEN /* --- Symbol ist Userdefined Function --- */
	     nxt_free_place = ufs + symbol[ufs].info + 1;
	     IF nxt_free_place >= ufe
		THEN /* --- kein Platz mehr: --- */
		     return 0;
		ELSE strcpy(symbol[nxt_free_place].name,symb);
		     symbol[nxt_free_place].info = 0;
		     symbol[ufs].info++;
		     return nxt_free_place;
	     ENDIF
	ELSE /* --- Symbol ist Userdefined Identifier --- */
	     nxt_free_place = uis + symbol[uis].info + 1;
	     IF nxt_free_place >= uie
		THEN /*** kein Platz mehr: ***/
		     return 0;
		ELSE strcpy(symbol[nxt_free_place].name,symb);
		     symbol[nxt_free_place].val = 0.0;
		     symbol[uis].info++;
		     return nxt_free_place;
	     ENDIF
     ENDIF
  ENDIF
  return 0;
}
