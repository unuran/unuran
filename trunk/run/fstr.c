
#include <ctype.h>
#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

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


/* --- Prototypen: --- */
void  show_symb_tab      (void);
void clear_symbol_area(int start, int end);
static void get_ds             (char *function, int  *scanpos, char *ds);
static void get_sf             (char *function, int  *scanpos, char *sf);
static int   get_uc_symbol      (char *function, int *scanpos, char *uc);
static int   nxt_symbol         (char function[], int *scanpos, char symb[],
				 int *tokenp, int *symbkindp, int *errposp);
static int   get_id_symbol      (char *function, int  *scanpos, char *id);
static int   get_ro_symbol      (char *function, int  *scanpos, char *ro);
static int   find_kind          (int token);
static int   find_index         (char *symb, int start, int end, int nxt_c);

/* --- Bereichs-Start- und -Endemarkierungen in der Symboltabelle: --- */
static int ros, roe, aos, aoe, mos, moe, hos, hoe, oss, ose, scs, sce;
static int ucs, uce, uis, uie, ufs, ufe, sfs, sfe;



/*---------------------------------------------------------------------------*/
/* functions for evaluating function tree                                    */
static double v_dummy  (int t, double l, double r);
static double v_less   (int t, double l, double r);
static double v_equal  (int t, double l, double r);
static double v_greater(int t, double l, double r);
static double v_less_or(int t, double l, double r);
static double v_unequal(int t, double l, double r);
static double v_grtr_or(int t, double l, double r);
static double v_or     (int t, double l, double r);
static double v_xor    (int t, double l, double r);
static double v_plus   (int t, double l, double r); 
static double v_minus  (int t, double l, double r); 
static double v_mul    (int t, double l, double r);
static double v_and    (int t, double l, double r);
static double v_div    (int t, double l, double r);
static double v_mod    (int t, double l, double r);
static double v_power  (int t, double l, double r);
static double v_not    (int t, double l, double r);
static double v_sconst (int t, double l, double r);
static double v_uconst (int t, double l, double r);
static double v_uident (int t, double l, double r);
static double v_ufuncs (int t, double l, double r);
static double v_exp    (int t, double l, double r);
static double v_ln     (int t, double l, double r);
static double v_log    (int t, double l, double r);
static double v_sin    (int t, double l, double r); 
static double v_cos    (int t, double l, double r);
static double v_tan    (int t, double l, double r);
static double v_sec    (int t, double l, double r); 
static double v_sqr    (int t, double l, double r);
static double v_abs    (int t, double l, double r);

/*---------------------------------------------------------------------------*/
/* functions for compting derivatives                                        */

/** TODO: incomplete prototype **/
static char *d_exp(), *d_dummy ();
static char *d_ln (), *d_add   ();
static char *d_log(), *d_mul   ();
static char *d_sin(), *d_div   ();
static char *d_cos(), *d_power ();
static char *d_tan(), *d_const ();
static char *d_sec();
/** static char *d_par   (); **/
static char *d_sqr(), *d_ufuncs();
static char *d_abs();


/* --- Import aus SCANNER.C: --- */ 
/** static void init_scanner(void); **/

/*long _STKSIZ = 500000; grosser Stack f. Rekursion in derive_expression*/

/************************************************************************/

/* --- Importe aus scanner.c und parser.c: --- */

static struct treenode *string2tree(char *function, int *errcodep, 
                                                    int *errposp);
static void pt_error    (char *function, int errcode, int errpos);
static void show_tree   (struct treenode *node, int level, int location);

char *readln     (char *s);

/************************************************************************/ 

/* --- Prototypen: --- */

static char *tree2string(struct treenode *tree_root, char *ret_str);
static char *tree2Cstring(struct treenode *tree_root, char *ret_str);
static char *strcpy3    (char *s1, char *s2, char *s3, char *s4); 
static char *strcpy5    (char *s1, char *s2, char *s3, char *s4, char *s5, char *s6); 
static double tree2float (struct treenode *E_root); 
/** static void gradient(struct treenode *root); **/
static void nxt_part_derivation(struct treenode *DP_root,struct treenode *root);
static struct treenode *part_deriv(struct treenode *root,char *par,char *ret_str); 
static char *derive_expression(struct treenode *E_root,char *par,char *ret_str);
static char *strcpy6(char *s1, char *s2, char *s3,
	      char *s4, char *s5, char *s6, char *s7);


static char *_unur_prepare_string( const char *str );
/*---------------------------------------------------------------------------*/
/* Prepare string for processing.                                            */
/*---------------------------------------------------------------------------*/



/************************************************************************/
/* --- SYMBOLTABELLE, mit Funktionsadressen --- */

static struct symbols symbol[] = {   
  {"_ROS",0, .0,v_dummy,  d_dummy, NULL},            /* RelationalOperators */
  {"<",   1, .0,v_less,   d_const, NULL},
  {"=",   1, .0,v_equal,  d_const, NULL},     /** TODO: == statt = **/
  {">",   1, .0,v_greater,d_const, NULL},
  {"<=",  1, .0,v_less_or,d_const, NULL},
  {"<>",  1, .0,v_unequal,d_const, NULL},
  {">=",  1, .0,v_grtr_or,d_const, NULL},
  {"_ROE",0, .0,v_dummy,  d_dummy, NULL},

  {"_AOS",0, .0,v_dummy,  d_dummy, NULL},                /* AddingOperators */
  {"or",  2, .0,v_or,     d_const, NULL},
  {"xor", 2, .0,v_xor,    d_const, NULL},
  {"+",   2, .0,v_plus,   d_add,   NULL},
  {"-",   3, .0,v_minus,  d_add,   NULL},
  {"_AOE",0, .0,v_dummy,  d_dummy, NULL},

  {"_MOS",0, .0,v_dummy,  d_dummy, NULL},           /* MultiplyingOperators */
  {"*",   4, .0,v_mul,    d_mul,   NULL},
  {"and", 4, .0,v_and,    d_const, NULL},
  {"/",   4, .0,v_div,    d_div,   NULL},
  {"mod", 4, .0,v_mod,    d_const, NULL},
  {"_MOE",4, .0,v_dummy,  d_dummy, NULL},

  {"_HOS",0, .0,v_dummy,  d_dummy, NULL},          /* HighPriorityOperators */
  {"^",   5, .0,v_power,  d_power, NULL},
  {"not", 6, .0,v_not,    d_const, NULL},
  {"_HOE",0, .0,v_dummy,  d_dummy, NULL},

  {"_OSS",0, .0,v_dummy,  d_dummy, NULL},                   /* OtherSymbols */
  {"(",   0, .0,v_dummy,  d_dummy, NULL},
  {")",   0, .0,v_dummy,  d_dummy, NULL},
  {",",   0, .0,v_dummy,  d_dummy, NULL},
  {"_OSE",0, .0,v_dummy,  d_dummy, NULL},

  {"_SCS",0, .0,v_dummy,  d_dummy, NULL},                /* SystemConstants */
  {"0",   0,0.0,v_sconst, d_const, NULL},
  {"1",   0,1.0,v_sconst, d_const, NULL},
/*    {"2",  0,2.0,v_sconst, d_const, NULL}, */
/*    {"3",  0,3.0,v_sconst, d_const, NULL}, */
/*    {"4",  0,4.0,v_sconst, d_const, NULL}, */
/*    {"5",  0,5.0,v_sconst, d_const, NULL}, */
/*    {"6",  0,6.0,v_sconst, d_const, NULL}, */
/*    {"7",  0,7.0,v_sconst, d_const, NULL}, */
/*    {"8",  0,8.0,v_sconst, d_const, NULL}, */
/*    {"9",  0,9.0,v_sconst, d_const, NULL}, */
  {"pi",  0, M_PI,v_sconst, d_const, NULL},
  {"e",   0, M_E ,v_sconst, d_const, NULL},
  {"_SCE",0, .0,v_dummy,  d_dummy, NULL},

  {"_UCS",0, .0,v_dummy,  d_dummy, NULL},           /* UserdefinedConstants */
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"",    0, .0,v_uconst, d_const, NULL},
  {"_UCE",0, .0,v_dummy,  d_dummy, NULL},

  {"_UIS",0, .0,v_dummy,  d_dummy, NULL},         /* UserdefinedIdentifiers */
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"",    0, .0,v_uident, d_const, NULL},
  {"_UIE",0, .0,v_dummy,  d_dummy, NULL},

  {"_UFS",0, .0,v_dummy,  d_dummy, NULL},           /* UserdefinedFunctions */
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"",    0, .0,v_ufuncs, d_ufuncs,NULL},
  {"_UFE",0, .0,v_dummy,  d_dummy, NULL},

  {"_SFS",0, .0,v_dummy,  d_dummy, NULL},                /* SystemFunctions */
  {"exp", 1, .0,v_exp,    d_exp,   NULL},
  {"ln",  1, .0,v_ln,     d_ln,    NULL},
  {"log", 2, .0,v_log,    d_log,   NULL},
  {"sin", 1, .0,v_sin,    d_sin,   NULL},
  {"cos", 1, .0,v_cos,    d_cos,   NULL},
  {"tan", 1, .0,v_tan,    d_tan,   NULL},
  {"sec", 1, .0,v_sec,    d_sec,   NULL},
  {"sqr", 1, .0,v_sqr,    d_sqr,   NULL},
  {"abs", 1, .0,v_abs,    d_abs,   NULL},
  {"_SFE",0, .0,v_dummy,  d_dummy, NULL},

  {"_END",0, .0,v_dummy,  d_dummy, NULL},
};

/* --- Prototypen: --- */ 
static struct treenode *FuncDefinition   (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *DefFuncDesignator(char *fstr, int *scanpos, int *ecp, int *epp, int *ftp);
static struct treenode *DefParameterlist (char *fstr, int *scanpos, int *ecp, int *epp, int *paranzp);
static struct treenode *Expression       (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *SimpleExpression (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *VTerm            (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *Term             (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *Factor           (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *bas_exp          (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *FuncDesignator   (char *fstr, int *scanpos, int *ecp, int *epp);
static struct treenode *ActualParameterlist(char *fstr, int *scanpos, int *ecp, int *epp, int corr_panz);

static struct treenode *create_node(char *symb, int token, int symbkind,
                double val, struct treenode *left, struct treenode *right);

static struct treenode *set_err (int err_nr, int *errcodep);
char            *readln  (char *s);
static int             simplification(char *symb, int t, struct treenode *l, 
                                                      struct treenode *r);
static void            check_reorg (struct treenode *root);
               
static char *errorstrings[] = { 
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
     if         (strcmp(s,"_ROS") == 0) ros = i;
	else if (strcmp(s,"_ROE") == 0) roe = i;
	else if (strcmp(s,"_AOS") == 0) aos = i;
	else if (strcmp(s,"_AOE") == 0) aoe = i;
	else if (strcmp(s,"_MOS") == 0) mos = i;
	else if (strcmp(s,"_MOE") == 0) moe = i;
	else if (strcmp(s,"_HOS") == 0) hos = i;
	else if (strcmp(s,"_HOE") == 0) hoe = i;
	else if (strcmp(s,"_OSS") == 0) oss = i;
	else if (strcmp(s,"_OSE") == 0) ose = i;
	else if (strcmp(s,"_SCS") == 0) scs = i;
	else if (strcmp(s,"_SCE") == 0) sce = i;
	else if (strcmp(s,"_UCS") == 0) ucs = i;
	else if (strcmp(s,"_UCE") == 0) uce = i;
	else if (strcmp(s,"_UIS") == 0) uis = i;
	else if (strcmp(s,"_UIE") == 0) uie = i;
	else if (strcmp(s,"_UFS") == 0) ufs = i;
	else if (strcmp(s,"_UFE") == 0) ufe = i;
	else if (strcmp(s,"_SFS") == 0) sfs = i;
	else if (strcmp(s,"_SFE") == 0) sfe = i;
  } while (!(strcmp(s,"_END") == 0));
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
  } while (!(strcmp(symbol[++i].name,"_END") == 0));
  printf("\nros=%4d hos=%4d ucs=%4d sfs=%4d",ros,hos,ucs,sfs);
  printf("\nroe=%4d hoe=%4d uce=%4d sfe=%4d",roe,hoe,uce,sfe);
  printf("\naos=%4d oss=%4d uis=%4d"        ,aos,oss,uis    );
  printf("\naoe=%4d ose=%4d uie=%4d"        ,aoe,ose,uie    );
  printf("\nmos=%4d scs=%4d ufs=%4d"        ,mos,scs,ufs    );
  printf("\nmoe=%4d sce=%4d ufe=%4d"        ,moe,scs,ufe    );
  printf("\n");
}
/************************************************************************/

int nxt_symbol (char *fstr, int *scanpos, char *symb, int *tokenp,
		int *symbkindp, int *errposp)

/*  Liefert aus dem String 'fstr' das naechste Symbol ab der
 *  Position 'scanpos' und speichert es in 'symb'. 
 *  Nach der Bearbeitung zeigt *scanpos auf das Zeichen, 
 *  das dem ermittelten Symbol unmittelbar folgt. Leerzeichen
 *  werden ignoriert. In *tokenp wird der Index des Symbols in der
 *  Symboltabelle. geliefert. *symbkindp ist die Art des Symbols
 *  (UIDENT, SCONST etc.), *errposp zeigt auf das Zeichen, auf das
 *  *scanpos beim Aufruf zeigte.
 *  Der zurueckgelieferte Funktionswert ist der Fehlercode.
 */

{
  char c;
  int  dummy = 0, nxt_c, errcode = 0;
  
  *errposp = *scanpos;
  
  if (*scanpos > strlen(fstr)) return 1;
  
  c = fstr[*scanpos];
  
  if ( (c >= '0' && c <= '9') || c == '.') {           /* UnsignedConstant */
    errcode = get_uc_symbol(fstr,scanpos,symb);
    if ((*tokenp = find_index(symb,scs,uce,dummy)) <= 0) return 2; }

  else if (c >=  'a' && c <= 'z') {                       /* Identifier */
    errcode = get_id_symbol(fstr,scanpos,symb);
    nxt_c = fstr[*scanpos];
    if ((*tokenp = find_index(symb,aos,sfe,nxt_c)) <= 0) return 3; }

  else if ( c == '<' || c == '>' ) {              /* relation symbol */
    errcode = get_ro_symbol(fstr,scanpos,symb);
    *tokenp = find_index(symb,ros,roe,dummy);                      }

  else {
    symb[0] = c; symb[1] = '\0';           /* alle anderen Zeichen */
    (*scanpos)++;
    *tokenp = find_index(symb,ros,sfe,dummy);
  }

  *symbkindp = find_kind(*tokenp);         /* Art des Symbols ermitteln */
  if (*symbkindp == 0) errcode = 4;        /*       unbekanntes Zeichen */

  return errcode;
}

/*********************************** ************************************/

/*********************************** ************************************/

static int get_uc_symbol(char *fstr, int *scanpos, char *uc)
     /* get unsigned constant */
/* UnsignedConstant ::= UnsignedInteger | UnsignedReal
 * UnsignedInteger  ::= DigitSequence
 * UnsignedReal ::= UnsignedInteger '.' DigitSequence ['E' ScaleFactor] |
 *                  UnsignedInteger 'E' ScaleFactor
 */

{
  /* store scan position */
  int startpos = *scanpos;

  /* copy digit sequence into uc */
  get_ds(fstr,scanpos,uc);

  if( fstr[*scanpos] == '.' ) {
    /* decimal point  --> copy point and following digits */
    *(uc + (*scanpos-startpos)) = '.';
    (*scanpos)++;
    get_ds(fstr,scanpos,uc+(*scanpos-startpos));
  }

  if( fstr[*scanpos] == 'e' ) {
    /* exponent --> copy indicator 'E' and following [sign and] digits */
    *(uc + (*scanpos-startpos)) = 'e';
    (*scanpos)++;
    get_sf(fstr,scanpos,uc+(*scanpos-startpos));
  }

  return 0;
}

/*********************** *********************** ************************/

static void get_ds(char *fstr, int *scanpos, char *ds)
     /* get digital sequence */
/*  DigitSequence   ::= Digit [ Digit ]
 *  Digit           ::= '0' | '1' | '2' | ... | '8' | '9'
 */

{
  /* copy digit */
  while ( (*ds = fstr[*scanpos]) >= '0' && *ds <= '9' ) {
     ds++;
     (*scanpos)++;
  }
  /* terminate string */
  *ds = '\0';
}

/*********************** *********************** ************************/

static void get_sf(char *fstr, int *scanpos, char *sf)
     /* get ??????  */
/*  ScaleFactor ::= [Sign] DigitSequence
 *  Sign        ::= '+' | '-'
 */

{
  /* copy sign */
  if( (sf[0]=fstr[*scanpos]) == '+' || sf[0] == '-' ){
     sf++;
     (*scanpos)++;
  }
  /* copy digital sequence (after sign) */
  get_ds(fstr,scanpos,sf);
}

/*********************************** ************************************/

static int get_id_symbol(char *fstr, int *scanpos, char *id)

/*  Identifier ::= Letter [ Letter | Digit ]
 *  Letter     ::= 'a' | 'b' | ... | 'z' | '_'
 *  Digit      ::= '0' | '1' | ... | '9'
 */

{
  /* copy word */
  while ( ((*id = fstr[(*scanpos)]) >= 'a' && *id <= 'z')
	  || *id == '_' 
	  || ( *id >= '0' && *id <= '9')) {
    id++;
    (*scanpos)++;
  }
  /* terminate string */
  *id = '\0';

  return 0;
}

/*********************************** ************************************/

static int get_ro_symbol(char *fstr, int *scanpos, char *ro)
     /*  *le='<','<=','<>'. */   /*  *ge= '>' oder '>='.*/
{
  /* copy relation operator */
  while ((*ro = fstr[*scanpos]) == '<' || *ro == '=' || *ro == '>') {
    ro++;
    (*scanpos)++;
  }
  /* terminate string */
  *ro = '\0';

  return 0;
}
/*********************************** ************************************/

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
 /** printf("\n %d: %s - %s -> %s",i,symb, symbol[i].name, symbol[64].name); **/
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
 
/************************************************************************/ 

struct treenode *
string2tree (char *functstr, int *errcodep, int *errposp) 

/*  Wandelt den String 'functstr' in einen Parse-Baum um. 'functstr' muss 
 *  in der Form f(a,b,...)=Expression vorliegen. In *errcodep und *errposp 
 *  wird im Fehlerfall eine Fehlernummer und die -position im String 
 *  'functstr' zurueckgegeben, sonst ist *errcodep Null. 
 */ 

{ 
  struct treenode *root;
  char *fstr;
  int scanpos = 0; 

  *errcodep = 0; 

  fstr = _unur_prepare_string(functstr);

  root = FuncDefinition(fstr, &scanpos, errcodep, errposp); 

  if (scanpos != strlen(fstr)) {
    free (fstr);
    return set_err(11,errcodep); 
  }

  free (fstr);
  return root; 
} 

/*********************************** ************************************/ 

static struct treenode *
FuncDefinition (char *fstr,int *scanpos,int *ecp,int *epp) 

/*  FuncDefinition ::= DefFuncDesignator '=' Expression 
 * 
 *                    '=' 
 *                   /   \ 
 *  DefFuncDesignator     Expression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, ft; 
  double           v = 0.0; 

  l = DefFuncDesignator(fstr,scanpos,ecp,epp,&ft);
  if (*ecp) return NULL; 

  *ecp = nxt_symbol(fstr,scanpos,symb,&t,&sk,epp); 

  if (strcmp(symb,"=") != 0) return set_err(12,ecp); 
  r = Expression(fstr,scanpos,ecp,epp);
  if (*ecp) return NULL; 
  root = create_node(symb,t,sk,v,l,r); 
  symbol[ft].tree = root; 
  return root; 
} 

/*********************** *********************** ************************/ 

static struct treenode *
DefFuncDesignator (char *fstr, int *spp, int *ecp, int *epp, int *ftp) 

/*  DefFuncDesignator ::= Identifier '(' DefParameterlist ')' 
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

  *ecp = nxt_symbol(fstr,spp,fsymb,ftp,&fsk,epp); 
  if (fsk != UFUNCS) return set_err(13,ecp); 
  if( symbol[*ftp].info != 0 ){ 
     /* --- Funktion war schon vorhanden => alten Baum loeschen: --- */ 
     _unur_fstr_free(symbol[*ftp].tree); 
  } 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,"(") != 0) return set_err(14,ecp); 

  r = DefParameterlist(fstr,spp,ecp,epp,&paranz); if (*ecp) return NULL; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,")") != 0) return set_err(15,ecp); 

  symbol[*ftp].info = paranz; 
  root = create_node(fsymb,*ftp,fsk,fv,NULL,r); 
  return root; 
} 

/**************** ****************** ****************** *****************/ 

static struct treenode *
DefParameterlist(char *fstr, int *spp, int *ecp, int *epp, int *paranzp) 

/*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')' 
 * 
 *       Identifier                                    ',' 
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

  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if (sk != UIDENT) return set_err(16,ecp); 
  *paranzp = 1; 
  root = create_node(symb,t,sk,v,NULL,NULL); 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&ct,&sk,epp); 
  while (strcmp(symb,",") == 0) {
     l = root; 
     *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
     if (sk != UIDENT) return set_err(17,ecp); 
     r = create_node(symb,t,sk,v,NULL,NULL); 
     (*paranzp)++; 
     root = create_node(",",ct,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(fstr,spp,symb,&ct,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 

/*********************** *********************** ************************/ 

static struct treenode *
set_err (int err_nr, int *errcodep) 

/*  Setzt *errcodep gleich err_nr und gibt NULL zurueck. */ 

{ 
  printf("\n%d %s",err_nr,errorstrings[err_nr]); 
  if (*errcodep == 0) *errcodep = err_nr;   /* nur den 1. Fehler melden */ 
  return NULL; 
} 

/*********************** *********************** ************************/ 

static struct treenode *
Expression(char *fstr, int *spp, int *ecp, int *epp) 

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

  l = SimpleExpression(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if( sk == REL_OP 
     ){ r = SimpleExpression(fstr,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ *spp = scanpos; 
          root = l; 
  } 
  return root; 
} 

/**************** ****************** ****************** *****************/ 

static struct treenode *
SimpleExpression(char *fstr, int *spp, int *ecp, int *epp) 

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

  root = VTerm(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  while (sk == ADD_OP) {
     l = root; 
     r = Term(fstr,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(symb,t,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 

/************ ************** *************** ************** *************/ 

static struct treenode *
VTerm(char *fstr, int *spp, int *ecp, int *epp) 

/*  Vterm ::= [ '+' | '-' ] Term 
 * 
 *                        '-' 
 *  Term  oder:          /   \ 
 *                    '0'     Term 
 *                   /   \ 
 *               NULL     NULL 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos = *spp; 
  double           v = 0.0; 

  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if( strcmp(symb,"-") == 0 
     ){ /* --- Term hat neg. Vorzeichen =>     --- */ 
          /* --- Vorzeichen in Operator wandeln: --- */ 
          l = create_node("0",scs+1,SCONST,0.0,NULL,NULL); 
          r = Term(fstr,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ /* --- Term hat pos. oder kein Vorzeichen: --- */ 
          if (strcmp(symb,"+") != 0) *spp = scanpos;  /* "+" ignorieren */ 
          root = Term(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  } 
  return root; 
} 

/********** *********** ************ ************ *********** ***********/ 

static struct treenode *
Term(char *fstr, int *spp, int *ecp, int *epp) 

/*  Term ::= Factor [ MultiplyingOperator Factor ]
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

  root = Factor(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  while (sk == MUL_OP) {
     l = root; 
     r = Factor(fstr,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(symb,t,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  }
  *spp = scanpos; 
  return root; 
} 
/******** ********* ********** *********** ********** ********* *********/ 

static struct treenode *
Factor(char *fstr, int *spp, int *ecp, int *epp) 

/*  Factor ::= Base [ '^' Exponent ] 
 * 
 *                          '^' 
 *  bas_exp  oder:         /   \ 
 *                  bas_exp     bas_exp 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos; 
  double           v = 0.0; 

  l = bas_exp(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if( strcmp(symb,"^") == 0 
     ){ r = bas_exp(fstr,spp,ecp,epp); if (*ecp) return NULL; 
          root = create_node(symb,t,sk,v,l,r); 
     }else{ *spp = scanpos; 
          root = l; 
  } 
  return root; 
} 

/******* ******** ******** ********* ********* ******** ******** ********/ 

static struct treenode *
bas_exp(char *fstr, int *spp, int *ecp, int *epp) 

/*  Base ::= Exponent ::= UnsignedConstant | Identifier | FuncDesignator | 
 *                        "not" Base | '(' Expression ')' 
 * 
 *       UnsignedConstant                Identifier 
 *      /                \     oder     /          \     oder 
 *  NULL                  NULL      NULL            NULL 
 * 
 *                                         "not" 
 *        FuncDesignator       oder       /     \        oder Expression 
 *                                    NULL       bas_exp 
 */ 

{ 
  struct treenode *root, *r; 
  struct treenode *bas_exp(); 
  char            symb[SYMBLENGTH]; 
  int             sk, t, scanpos = *spp; 
  double           v = 0.0; 

  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if( sk==SCONST || sk==UCONST || sk==UIDENT || strcmp(symb,"not")==0 
     ){ /* --- neuen Knoten bilden: --- */ 
          if( strcmp(symb,"not") == 0 
             ){ r = bas_exp(fstr,spp,ecp,epp); if (*ecp) return NULL; 
             }else{ r = NULL; 
          } 
          root = create_node(symb,t,sk,v,NULL,r); 
     }else{ /* --- neuer Knoten nicht noetig: --- */ 
          if( sk == UFUNCS || sk == SFUNCS 
             ){ *spp = scanpos; 
                  root = FuncDesignator(fstr,spp,ecp,epp); 
                  if (*ecp) return NULL; 
             }else{ if( strcmp(symb,"(") == 0 
                     ){ root = Expression(fstr,spp,ecp,epp); 
                          if (*ecp) return NULL; 
                          *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
                          if (strcmp(symb,")")!=0) return set_err(18,ecp); 
                     }else{ *spp = scanpos; 
                          return set_err(19,ecp); 
                  } 
          } 
  } 
  return root; 
} 

/******* ******* ******* ******* ******* ******* ******* ******* ********/ 

static struct treenode *
FuncDesignator(char *fstr, int *spp, int *ecp, int *epp) 
/*  FuncDesignator ::= FuncIdentifier '(' ActualParameterlist ')' 
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

  *ecp = nxt_symbol(fstr,spp,fsymb,&ft,&fsk,epp); 
  /*if (fsk != SFUNCS && fsk != UFUNCS) return set_err(20,ecp);*/ 
  if (fsk != SFUNCS) return set_err(20,ecp); 
  corr_panz = symbol[ft].info; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,"(") != 0) return set_err(21,ecp); 
  r = ActualParameterlist(fstr,spp,ecp,epp,corr_panz); if (*ecp) return NULL; 
  *ecp = nxt_symbol(fstr,spp,symb,&t,&sk,epp); 
  if (strcmp(symb,")") != 0) return set_err(22,ecp); 
  root = create_node(fsymb,ft,fsk,fv,NULL,r); 
  return root; 
} 

/****** ****** ****** ****** ******* ******* ****** ****** ****** *******/ 

static struct treenode *
ActualParameterlist(char *fstr, int *spp, int *ecp, int *epp, int corr_panz) 
/*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ] 
 * 
 *                                           ',' 
 *  Expression  oder:                       /   \ 
 *                     more expressions tree     Expression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             sk, ct, scanpos, panz; 
  double           v = 0.0; 

  root = Expression(fstr,spp,ecp,epp); if (*ecp) return NULL; 
  panz = 1; 
  scanpos = *spp; 
  *ecp = nxt_symbol(fstr,spp,symb,&ct,&sk,epp); 
  while (strcmp(symb,",") == 0) {
     panz++; 
     if (panz > corr_panz) return set_err(23,ecp); 
     l = root; 
     r = Expression(fstr,spp,ecp,epp); if (*ecp) return NULL; 
     root = create_node(",",ct,sk,v,l,r); 
     scanpos = *spp; 
     *ecp = nxt_symbol(fstr,spp,symb,&ct,&sk,epp); 
  }
  if (panz < corr_panz) return set_err(24,ecp); 
  *spp = scanpos; 
  return root; 
} 
/*********************** *********************** ************************/ 

static struct treenode *
create_node (char *symb, int token, int symbkind,
	     double val, struct treenode *left, struct treenode *right) 

/*  Setzt im Knoten mit der Wurzel root den String, das Token und die Art 
 *  des Symbols (REL_OP, ADD_OP etc.), 
 */ 

{ 
  struct treenode *root; 

  if ( left != NULL && right != NULL 
       && simplification(symb,token,left,right) ) { 
    /* node has been simplified  -->  use left node, remove right node */
    root = left;     
    free(right);   /* free memory for leave */
  }

  else {
    /* make new node */
    root = _unur_malloc(sizeof(struct treenode)); 
    strcpy(root->symb,symb); 
    root->token    = token; 
    root->symbkind = symbkind; 
    root->val      = val; 
    root->left     = left; 
    root->right    = right; 
  } 

  /* try to reorganize tree */
  if ( ( root->symb[0] == '+' ||
	 root->symb[0] == '-' ||
	 root->symb[0] == '*' ||
	 root->symb[0] == '/' )
       && root->symb[1] == '\0' )
    check_reorg(root); 

  /* return node */
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
  int and      = strcmp(   symb,"and") == 0; 
  int mod      = strcmp(   symb,"mod") == 0; 
  int l_0      = strcmp(l->symb,"0"  ) == 0; 
  int r_0      = strcmp(r->symb,"0"  ) == 0; 
  int l_1      = strcmp(l->symb,"1"  ) == 0; 
  int r_1      = strcmp(r->symb,"1"  ) == 0; 
  int leaves   = l->left==NULL && l->right==NULL && r->left==NULL && r->right==NULL; 
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

void check_reorg(struct treenode *root) 

/*  Ueberprueft, ob bei Addition, Subtraktion, Multiplikation oder 
 *  Division durch Umstellen der Knoten Vereinfachung moeglich ist und 
 *  fuehrt sie dann durch. 
 */ 

{ 

  /** TODO: root->left != 0 und root->right != NULL  !!!!! **/


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
     free(root->right->left);/* Speicher f."0"-Blatt wieder frei */ 
     temp        = root->right; 
     root->right = root->right->right; 
     free(temp);         /* Speicher fuer "-"-Knoten wieder frei */ 
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
     free(root->left->left);/* Speicher f."0"-Blatt  wieder frei */ 
     temp        = root->left->right; 
     free(root->left); /* Speicher fuer f."-"-Knoten wieder frei */ 
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

void _unur_fstr_debug_tree( struct treenode *root )
{
  printf("parse tree:  (left nodes above right nodes)\n\n"); 

  show_tree(root,0,0);
  printf("\n");

} /* end of _unur_fstr_debug_tree() */

/************************************************************************/ 

void show_tree(struct treenode *node, int level, int location)
                                 /* Gibt rekursiv einen Parse-Baum aus. */ 

     /* node ... pointer to node */
     /* level ... level of node in tree*/
     /* location ... indicates location of node in tree by bit array */
{ 
  int i, mask; 

  /* draw lines for tree */
  for (i = 0, mask = 1; i < level; i++, mask <<= 1) 
    if (mask & location) 
      printf("|   "); 
    else 
      printf("    "); 

  /* print node */
  if( node != NULL ) {
    printf("+--'%s'", node->symb);

    printf("\t\t|token=%d|sk=%d|val=%g",
	   node->token, node->symbkind, node->val );

    printf("\n"); 

    /* print left and right node */
    if ( node->left || node->right) {
      /* ... unless both leaves are empty */
      show_tree(node->left,level+1,location|(mask<<1)); 
      show_tree(node->right,level+1,location); 
    }
  }

  else {  /* empty leave */
    printf("+--(void)\n"); 
  } 

} 

/************************************************************************/ 

void _unur_fstr_free(struct treenode *root)  
                          /* Gibt Speicher fuer schon exist. Baum frei. */ 
{ 
  if( root != NULL ) {
    _unur_fstr_free(root->right); 
    _unur_fstr_free(root->left); 
    free(root); 
  } 
} /* end of _unur_fstr_free() */

/************************************************************************/ 

/************************************************************************
 * MODUL pars.c                                                         *
 *                                                                      *
 *                                                                      * 
 ************************************************************************/


/************************************************************************/ 


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

  

  /** TODO: initialize to avoid warning ???? **/
  t_m = t_l = t_r = 0;
  p_m = p_l = p_r = 0;
  sk_m = sk_l = sk_r = 0;


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
	if( strcmp(symb_m,"/") == 0 || strcmp(symb_m,"mod") == 0 ||
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
     if( strcmp(symb_m,"xor") == 0 || strcmp(symb_m,"and") == 0 ||
	strcmp(symb_m,"or" ) == 0 || strcmp(symb_m,"mod") == 0 ){
	 strcpy3(temp_str," ",symb_m," ");
	 strcpy(symb_m,temp_str);
     }

     /* --- Rueckgabestring zusammensetzen aus linkem String,
			  Verknuepfungsoperator und rechtem String: --- */
     return strcpy3(ret_str,str_l,symb_m,str_r);
  }

  /*-- bei "not" klammern, wenn Operator geringerer Prioritaet folgt: --*/
  if( strcmp(symb_m,"not") == 0 ){
     if( sk_r==REL_OP || sk_r==ADD_OP||sk_r==MUL_OP||strcmp(symb_r,"^")==0
	){ strcpy3(ret_str,"not(",str_r,")");
	     return ret_str;
	}else{ /* --- sonst nicht klammern: --- */
	     return strcpy3(ret_str,"not ",str_r,"");
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

  /** TODO: initialize to avoid warning ???? **/
  t_m = t_l = t_r = 0;
  p_m = p_l = p_r = 0;
  sk_m = sk_l = sk_r = 0;




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
	if( strcmp(symb_m,"/") == 0 || strcmp(symb_m,"mod") == 0 ||
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
     if( strcmp(symb_m,"xor") == 0 || strcmp(symb_m,"and") == 0 ||
	strcmp(symb_m,"or" ) == 0 || strcmp(symb_m,"mod") == 0 ){
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

  /*-- bei "not" klammern, wenn Operator geringerer Prioritaet folgt: --*/
  if( strcmp(symb_m,"not") == 0 ){
     if( sk_r==REL_OP || sk_r==ADD_OP||sk_r==MUL_OP||strcmp(symb_r,"^")==0
	){ strcpy3(ret_str,"not(",str_r,")");
	     return ret_str;
	}else{ /* --- sonst nicht klammern: --- */
	     return strcpy3(ret_str,"not ",str_r,"");
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

  /** TODO: initialize to avoid warning ???? **/
  val_l = val_r = 0;

  /* --- Werte der Wurzel: --- */
  if (E_root->left  != NULL)  val_l = tree2float(E_root->left);
  if (E_root->right != NULL)  val_r = tree2float(E_root->right);
  /**  printf("ergebnis: %f - %f - %s - %f \n",val_l, val_r,
       symbol[E_root->token].name , symbol[E_root->token].val      ); **/
  return (*symbol[E_root->token].vcalc)(E_root->token,val_l,val_r);
  
}
/*********************************** ************************************/

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
/** TODO: verschachtelte userdefinierte Funktionen noch nicht implementiert: **/
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

/************************************************************************
 * Prozeduren zur analytischen Ableitung eines Parse-Baumes:            *
 ************************************************************************/

#if 0
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
#endif

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
		       strcpy6(s,l,"^",r,"*(",dr,"*ln(");
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
	   { return strcpy3(s,dr,"*exp",r); }

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
	     return strcpy6(s,dr,"/(",r,"*ln(",l,"))");
	   }

char *d_sin(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SIN(r))' = r'*COS(r)        */
	   { return strcpy3(s,dr,"*cos",r); }

char *d_cos(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (COS(r))' = -r'*SIN(r)       */
	   { return strcpy6(s,"-",dr,"*sin",r,"",""); }

char *d_tan(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (TAN(r))' = r'*(SEC(r))^2    */
	   { return strcpy6(s,dr,"*(sec",r,")^2","",""); }

char *d_sec(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SEC(r))' = r'*TAN(r)*SEC(r) */
	   { return strcpy6(s,dr,"*tan",r,"*sec",r,""); }

char *d_sqr(char *par,struct treenode *w,char *l,char *r,char *dl,
	    char *dr,char *s)           /* (SQR(r))' = r'/(2*SQR(r))    */
	   { return strcpy6(s,dr,"/(2*sqr(",r,"))" ,"",""); }

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



/***************************************************************************************/
struct treenode *_unur_fstr2tree(char *function, int *errcodep, int *errposp)

{ char            *fvonx;
  struct treenode *root;

 
  fvonx= malloc(MAXLENGTH*sizeof(char));     
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

  xtok=find_index("x",uis,ufe,0);
  ftok=find_index("f",uis,ufe,0);
  froot=(*symbol[ftok].tree).right;             /* Achtung Fehler in Beschreibung !!! */
  /**   froot=symbol[ftok].tree;  **/
  symbol[xtok].val= argument;
  result=tree2float(froot);
  return result;
  }
/***************************************************************************************/
char *Ntree2string(struct treenode *tree_root, char *ret_str)

{
  struct treenode  *froot;
  int              ftok;

  ftok=find_index("f",uis,ufe,0);
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
  
   strcpy(x,"x"); 
   parsetreeh=part_deriv(root,x,ret_str);
   return parsetreeh;
}

/***************************************************************************************/
double  _unur_fstr_dev_eval_tree(struct treenode *E_root, double argument)
{
  double           result;
  struct treenode *froot;
  int             ftok,xtok;

  xtok=find_index("x",uis,ufe,0);
  ftok=find_index("f_x",uis,ufe,0);
  froot=(*symbol[ftok].tree).right;            
  symbol[xtok].val= argument;
  result=tree2float(froot);
  return result;
  }


/*****************************************************************************/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

char *
_unur_prepare_string( const char *str )
     /*----------------------------------------------------------------------*/
     /* Prepare string for processing:                                       */
     /*   Copy into working strint.                                          */
     /*   Remove all white spaces and convert to lower case letters.         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   str      ... pointer to string                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to working string.                                         */
     /*                                                                      */
     /* as a side effect, a new string is allocated.                         */
     /*----------------------------------------------------------------------*/
{
  char *tmp, *ptr;
  char *new;       /* pointer to working copy of string */
  size_t len;      /* length of string */

  /* length of string */
  len = strlen(str)+1;
  /* allocate memory for copy */
  new = _unur_malloc( len * sizeof(char) );
  /* copy memory */
  ptr = memcpy(new,str,len);

  /* copy characters but skip all white spaces */
  for (tmp = ptr; *tmp != '\0'; tmp++)
    if ( !isspace(*tmp) ) {
      *ptr = tolower(*tmp);
      ptr++;
    }

  /* terminate string */
  *ptr = '\0';

  /* return pointer to working copy */
  return new;

} /* end of _unur_prepare_string() */

/*---------------------------------------------------------------------------*/



/*****************************************************************************/
/**  End                                                                    **/
/*****************************************************************************/
