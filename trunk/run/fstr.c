
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


/*---------------------------------------------------------------------------*/

/* structure for storing data while parsing the function string */
struct parser_data {
  char *fstr;           /* pointer to function string                        */
  int   scanpos;        /* pointer to location in string                     */
  int   lastpos;        /* position in string before scanning next token     */
  int   errno;          /* error code                                        */
  int   token;          /* scanned token                                     */
  int   length_token;   /* length of scanned token                           */
};

/* --- Prototypen: --- */
void  show_symb_tab      (void);
static void get_ds             (char *function, int  *scanpos, char *ds);
static void get_sf             (char *function, int  *scanpos, char *sf);
static int   get_uc_symbol      (char *function, int *scanpos, char *uc);
static int   nxt_symbol         (struct parser_data *pdata, char *symb);
static int   get_id_symbol      (char *function, int  *scanpos, char *id);
static int   get_ro_symbol      (char *function, int  *scanpos, char *ro);
// static int   find_index         (char *symb, int start, int end, int nxt_c);
static int find_symbol_in_table(char *symb, int start, int end);

int find_user_defined(char *symb, int nxt_c);

/* --- Bereichs-Start- und -Endemarkierungen in der Symboltabelle: --- */
static int ros, roe, aos, aoe, mos, moe, hos, hoe, oss, ose, scs, sce, sfs, sfe;

static int n_symbols;

static int uconst, ufunct, uident;

static int xxxminus, xxxplus, xxxmul, xxxdiv;
static int xxxcomma;


/*---------------------------------------------------------------------------*/
/* functions for evaluating function tree                                    */
static double v_dummy  (double l, double r);
static double v_less   (double l, double r);
static double v_equal  (double l, double r);
static double v_greater(double l, double r);
static double v_less_or(double l, double r);
static double v_unequal(double l, double r);
static double v_grtr_or(double l, double r);
static double v_or     (double l, double r);
static double v_xor    (double l, double r);
static double v_plus   (double l, double r); 
static double v_minus  (double l, double r); 
static double v_mul    (double l, double r);
static double v_and    (double l, double r);
static double v_div    (double l, double r);
static double v_mod    (double l, double r);
static double v_power  (double l, double r);
static double v_not    (double l, double r);
static double v_const  (double l, double r);
// static double v_uident (double l, double r);
// static double v_ufuncs (double l, double r);
static double v_exp    (double l, double r);
static double v_ln     (double l, double r);
static double v_log    (double l, double r);
static double v_sin    (double l, double r); 
static double v_cos    (double l, double r);
static double v_tan    (double l, double r);
static double v_sec    (double l, double r); 
static double v_sqr    (double l, double r);
static double v_abs    (double l, double r);

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
// static char *d_ufuncs();
static char *d_sqr();
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
static double tree2float (struct treenode *E_root, double x); 
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

static char variable[SYMBLENGTH];     /* name of variable */

/** IMPORTANT: strings that do not start with a digit, letter or `.' 
               must constist of a single character!!

other symbols MUST consist of single characters!!! 

there MUST not be more than two arguments for a system function !!!!

**/

static struct symbols symbol[] = {   
  {"_ROS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},            /* RelationalOperators */
  {"<",   REL_OP,1, .0,v_less,   d_const, NULL},
  {"=",   REL_OP,1, .0,v_equal,  d_const, NULL},     /** TODO: == statt = **/
  {">",   REL_OP,1, .0,v_greater,d_const, NULL},
  {"<=",  REL_OP,1, .0,v_less_or,d_const, NULL},
  {"<>",  REL_OP,1, .0,v_unequal,d_const, NULL},
  {">=",  REL_OP,1, .0,v_grtr_or,d_const, NULL},
  {"_ROE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  {"_AOS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},                /* AddingOperators */
  {"or",  ADD_OP,2, .0,v_or,     d_const, NULL},
  {"xor", ADD_OP,2, .0,v_xor,    d_const, NULL},
  {"+",   ADD_OP,2, .0,v_plus,   d_add,   NULL},
  {"-",   ADD_OP,3, .0,v_minus,  d_add,   NULL},
  {"_AOE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  {"_MOS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},           /* MultiplyingOperators */
  {"*",   MUL_OP,4, .0,v_mul,    d_mul,   NULL},
  {"and", MUL_OP,4, .0,v_and,    d_const, NULL},
  {"/",   MUL_OP,4, .0,v_div,    d_div,   NULL},
  {"mod", MUL_OP,4, .0,v_mod,    d_const, NULL},
  {"_MOE",NOSYMB,4, .0,v_dummy,  d_dummy, NULL},

  {"_HOS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},          /* HighPriorityOperators */
  {"^",   HPR_OP,5, .0,v_power,  d_power, NULL},
  {"not", HPR_OP,6, .0,v_not,    d_const, NULL},
  {"_HOE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  {"_OSS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},                   /* OtherSymbols */
  {"(",   OTHERS,0, .0,v_dummy,  d_dummy, NULL},
  {")",   OTHERS,0, .0,v_dummy,  d_dummy, NULL},
  {",",   OTHERS,0, .0,v_dummy,  d_dummy, NULL},
  {"_OSE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  {"_SCS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},                /* SystemConstants */
  {"pi",  SCONST,0, M_PI,v_const, d_const, NULL},
  {"e",   SCONST,0, M_E ,v_const, d_const, NULL},
  {"_SCE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  /* UserdefinedConstants */
  {"UCONST",UCONST, 0, .0,v_const, d_const, NULL},

  /* UserdefinedIdentifiers */
/*    {"",    UIDENT,0, .0,v_uident, d_const, NULL}, */
  {"",    UIDENT,0, .0,v_dummy, d_const, NULL},
  {"VAR", UIDENT,0, .0,v_dummy, d_const, NULL},

  /* UserdefinedFunctions */
  {"UFUNCT",UFUNCS,0, .0,v_dummy, d_dummy,NULL},

  {"_SFS",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},                /* SystemFunctions */
  {"exp", SFUNCS,1, .0,v_exp,    d_exp,   NULL},
  {"ln",  SFUNCS,1, .0,v_ln,     d_ln,    NULL},
  {"log", SFUNCS,2, .0,v_log,    d_log,   NULL},
  {"sin", SFUNCS,1, .0,v_sin,    d_sin,   NULL},
  {"cos", SFUNCS,1, .0,v_cos,    d_cos,   NULL},
  {"tan", SFUNCS,1, .0,v_tan,    d_tan,   NULL},
  {"sec", SFUNCS,1, .0,v_sec,    d_sec,   NULL},
  {"sqr", SFUNCS,1, .0,v_sqr,    d_sqr,   NULL},
  {"abs", SFUNCS,1, .0,v_abs,    d_abs,   NULL},
  {"_SFE",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},

  {"_END",NOSYMB,0, .0,v_dummy,  d_dummy, NULL},
};

/* --- Prototypen: --- */ 
static struct treenode *FuncDefinition   (struct parser_data *pdata);
static struct treenode *_unur_DefFunctDesignator(struct parser_data *pdata, int *ftp);
static struct treenode *_unur_DefParameterlist(struct parser_data *pdata, int *n_params);
static struct treenode *_unur_Expression       (struct parser_data *pdata);
static struct treenode *SimpleExpression (struct parser_data *pdata);
static struct treenode *VTerm            (struct parser_data *pdata);
static struct treenode *Term             (struct parser_data *pdata);
static struct treenode *Factor           (struct parser_data *pdata);
static struct treenode *bas_exp          (struct parser_data *pdata);
static struct treenode *_unur_FunctDesignator   (struct parser_data *pdata);
static struct treenode *_unur_ActualParameterlist(struct parser_data *pdata, int n_params);

static struct treenode *create_node(char *symb, int token,
                struct treenode *left, struct treenode *right);

static struct treenode *set_err (int err_nr,struct parser_data *pdata); 
char            *readln  (char *s);
static int             simplification(char *symb, int t, struct treenode *l, 
                                                      struct treenode *r);
static char *errorstrings[] = { 
     "OK",                                                /* 0*/ 
     "scan pointer too big",                              /* 1*/ 
     "unsigned constant area full",                       /* 2*/ 
     "identifier area full",                              /* 3*/ 
     "unkown symbol",                                     /* 4*/
     "5", "6", "7", 
     "8", "9", "10", 
     "operator expected        in string2tree",           /*11*/ 
     "expected symbol: '='     in FunctionDefinition",    /*12*/ 
     "user function expected   in DefFunctDesignator",    /*13*/ 
     "expected symbol: '('     in DefFunctDesignator",    /*14*/ 
     "expected symbol: ')'     in DefFunctDesignator",    /*15*/ 
     "user identifier expected in DefParameterlist",      /*16*/ 
     "user identifier expected in DefParameterlist",      /*17*/ 
     "expected symbol: ')'     in bas_exp",               /*18*/ 
     "unknown symbol           in bas_exp",               /*19*/ 
     "function expected        in FunctDesignator",       /*20*/ 
     "expected symbol: '('     in FunctDesignator",       /*21*/ 
     "expected symbol: ')'     in FunctDesignator",       /*22*/ 
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
	else if (strcmp(s,"_SFS") == 0) sfs = i;
	else if (strcmp(s,"_SFE") == 0) sfe = i;
  } while (!(strcmp(s,"_END") == 0));

  /* size of symbol table */
  n_symbols = i;

  /* find location of marker for user defined constant */
  uconst = find_symbol_in_table("UCONST",0,n_symbols);

  /* find location of marker for user defined function */
  ufunct = find_symbol_in_table("UFUNCT",0,n_symbols);

  uident = find_symbol_in_table("VAR",0,n_symbols);

  xxxminus = find_symbol_in_table("-",0,n_symbols);
  xxxplus = find_symbol_in_table("+",0,n_symbols);
  xxxmul = find_symbol_in_table("*",0,n_symbols);
  xxxdiv = find_symbol_in_table("/",0,n_symbols);

  xxxcomma = find_symbol_in_table(",",0,n_symbols);

  /** error check **/


  /* --- Benutzerdefinierte Eintraege loeschen: --- */

  memset(variable,'\0',SYMBLENGTH);

}

/*********************************** ************************************/



/************************************************************************/

void show_symb_tab(void)         /* Gibt die aktuelle Symboltabelle aus */
{
  int i = 0;

  printf("\nSymboltabelle:\n");
  do {
    printf("%4d: %7s inf:%3d val:%f \n",
	   i, symbol[i].name,symbol[i].info,symbol[i].val);
  } while (!(strcmp(symbol[++i].name,"_END") == 0));

}
/************************************************************************/

int nxt_symbol (struct parser_data *pdata, char *symb)

/*  int nxt_symbol (char *fstr, int *scanpos, char *symb, int *errcodep, int *errposp) */

/*  Liefert aus dem String 'fstr' das naechste Symbol ab der
 *  Position 'scanpos' und speichert es in 'symb'. 
 *  Nach der Bearbeitung zeigt *scanpos auf das Zeichen, 
 *  das dem ermittelten Symbol unmittelbar folgt. Leerzeichen
 *  werden ignoriert. In *tokenp wird der Index des Symbols in der
 *  Symboltabelle. geliefert. 
 *  *errposp zeigt auf das Zeichen, auf das *scanpos beim Aufruf zeigte.
 *  Der zurueckgelieferte Funktionswert der index des symbols in der Tabelle.
 */

{
  int token;
  int errcode = 0;
  char c;


#if 0
  if (pdata->lastpos == pdata->scanpos) {
    /* we already have read the next token */
    /* update scan position */
    pdata->scanpos = pdata->lastpos + pdata->length_token;
    /* return last scanned token */
    return pdata->token;
  }
#endif

#if 1
  int lastpos = pdata->lastpos;
  int newpos = pdata->lastpos + pdata->length_token;
  int read = (pdata->lastpos == pdata->scanpos);
#endif

  /* store position of pointer */
  pdata->lastpos = pdata->scanpos;
  
  printf("scan = %d, len = %d\n",pdata->scanpos,strlen(pdata->fstr));
  printf("%s<<\n",pdata->fstr + pdata->scanpos);

  if (pdata->scanpos > strlen(pdata->fstr))
    /* end of string */
    return 1;
  
  c = pdata->fstr[pdata->scanpos];
  
  if ( (c >= '0' && c <= '9') || c == '.') {           /* UnsignedConstant */
    errcode = get_uc_symbol(pdata->fstr,&(pdata->scanpos),symb);
    token = uconst;
  }

  else if (c >=  'a' && c <= 'z') {                       /* Identifier */
    errcode = get_id_symbol(pdata->fstr,&(pdata->scanpos),symb);

    if ( ( (token = find_symbol_in_table(symb,aos,sfe)) <= 0 ) && 
	 ( (token = find_user_defined(symb,pdata->fstr[pdata->scanpos])) <= 0 ) )
      errcode = 4;
  }

  else if ( c == '<' || c == '>' ) {              /* relation symbol */
    errcode = get_ro_symbol(pdata->fstr,&(pdata->scanpos),symb);
    if ((token = find_symbol_in_table(symb,ros,roe)) <= 0 )
      errcode = 4;
  }

  else {
    symb[0] = c; symb[1] = '\0';           /* alle anderen Zeichen */
    (pdata->scanpos)++;
    if ((token = find_symbol_in_table(symb,ros,sfe)) <= 0 )
      errcode = 4;
  }

  /* set errorcode */
  pdata->errno = errcode;

  if (read) printf("READ!\n");
  printf("last = %d, new = %d, token = %d\n", lastpos, newpos, pdata->token);
  printf("          new = %d, token = %d\n\n", pdata->scanpos, token);

  /* store token and length of token */
  pdata->token = (errcode) ? 0 : token;
  pdata->length_token = pdata->scanpos - pdata->lastpos;

  return pdata->token;
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
  if ( (sf[0]=fstr[*scanpos]) == '+' || sf[0] == '-' ) {
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

int find_symbol_in_table(char *symb, int start, int end)
     /* find symbol in table between postion start+1 and end-1 */
     /* return 0 of not found */
{
  int i;

  /* search for symbol in table */
  for (i = start + 1; i < end; i++)
    if (strcmp(symb,symbol[i].name) == 0) break;

  /* return location if symbol found and 0 otherwise */
  return ((i < end ) ? i : 0);

}

/*********************************** ************************************/

int find_user_defined(char *symb, int nxt_c)
{
  /* symbol not in table */

  if (nxt_c == '(') {
    /* --- Symbol ist Userdefined Function --- */
    /* we do not store symbol but use "UFUNCT" instead */
    return ufunct;
  }
    
  else {
    /* --- Symbol is variable  (Userdefined Identifier) --- */
    
    if (variable[0] == '\0')
      /* new variable --> store name */
      strcpy(variable,symb);   /** TODO: array-grenzen **/
    else
      /* the identifier name must match with variable name */
      if (strcmp(variable,symb) != 0) return 0;
    
    /* we use the marker for user defined identifier for variables */
    return uident;
  }
}

/*********************************** ************************************/

/************************************************************************/ 

struct treenode *
string2tree (char *functstr, int *errcodep, int *errposp) 

/*  Wandelt den String 'functstr' in einen Parse-Baum um. 'functstr' muss 
 *  in der Form f(a,b,...)=Expression vorliegen. In *errcodep und *errposp 
 *  wird im Fehlerfall eine Fehlernummer und die -position im String 
 *  'functstr' zurueckgegeben, sonst ist *errcodep Null. 
 */ 

{ 
  struct parser_data pdata;
  struct treenode *root;
  char *fstr;

  fstr = _unur_prepare_string(functstr);

  pdata.fstr = fstr;
  pdata.scanpos = 0;
  pdata.lastpos = -1;
  pdata.errno = 0;
  pdata.token = 0;
  pdata.length_token = -1;

  root = FuncDefinition(&pdata);

  if (pdata.scanpos != strlen(fstr)) {
    free (fstr);
    set_err(11,&pdata); 
    *errcodep = pdata.errno;
    *errposp  = pdata.lastpos;
    return NULL;
  }

  *errcodep = pdata.errno;
  *errposp  = pdata.lastpos;

  free (fstr);
  return root; 
} 

/*********************************** ************************************/ 

static struct treenode *
FuncDefinition (struct parser_data *pdata)

/*  FuncDefinition ::= DefFunctDesignator '=' Expression 
 * 
 *                    '=' 
 *                   /   \ 
 *  DefFunctDesignator     Expression 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             t, ft; 

  /* left hand side: DefFunctDesignator */
  l = _unur_DefFunctDesignator(pdata,&ft);

  /* we do not use this information yet --> ignore any error message */
  if (pdata->errno) return NULL;

  t = nxt_symbol(pdata,symb);

  if (strcmp(symb,"=") != 0) 
    return set_err(12,pdata); 

  r = _unur_Expression(pdata);

  if (pdata->errno) return NULL;

  root = create_node(symb,t,l,r); 

  symbol[ft].tree = root; 

  return root; 
} 

/*********************** *********************** ************************/ 

struct treenode *
_unur_DefFunctDesignator (struct parser_data *pdata, int *ftp) 

/*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')' 
 * 
 *       Identifier 
 *      /          \ 
 *  NULL            DefParameterlist 
 */ 

{ 
  struct treenode *node, *params; 
  char            fsymb[SYMBLENGTH]; 
  int             n_params; 
  int  funct;

  /* get function identifier */
  funct = nxt_symbol(pdata,fsymb);
  if (pdata->errno) return NULL;

  /* this must be a user defined function of course */
  if (symbol[funct].type != UFUNCS) return set_err(13,pdata); 

  /* read opening parenthesis '(' */
  if ( pdata->fstr[(pdata->scanpos)++] != '(' )  return set_err(14,pdata);

  /* read the parameter list */
  params = _unur_DefParameterlist(pdata,&n_params);
  if (pdata->errno) return NULL;

  /* read closing parenthesis ')' */
  if ( pdata->fstr[(pdata->scanpos)++] != ')' )  return set_err(15,pdata);

  /* store function header in node */
  node = create_node(fsymb,funct,NULL,params); 

  /* store index */
  *ftp = funct;
  
  /* return pointer to function desigatior */
  return node;
} /* end of _unur_DefFunctDesignator() */

/**************** ****************** ****************** *****************/ 

struct treenode *
_unur_DefParameterlist(struct parser_data *pdata, int *n_params) 

/*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')' 
 * 
 *       Identifier                                    ',' 
 *      /          \      oder:                       /   \ 
 *  NULL            NULL         more identifiers tree     Identifier 
 *                                                        /          \ 
 *                                                    NULL            NULL 
 */ 

{ 
  struct treenode *node, *left, *right; 
  char            symb[SYMBLENGTH]; 
  int token;

  /* read next symbol from string */
  token = nxt_symbol(pdata,symb);
 
  /* this must be a user defined identifier (i.e., a variable) of course */ 
  if (symbol[token].type != UIDENT) return set_err(16,pdata); 

  /* make node for first parameter of function and set */
  /* counter for parameters to 1                       */
  node = create_node(symb,token,NULL,NULL); 
  *n_params = 1;

  /* scan string while the character following the     */
  /* scanposition is list separator `,'                */
  while ( pdata->fstr[pdata->scanpos] == ',' ) {
    /* increment scan position */
    ++(pdata->scanpos);

    /* old node becomes left node of `,' node */
    left = node; 

    /* get next variable */
    token = nxt_symbol(pdata,symb);
    if (symbol[token].type != UIDENT) return set_err(17,pdata); 

    /* make node for next variable (becomes right node) */
    /* and update counter for parameters                */
    right = create_node(symb,token,NULL,NULL); 
    (*n_params)++; 

    /* make node for `,' separator */
    node = create_node(",",xxxcomma,left,right); 
  }
  
  /* return pointer to parameter list */
  return node; 
} /* end of _unur_DefParameterlist() */ 

/*********************** *********************** ************************/ 

static struct treenode *
set_err (int err_nr, struct parser_data *pdata)

/*  Setzt *errcodep gleich err_nr und gibt NULL zurueck. */ 

{ 
  printf("\n%d %s",err_nr,errorstrings[err_nr]); 
  if (pdata->errno) pdata->errno = err_nr;   /* nur den 1. Fehler melden */ 
  return NULL; 
} 

/*********************** *********************** ************************/ 

static struct treenode *
_unur_Expression(struct parser_data *pdata)

/* Expression ::= SimpleExpression [ RelationalOperator SimpleExpression ] 
 * 
 *                                     RelationalOperator 
 * SimpleExpression  oder:            /                  \ 
 *                    SimpleExpression                    SimpleExpression 
 */ 

{ 
  struct treenode *node, *left, *right; 
  char            symb[SYMBLENGTH]; 
  int             token, last_scanpos; 

  /* read simple expression from function string */
  left = SimpleExpression(pdata);
  if (pdata->errno) return NULL; 

  /* save position and read next symbol */ 
  last_scanpos = pdata->scanpos;
  token = nxt_symbol(pdata,symb);

  if( symbol[token].type == REL_OP ) {
    /* relation operator  --> read r.h.s.*/
    right = SimpleExpression(pdata);
    if (pdata->errno) return NULL; 
    /* create a new node */
    node = create_node(symb,token,left,right); 
  }

  else {
    /* we only have a simple expression */
    pdata->scanpos = last_scanpos; 
    node = left; 
  } 

  /* return pointer to Expression */
  return node; 
} /* end of _unur_Expression() */

/**************** ****************** ****************** *****************/ 

static struct treenode *
SimpleExpression(struct parser_data *pdata)

/*  SimpleExpression ::= VTerm { AddingOperator Term } 
 * 
 *                                AddingOperator 
 *  VTerm  oder:                 /              \ 
 *                more terms tree                Term 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             t, xscanpos; 

  root = VTerm(pdata);
  if (pdata->errno) return NULL; 
  xscanpos = pdata->scanpos;
  t = nxt_symbol(pdata,symb);

  while (symbol[t].type == ADD_OP) {
     l = root; 
     r = Term(pdata);
     if (pdata->errno) return NULL; 
     root = create_node(symb,t,l,r); 
     xscanpos = pdata->scanpos; 
     t = nxt_symbol(pdata,symb);
  }
  pdata->scanpos = xscanpos; 
  return root; 
} 

/************ ************** *************** ************** *************/ 

static struct treenode *
VTerm(struct parser_data *pdata)

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
  int             t, xscanpos = pdata->scanpos;

  t = nxt_symbol(pdata,symb);

  if( strcmp(symb,"-") == 0 ) {
    /* --- Term hat neg. Vorzeichen =>     --- */ 
    /* --- Vorzeichen in Operator wandeln: --- */ 
    l = create_node("0.0",uconst,NULL,NULL); 
    r = Term(pdata);
    if (pdata->errno) return NULL; 
    root = create_node(symb,t,l,r); 
  }
  else {
    /* --- Term hat pos. oder kein Vorzeichen: --- */ 
    if (strcmp(symb,"+") != 0) pdata->scanpos = xscanpos;  /* "+" ignorieren */ 
    root = Term(pdata);
    if (pdata->errno) return NULL; 
  } 
  return root; 
} 

/********** *********** ************ ************ *********** ***********/ 

static struct treenode *
Term(struct parser_data *pdata)

/*  Term ::= Factor [ MultiplyingOperator Factor ]
 * 
 *                                   MultiplyingOperator 
 *  Factor  oder:                   /                   \ 
 *                 more factors tree                     Factor 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             t, xscanpos; 

  root = Factor(pdata);
  if (pdata->errno) return NULL;
  xscanpos = pdata->scanpos; 

  t = nxt_symbol(pdata,symb);

  while (symbol[t].type == MUL_OP) {
     l = root; 
     r = Factor(pdata);
     if (pdata->errno) return NULL;
     root = create_node(symb,t,l,r); 
     xscanpos = pdata->scanpos; 
     t = nxt_symbol(pdata,symb);
  }
  pdata->scanpos = xscanpos; 
  return root; 
} 

/******** ********* ********** *********** ********** ********* *********/ 

static struct treenode *
Factor(struct parser_data *pdata)

/*  Factor ::= Base [ '^' Exponent ] 
 * 
 *                          '^' 
 *  bas_exp  oder:         /   \ 
 *                  bas_exp     bas_exp 
 */ 

{ 
  struct treenode *root, *l, *r; 
  char            symb[SYMBLENGTH]; 
  int             t, xscanpos; 

  l = bas_exp(pdata);
  if (pdata->errno) return NULL;
  xscanpos = pdata->scanpos; 

  t = nxt_symbol(pdata,symb);

  if( strcmp(symb,"^") == 0 ) {
    r = bas_exp(pdata);
    if (pdata->errno) return NULL;
    root = create_node(symb,t,l,r); 
  }
  else {
    pdata->scanpos = xscanpos; 
    root = l; 
  } 
  return root; 
} 

/******* ******** ******** ********* ********* ******** ******** ********/ 

static struct treenode *
bas_exp(struct parser_data *pdata)

/*  Base ::= Exponent ::= UnsignedConstant | Identifier | FunctDesignator | 
 *                        "not" Base | '(' Expression ')' 
 * 
 *       UnsignedConstant                Identifier 
 *      /                \     oder     /          \     oder 
 *  NULL                  NULL      NULL            NULL 
 * 
 *                                         "not" 
 *        FunctDesignator      oder       /     \        oder Expression 
 *                                    NULL       bas_exp 
 */ 

{ 
  struct treenode *node, *r; 
  char            symb[SYMBLENGTH]; 
  int             xscanpos = pdata->scanpos; 

  int token;


  /* get next symbol from string */
  token = nxt_symbol(pdata,symb);

  /* constants and variables */
  if( symbol[token].type==UCONST ||
      symbol[token].type==UIDENT ||
      symbol[token].type==SCONST ) {
    /* make a new end node */
    node = create_node(symb,token,NULL,NULL); 
  }

  /* "not" operator */
  else if( strcmp(symb,"not")==0 ) {
    r = bas_exp(pdata);
    if (pdata->errno) return NULL;
    node = create_node(symb,token,NULL,r); 
  }

  /* system function */
  else if( symbol[token].type == SFUNCS ) {
    pdata->scanpos = xscanpos; 
    node = _unur_FunctDesignator(pdata);
    if (pdata->errno) return NULL;
  }
  
  else if( strcmp(symb,"(") == 0 ) {
    node = _unur_Expression(pdata); 
    if (pdata->errno) return NULL;
    token = nxt_symbol(pdata,symb);
    if (strcmp(symb,")")!=0) return set_err(18,pdata); 
  }
  
  else {
/*      pdata->scanpos = scanpos;  */
    return set_err(19,pdata); 
  } 

  /* return pointer to base or exponent of an expression */
  return node; 
} 

/******* ******* ******* ******* ******* ******* ******* ******* ********/ 

struct treenode *
_unur_FunctDesignator(struct parser_data *pdata)
/*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')' 
 * 
 *       Identifier 
 *      /          \ 
 *  NULL            ActualParameterlist 
 */ 

{ 
  struct treenode *node, *params; 
  char            fsymb[SYMBLENGTH]; 
  int             n_params; 
  int  funct;

  /* get function identifier */
  funct = nxt_symbol(pdata,fsymb);
  if (pdata->errno) return NULL;

  /* this must be a system function */
  if (symbol[funct].type != SFUNCS) return set_err(20,pdata); 

  /* get number of parameter for this function */
  n_params = symbol[funct].info;

  /* read opening parenthesis '(' */
  if ( pdata->fstr[(pdata->scanpos)++] != '(' )
    return set_err(21,pdata);

  /* read the parameter list */
  params = _unur_ActualParameterlist(pdata,n_params);
  if (pdata->errno) return NULL;

  /* read closing parenthesis ')' */
  if ( pdata->fstr[(pdata->scanpos)++] != ')' )
    return set_err(15,pdata);
  
  switch (n_params) {
  case 1:    /* one argument */
    /* store function in new node */
    node = create_node(fsymb,funct,NULL,params); 
    break;
  case 2:    /* two arguments */
    /* change parameter node (which must be a `,' node) to function node */
    if (params->token != xxxcomma) return NULL;  /** ERROR **/
    node = params;
    node->token    = funct;
    node->symbkind = symbol[funct].type;
    node->symb     = symbol[funct].name;
    break;
  default:
    /** Error: should not happen **/
    return NULL;
  }

  /* return pointer to function */
  return node; 
} /* end of _unur_FunctDesignator() */

/****** ****** ****** ****** ******* ******* ****** ****** ****** *******/ 

struct treenode *
_unur_ActualParameterlist(struct parser_data *pdata, int n_params) 
/*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ] 
 * 
 *                                           ',' 
 *  Expression  oder:                       /   \ 
 *                     more expressions tree     Expression 
 */ 

{ 
  struct treenode *node, *left, *right; 
  int c_params;   /* counter for parameters */

  /* read first parameter from string ...  */
  node = _unur_Expression(pdata);
  if (pdata->errno) return NULL;
  /* .. and set counter for parameters to 1 */
  c_params = 1; 

  /* scan string while the character following the     */
  /* scanposition is list separator `,'                */
  while ( pdata->fstr[pdata->scanpos] == ',' ) {
    /* increment scan position */
    ++(pdata->scanpos);

    /* update counter for parameters */
    c_params++; 
    if (c_params > n_params) return set_err(23,pdata); 

    /* old node becomes left node of `,' node */
    left = node; 

    /* make node for next variable (becomes right node) */
    right = _unur_Expression(pdata);
    if (pdata->errno) return NULL;

    /* make node for `,' separator */
    node = create_node(",",xxxcomma,left,right); 
  }
  /* check number of parameters */
  if (c_params < n_params) return set_err(24,pdata); 

  /* return pointer to parameter list */
  return node; 
} /* end of _unur_ActualParameterlist() */

/*********************** *********************** ************************/ 

static struct treenode *
create_node (char *symb, int token,
	     struct treenode *left, struct treenode *right) 

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
    root->symb = symbol[token].name; 
    root->token    = token; 
    root->symbkind = symbol[token].type; 
    root->left     = left; 
    root->right    = right; 

    /* compute and/or store constants in val field */
    if (root->symbkind == UCONST)
      /* user defined constant, i.e. a number */
      root->val    = atof(symb);
    else if (root->symbkind == SCONST)
      /* system constant */
      root->val    = symbol[token].val;
    else
      root->val    = 0.;
  } 

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
  int comma    = t == xxxcomma;
  int product  = strcmp(   symb,"*"  ) == 0; 
  int quotient = strcmp(   symb,"/"  ) == 0; 
  int power    = strcmp(   symb,"^"  ) == 0; 
  int plus     = strcmp(   symb,"+"  ) == 0; 
  int minus    = strcmp(   symb,"-"  ) == 0; 
  int and      = strcmp(   symb,"and") == 0; 
  int mod      = strcmp(   symb,"mod") == 0; 
  int l_0      = ((l->symbkind == UCONST || l->symbkind == SCONST) && l->val == 0.);
  int r_0      = ((r->symbkind == UCONST || r->symbkind == SCONST) && r->val == 0.);
  int l_1      = ((l->symbkind == UCONST || l->symbkind == SCONST) && l->val == 1.);
  int r_1      = ((r->symbkind == UCONST || r->symbkind == SCONST) && r->val == 1.);
  //  int leaves   = l->left==NULL && l->right==NULL && r->left==NULL && r->right==NULL; 
  //  int eq_leaves= leaves && strcmp(l->symb,r->symb) == 0; 
  int l_const  = l->symbkind == SCONST || l->symbkind == UCONST; 
  int r_const  = r->symbkind == SCONST || r->symbkind == UCONST; 

  /* Do not remove node with constant 0 */
  if (l_const && l_0) l_const = FALSE;

  /* --- Transform: x/x or x^0 or 1^x  ==>  1 --- */ 
  if ( 0                 // (quotient && eq_leaves)    /** eq_leaves stimmt so nicht !!! **/
       || (power && (r_0 || l_1)) ){ 
    l->symb = NULL;
    l->token    = uconst; 
    l->symbkind = UCONST; 
    l->val      = 1.0; 
    l->left     = NULL; 
    l->right    = NULL; 
    return TRUE; 
  } 

  /*--- Transform: 0*x or x*0 or 0ANDx or xAND0 or 0/x or 0^x or 0MODx  ==>  0 ---*/ 
  if ( ((product||and) && (l_0||r_0)) || (l_0 && (quotient||power||mod)) ){ 
    l->symb = NULL;
    l->token    = uconst; 
    l->symbkind = UCONST; 
    l->val      = 0.0; 
    l->left     = NULL; 
    l->right    = NULL; 
    return TRUE; 
  } 

  /*--- Transform: x+0 or x-0 or x*1 or x/1 or xMOD1 or x^1   ==> x (nothing to do) ---*/ 
  if ( (r_0 && (plus||minus)) || (r_1 && (product||quotient||mod||power)) ) { 
    return TRUE; 
  } 

  /*--- Transform: 0+x or 1*x  ==>  x (return right node) ---*/ 
  if ( (l_0 && plus ) || (l_1 && product) ) { 
    l->symb = r->symb; 
    l->token    = r->token; 
    l->symbkind = r->symbkind; 
    l->val      = r->val; 
    l->left     = r->left; 
    l->right    = r->right; 
    return TRUE; 
  } 

  /*--- Transform: 0-x  ==>  -x (return right node) ---*/ 
  if ( l_0 && minus ) { 
    l->symb = r->symb; 
    l->token    = r->token; 
    l->symbkind = r->symbkind; 
    l->val      = - r->val; 
    l->left     = r->left; 
    l->right    = r->right; 
    return TRUE; 
  } 

  /*--- Transform: leaves are constants  -->  compute ---*/ 
  if ( l_const && r_const && !comma) { 
    l->val = (*symbol[t].vcalc)(l->val,r->val);      /* compute new value */
    l->token    = uconst;
    l->symbkind = UCONST;
    l->left     = NULL; 
    l->right    = NULL; 
    return TRUE; 
  } 

  /* no transformations */
  return FALSE; 
}

/**************** ****************** ****************** *****************/ 

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

  /* draw vertical lines in tree */
  for (i = 0, mask = 1; i < level; i++, mask <<= 1) 
    if (mask & location) 
      printf("|   "); 
    else 
      printf("    "); 

  /* print node */
  if( node != NULL ) {
    
    /* draw horizontal line in tree */
    (mask & location) ? printf("+--") : printf("\\__");

    /* print symbol */
    switch (node->symbkind) {
    case SCONST:
      printf("'%s'\t(const=%g)", node->symb,node->val);  break;
    case UCONST:
      printf("'%g'\t(const)", node->val);  break;
    case UIDENT:
      printf("'%s'\t(variable)", variable);  break;
    case UFUNCS:
      printf("'FUNCTION'");  break;
    default:
      printf("'%s'", node->symb);
    }
    /* end of line */
    printf("\n");

    /* print left and right node */
    if ( node->left || node->right) {
      /* ... unless both leaves are empty */
      show_tree(node->left,level+1,location|(mask<<1)); 
      show_tree(node->right,level+1,location); 
    }
  }

  else {  /* empty leave */
    (mask & location) ? printf("+--") : printf("\\__");
    printf("(void)\n"); 
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

double tree2float(struct treenode *E_root, double x)

/*  Erwartet in E_root den Zeiger auf den Parsebaum einer Expression und
 *  liefert als Ergebnis den numerischen Wert des Baumes.
 */

{
  double val_l, val_r;

  switch (E_root->symbkind) {
  case UCONST:
  case SCONST:
    /* node contains constant */
    return E_root->val;

  case UIDENT:
    /* variable */
    return x;

  default:
    /* use evaluation function */
    /* compute values at leaves */
    val_l = (E_root->left  != NULL) ? tree2float(E_root->left, x) : 0. ;
    val_r = (E_root->right != NULL) ? tree2float(E_root->right,x) : 0. ;

    return (*symbol[E_root->token].vcalc)(val_l,val_r);
  }
}

/*********************************** ************************************/

double v_dummy  (double l, double r) { return 0.0; }
double v_less   (double l, double r) { return (l <  r); }
double v_equal  (double l, double r) { return (l == r); }
double v_greater(double l, double r) { return (l  > r); }
double v_less_or(double l, double r) { return (l <= r); }
double v_unequal(double l, double r) { return (l != r); }
double v_grtr_or(double l, double r) { return (l >= r); }
double v_or   (double l, double r) { return   l || r;  }
double v_xor  (double l, double r) { return !(l == r); }
double v_plus (double l, double r) { return   l +  r;  }
double v_minus(double l, double r) { return   l -  r;  }
double v_mul(double l, double r) { return l *  r; }
double v_and(double l, double r) { return (double)((int)l && (int)r); }
double v_div(double l, double r) { return l /  r; }
double v_mod(double l, double r) { return (double)((int)l % (int)r); }
double v_power(double l, double r) { return pow(l,r); }
double v_not  (double l, double r) { return (double)(!(int)r)   ; }
double v_const (double l, double r) { return 0.0; }   /** nothing TODO, value in node **/
double v_exp (double l, double r) { return exp (  r); }
double v_ln  (double l, double r) { return log (  r); }
double v_log (double l, double r) { return log (r)/log(l); }
double v_sin (double l, double r) { return sin (  r); }
double v_cos (double l, double r) { return cos (  r); }
double v_tan (double l, double r) { return tan (  r); }
double v_sec (double l, double r) { return 1/cos( r); }
double v_sqr (double l, double r) { return sqrt(  r); }
double v_abs (double l, double r) { return abs (  r); }


/************************************************************************
 * Prozeduren zur analytischen Ableitung eines Parse-Baumes:            *
 ************************************************************************/

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
  struct treenode *DFD_root=root->left;    /*zeigt auf DefFunctDesignator*/
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

#if 0
char *d_ufuncs(char *par,struct treenode *w,char *l,char *r,char *dl,
	       char *dr,char *s)    /* Abl. der benutzerdef. Funktionen */
	      {
		struct treenode *parsetree;
		parsetree = part_deriv(symbol[w->token].tree,par,s);
		return strcpy(s,parsetree->left->symb);
	      }
#endif

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

{ 
/*    char            *fvonx; */
  struct treenode *root;

 
/*    fvonx= malloc((10+MAXLENGTH)*sizeof(char));      */
/*    strcpy(fvonx,"f(x)=");   */
/*    strcat(fvonx,function); */
/*    function = fvonx;  */
/*    free(fvonx);  */

  root=string2tree(function,errcodep,errposp);

  return root;
  
}

/***************************************************************************************/
double  _unur_fstr_eval_tree(struct treenode *E_root, double x)
{  
  double          result;
  struct treenode *froot;
  int             ftok;

  ftok = ufunct;  /* index for "UFUNCT" in table */

  froot=(*symbol[ftok].tree).right;             /* Achtung Fehler in Beschreibung !!! */
  result=tree2float(froot,x);
  return result;
  }
/***************************************************************************************/
char *Ntree2string(struct treenode *tree_root, char *ret_str)

{
  struct treenode  *froot;
  int              ftok;

  ftok = ufunct;  /* index for "UFUNCT" in table */
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

