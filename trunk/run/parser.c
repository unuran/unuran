

#include "scanpars.h" 
#define CALC_TEST   TRUE

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
struct treenode *FuncDefinition   (char *f, int *spp, int *ecp, int *epp);
struct treenode *DefFuncDesignator(char *f, int *spp, int *ecp, int *epp, 
                                                                int *ftp);
struct treenode *DefParameterlist (char *f, int *spp, int *ecp, int *epp, 
                                                            int *paranzp);
struct treenode *Expression       (char *f, int *spp, int *ecp, int *epp);
struct treenode *SimpleExpression (char *f, int *spp, int *ecp, int *epp);
struct treenode *VTerm            (char *f, int *spp, int *ecp, int *epp);
struct treenode *Term             (char *f, int *spp, int *ecp, int *epp);
struct treenode *Factor           (char *f, int *spp, int *ecp, int *epp);
struct treenode *bas_exp          (char *f, int *spp, int *ecp, int *epp);
struct treenode *FuncDesignator   (char *f, int *spp, int *ecp, int *epp);
struct treenode *ActualParameterlist(char *f, int *spp, int *ecp, 
                                                 int *epp, int corr_panz);
struct treenode *create_node(char *symb, int token, int symbkind,
                double val, struct treenode *left, struct treenode *right);
struct treenode *tnalloc (void);
struct treenode *set_err (int err_nr, int *errcodep);
char            *readln  (char *s);
char            *erase_back_blancs (char *s); 
int             simplification(char *symb, int t, struct treenode *l, 
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
  if( quotient && eq_leaves || power && (r_0 || l_1) ){ 
     strcpy(l->symb,"1"); 
     l->token    = scs + 2; 
     l->symbkind = SCONST; 
     l->val      = 1.0; 
     l->left     = NULL; 
     l->right    = NULL; 
     return TRUE; 
  } 

  /*-- Ueberpruefen, ob 0*x,x*0,0ANDx,xAND0,0/x,0^x,0MODx => Blatt=0: --*/ 
  if( (product||and) && (l_0||r_0) || l_0 && (quotient||power||mod) ){ 
     strcpy(l->symb,"0"); 
     l->token    = scs + 1; 
     l->symbkind = SCONST; 
     l->val      = 0.0; 
     l->left     = NULL; 
     l->right    = NULL; 
     return TRUE; 
  } 

  /*- Ueberpruefen, ob x+0,x-0,x*1,x/1,x MOD 1,x^1 => l. Baum zurueck:- */ 
  if( r_0 && (plus||minus) || r_1 && (product||quotient||mod||power) ){ 
     return TRUE; 
  } 

  /* --- Ueberpruefen, ob 0+x, 1*x => rechten Baum zurueck: --- */ 
  if( l_0 && plus || l_1 && product ){ 
     strcpy(l->symb,r->symb); 
     l->token    = r->token; 
     l->symbkind = r->symbkind; 
     l->val      = r->val; 
     l->left     = r->left; 
     l->right    = r->right; 
     return TRUE; 
  } 

#if CALC_TEST

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

#endif

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
