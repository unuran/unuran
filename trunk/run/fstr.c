/*  FunctDefinition ::= DefFunctDesignator '=' Expression               */
     /*                    '='                                               */
     /*                   /   \                                              */
     /*  DefFunctDesignator     Expression                                   */

/*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'          */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            DefParameterlist                                    */
     /*                                                                      */

/*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'          */
     /*                                                                      */
     /*       Identifier                                 ','                 */
     /*      /          \      or:                      /   \                */
     /*  NULL            NULL      more identifiers tree     Identifier      */
     /*                                                     /          \     */
     /*                                                 NULL            NULL */

/* Expression ::= SimpleExpression [ RelationOperator SimpleExpression ] */ 
     /*                                                                      */
     /*                                    RelationOperator                  */
     /* SimpleExpression  or:             /                \                 */
     /*                   SimpleExpression                  SimpleExpression */
     /*                                                                      */

/*  SimpleExpression ::= STerm { AddingOperator Term }                  */
     /*                                                                      */
     /*                                AddingOperator                        */
     /*  STerm  or:                   /              \                       */
     /*                more terms tree                Term                   */
     /*                                                                      */

/*  STerm ::= [ '+' | '-' ] Term                                        */
     /*                                                                      */
     /*                        '-'                                           */
     /*  Term  or:            /   \                                          */
     /*                    '0'     Term                                      */
     /*                   /   \                                              */
     /*               NULL     NULL                                          */
     /*                                                                      */

/*  Term ::= Factor [ MultiplyingOperator Factor ]                      */
     /*                                                                      */
     /*                                   MultiplyingOperator                */
     /*  Factor  or:                     /                   \               */
     /*                 more factors tree                     Factor         */
     /*                                                                      */

/*  Factor ::= Base [ '^' Exponent ]                                    */
     /*                                                                      */
     /*                          '^'                                         */
     /*  Bas_Exp  or:           /   \                                        */
     /*                  Bas_Exp     Bas_Exp                                 */
     /*                                                                      */

/*  Base ::= Exponent                                                   */
     /*  Exponent ::= UnsignedConstant | Identifier | FunctDesignator |      */ 
     /*               "not" Base | '(' Expression ')'                        */ 
     /*                                                                      */
     /*       UnsignedConstant                Identifier                     */
     /*      /                \      or     /          \       or            */
     /*  NULL                  NULL      NULL            NULL                */
     /*                                                                      */
     /*                              "not"                                   */
     /*  FunctDesignator    or      /     \         or   Expression          */
     /*                         NULL       Bas_Exp                           */
     /*                                                                      */

/*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'      */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            ActualParameterlist                                 */
     /*                                                                      */

/*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]     */
     /*                                                                      */
     /*                                           ','                        */
     /*  Expression  or:                         /   \                       */
     /*                     more expressions tree     Expression             */
     /*                                                                      */

/*   UnsignedConstant ::= UnsignedInteger | UnsignedReal                */
/*   UnsignedInteger  ::= DigitSequence                                 */
/*   UnsignedReal ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor] */
/*   DigitSequence ::= Digit [ Digit [...] ]                            */
/*   Digit         ::= '0' | '1' | '2' | ... | '8' | '9'                */
/*   ScaleFactor ::= [Sign] DigitSequence                               */
/*   Sign        ::= '+' | '-'                                          */
/*   Identifier ::= Letter [ Letter | Digit [...] ]                     */
/*   Letter     ::= 'a' | 'b' | ... | 'z' | '_'                         */
/*   RelationOperator ::= RelationChar [ RelationChar ]                 */ 
/*   RelationChar     ::= '<' | '=' | '>'                               */

/** Simplification **/

  /*             Exp                    Exp
   *            /   \                  /   \ 
   *        NULL    ','       ==>     X     Y
   *               /   \ 
   *              X     Y
   */ 

  /*          Operator            Operator
   *            /   \      or      /   \        ==>     Const (result of computation)
   *       Const     Const     NULL     Const
   */ 

  /*          '+'               '*'
   *         /   \      or     /   \        ==>     X
   *        0     X           1     X
   */ 

  /*          '+'               '-'  
   *         /   \      or     /   \      or
   *        X     0           X     0
   *
   *          '*'               '/'         
   *         /   \      or     /   \      or
   *        X     1           X     1       
   *
   *          '^'               "mod"
   *         /   \      or     /     \    ==>     X
   *        X     1           X       1
   */ 

  /*          '*'               '*'         
   *         /   \      or     /   \      or
   *        0     X           X     0       
   *
   *          '/'               '^'         
   *         /   \      or     /   \      or
   *        0     X           0     X       
   *
   *         "and"             "and"         
   *         /   \      or     /   \      or
   *        0     X           X     0       
   *
   *         "mod"      
   *         /   \                        ==>     0
   *        0     X     
   */

  /*          '^'               '^'
   *         /   \      or     /   \        ==>     1
   *        X     0           1     X
   */ 

  /*              '/'              
   *          ___/   \___
   *         /           \ 
   *        X             X        ==>     1     
   *       / \           / \
   *   NULL   NULL   NULL   NULL
   */ 

/** Reorganization**/

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 

  /*            '-'                  '+'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *         '-'    Y     ==>      Y     X
   *        /   \
   *       0     X
   */ 

  /*            '*'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    0    '*'
   *              /   \                /   \
   *             0     Y              X     Y
   */ 


/*---------------------------------------------------------------------------*/

#include <ctype.h>
#include <source_unuran.h>

/*---------------------------------------------------------------------------*/

#define GENTYPE "FSTRING"      /* (pseudo) type of generator                 */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for evaluating nodes of the function tree                      **/
/*****************************************************************************/

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
static double v_exp    (double l, double r);
static double v_ln     (double l, double r);
static double v_log    (double l, double r);
static double v_sin    (double l, double r); 
static double v_cos    (double l, double r);
static double v_tan    (double l, double r);
static double v_sec    (double l, double r); 
static double v_sqr    (double l, double r);
static double v_abs    (double l, double r);
static double v_sgn    (double l, double r);

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
static char *d_sqr();
static char *d_abs();

/*****************************************************************************/
/** List of known symbols                                                   **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Symbol types for tokens                                                   */

enum {
  S_NOSYMBOL = 0,    /* token is not a symbol / is unknown                   */
  S_SFUNCT,          /* system function                                      */
  S_SCONST,          /* system constant                                      */
  S_UIDENT,          /* user defined identifier (variable)                   */
  S_UFUNCT,          /* user defined function                                */
  S_UCONST,          /* user defined constant                                */
  S_REL_OP,          /* relation operator                                    */
  S_ADD_OP,          /* addition operator                                    */
  S_MUL_OP,          /* multiplication operator                              */
  S_HPR_OP,          /* higher priority operator                             */
  S_OTHERS           /* other operator                                       */
};

/*  IMPORTANT: symbol types S_ADD_OP, S_MUL_OP and S_HPR_OP indicate         */
/*  priorities of the corresponding operators.                               */

/*---------------------------------------------------------------------------*/
/* structure for storing known symbols                                       */

#define SYMBLENGTH 20           /* Maximal length for symbols                */

/* structure for storing known symbols                                       */
struct symbols { 
  char   name[SYMBLENGTH];       /* Name des Symbols (z. B. "SIN")  */ 
  int    type;                   /* type of symbol */
  int    info;                   /* priority or 
				    number of argument for system function   */
  double val;                    /* value of constant (0. otherwise)         */

  double (*vcalc)(double l, double r);   /* pointer to function for         
					    computing value of node          */

  char            *(*dcalc)(char *par,struct ftreenode *w,
                            char *l, char *r, char *dl, char *dr,char *s);
                              /* Zeiger auf Ableitungsfunktion   */ 
};

/*---------------------------------------------------------------------------*/
/* List of known symbols                                                     */

/*  IMPORTANT: Symbols that do NOT start with a digit, a letter or `.'       */
/*  must not have more than ONE single character.                            */
/*  The only exception are relation operators.                               */
/*                                                                           */
/*  IMPORTANT: symbol types S_ADD_OP, S_MUL_OP and S_HPR_OP indicate         */
/*  priorities of the corresponding operators.                               */
/*                                                                           */
/*  In the current implementation system function must not have more than    */
/*  two arguments.                                                           */
/*                                                                           */
/*  The first entry in the list must be empty!                               */
/*                                                                           */
/*  There must be FOUR markers in the list (in this ordering!):              */
/*     _ROS   ... start of relation operators                                */
/*     _NAS   ... start of non-alphanumeric symbols (single characters!)     */
/*     _ANS   ... start of alphanumeric symbols (constants and functions)    */
/*     _END   ... end of list                                                */

static struct symbols symbol[] = {   
  {""    , S_NOSYMBOL, 0, 0.0 , v_dummy  , d_dummy }, /* void                */

  /* user defined symbols */                                             
  {"UCONST",S_UCONST , 0, 0.0 , v_const  , d_const }, /* constant            */
  {"UFUNCT",S_UFUNCT , 0, 0.0 , v_dummy  , d_dummy }, /* function            */
  {"VAR" , S_UIDENT  , 0, 0.0 , v_dummy  , d_const }, /* variable            */

  /* relation operators  */
  {"_ROS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_dummy }, /* marker for relation operators */
  {"<"   , S_REL_OP  , 1, 0.0 , v_less   , d_const },
  {"="   , S_REL_OP  , 1, 0.0 , v_equal  , d_const },
  {"=="  , S_REL_OP  , 1, 0.0 , v_equal  , d_const },
  {">"   , S_REL_OP  , 1, 0.0 , v_greater, d_const },
  {"<="  , S_REL_OP  , 1, 0.0 , v_less_or, d_const },
  {"<>"  , S_REL_OP  , 1, 0.0 , v_unequal, d_const },
  {">="  , S_REL_OP  , 1, 0.0 , v_grtr_or, d_const },

  /* special symbols */
  {"_NAS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_dummy }, /* marker for non-alphanumeric symbols */
  {"("   , S_OTHERS  , 0, 0.0 , v_dummy  , d_dummy },
  {")"   , S_OTHERS  , 0, 0.0 , v_dummy  , d_dummy },
  {","   , S_OTHERS  , 0, 0.0 , v_dummy  , d_dummy },

  /* arithmetic operators */
  {"+"   , S_ADD_OP  , 2, 0.0 , v_plus   , d_add   },
  {"-"   , S_ADD_OP  , 3, 0.0 , v_minus  , d_add   },
  {"*"   , S_MUL_OP  , 4, 0.0 , v_mul    , d_mul   },
  {"/"   , S_MUL_OP  , 4, 0.0 , v_div    , d_div   },
  {"^"   , S_HPR_OP  , 5, 0.0 , v_power  , d_power },
  {"_ANS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_dummy }, /* marker for alphanumeric symbols */
  {"mod" , S_MUL_OP  , 4, 0.0 , v_mod    , d_const },

  /* logical operators */
  {"and" , S_MUL_OP  , 2, 0.0 , v_and    , d_const },
  {"or"  , S_ADD_OP  , 4, 0.0 , v_or     , d_const },
  {"xor" , S_ADD_OP  , 4, 0.0 , v_xor    , d_const },
  {"not" , S_HPR_OP  , 6, 0.0 , v_not    , d_const },

  /* system constants */
  {"pi"  , S_SCONST  , 0, M_PI, v_const  , d_const },
  {"e"   , S_SCONST  , 0, M_E , v_const  , d_const },

  /* system functions */
  {"exp" , S_SFUNCT  , 1, 0.0 , v_exp    , d_exp   },
  {"ln"  , S_SFUNCT  , 1, 0.0 , v_ln     , d_ln    },
  {"log" , S_SFUNCT  , 2, 0.0 , v_log    , d_log   },
  {"sin" , S_SFUNCT  , 1, 0.0 , v_sin    , d_sin   },
  {"cos" , S_SFUNCT  , 1, 0.0 , v_cos    , d_cos   },
  {"tan" , S_SFUNCT  , 1, 0.0 , v_tan    , d_tan   },
  {"sec" , S_SFUNCT  , 1, 0.0 , v_sec    , d_sec   },
  {"sqr" , S_SFUNCT  , 1, 0.0 , v_sqr    , d_sqr   },
  {"abs" , S_SFUNCT  , 1, 0.0 , v_abs    , d_abs   },
  {"sgn" , S_SFUNCT  , 1, 0.0 , v_sgn    , d_const },

  /* marker for end-of-table */
  {"_END", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_dummy },
};

/*  location of special symbols in table                                     */
/** DO NOT CHANGE THEIR POSITION IN TABLE **/
static int s_uconst = 1;      /* user defined constant                       */
static int s_ufunct = 2;      /* user defined function                       */
static int s_uident = 3;      /* user defined variable                       */

/* location of special symbols */
static int s_comma, s_minus, s_plus, s_mul;

/* Marker for different regions in symbol table                              */
static int _ros_start, _ros_end;    /* relation symbols                      */
static int _nas_start, _nas_end;    /* non-alphanumeric symbols              */
static int _ans_start, _ans_end;    /* alphanumeric symbols                  */
static int _end;                    /* end of list                           */

/*****************************************************************************/
/** Other structures an definitions                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* structure for storing data while tokizing and parsing the function string */

struct parser_data {
  char  *fstr;          /* pointer to function string                        */
  int   *token;         /* array to store marker for each token in string    */
  char  *tstr;          /* working array for tokenized string                */
  char **tpos;          /* array to store pointers to each token in string   */
  int    tno;           /* pointer to token in list                          */
  int    n_tokens;      /* total number of tokens                            */
  char  *variable_name; /* name of user defined identifier                   */
  char  *function_name; /* name of user defined function                     */
  int    scanpos;       /* pointer to location in string                     */
  int    lastpos;       /* position in string before scanning next token     */
  int    len_fstr;      /* length of function string                         */
  int    errno;         /* error code                                        */
};

/*---------------------------------------------------------------------------*/
/* Error codes                                                               */

enum {
  SUCCESS = 0,          /* executed succesfully, no errors detected          */
  ERR_UNFINISHED,       /* incomplete. not all tokens parsed                 */
  ERR_UNKNOWN_SYMBOL,   /* unknown symbol in function string                 */
  ERR_EXPECT_EQUAL,     /* expected symbol: '='                              */
  ERR_EXPECT_OPEN_P,    /* expected symbol: '('                              */
  ERR_EXPECT_CLOSE_P,   /* expected symbol: ')'                              */
  ERR_INVALID_N_PARAMS, /* invalid number of parameters for function         */
  ERR_EXPECT_FUNCT,     /* function (name) expected                          */ 
  ERR_EXPECT_VAR,       /* user identifier (variable name) expected          */
  ERR_MISSING           /* more tokens expected                              */ 
};

/*****************************************************************************/
/** Prototypes                                                              **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Evaluate function tree                                                    */
/*---------------------------------------------------------------------------*/

static double _unur_fstr_eval_node (struct ftreenode *node, double x);
/*---------------------------------------------------------------------------*/
/* Evaluate function tree starting from `node' at x                          */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Initialize and destroy parser                                             */
/*---------------------------------------------------------------------------*/

static struct parser_data *_unur_fstr_parser_init (const char *fstr);
/*---------------------------------------------------------------------------*/
/* Create and initialize parser object for given function string.            */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_symbols_init (void);
/*---------------------------------------------------------------------------*/
/* Prepare table of known symbols for usage.                                 */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_parser_free (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Destroy parser object.                                                    */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Scan function string                                                      */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_tokenize (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Tokenize function string.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_next_symbol (struct parser_data *pdata, char *symb);
/*---------------------------------------------------------------------------*/
/* Get next symbol in function string.                                       */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc);
/*---------------------------------------------------------------------------*/
/* Get Unsigned Constant.                                                    */
/*  UnsignedConstant ::= UnsignedInteger | UnsignedReal                      */
/*  UnsignedInteger  ::= DigitSequence                                       */
/*  UnsignedReal ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor]   */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds);
/*---------------------------------------------------------------------------*/
/* Get Digital Sequence.                                                     */
/*  DigitSequence ::= Digit [ Digit [...] ]                                  */
/*  Digit         ::= '0' | '1' | '2' | ... | '8' | '9'                      */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf);
/*---------------------------------------------------------------------------*/
/* Get Scale Factor.                                                         */
/*  ScaleFactor ::= [Sign] DigitSequence                                     */
/*  Sign        ::= '+' | '-'                                                */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_Identifier (struct parser_data *pdata, char *id);
/*---------------------------------------------------------------------------*/
/* Get Identifier.                                                           */
/*  Identifier ::= Letter [ Letter | Digit [...] ]                           */
/*  Letter     ::= 'a' | 'b' | ... | 'z' | '_'                               */
/*  Digit      ::= '0' | '1' | ... | '9'                                     */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_RelationOperator (struct parser_data *pdata, char *ro);
/*---------------------------------------------------------------------------*/
/* Get Relation Operator.                                                    */
/*  RelationOperator ::= RelationChar [ RelationChar [...] ]                 */ 
/*  RelationChar     ::= '<' | '=' | '>'                                     */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_find_symbol (char *symb, int start, int end);
/*---------------------------------------------------------------------------*/
/* Find symbol in table between position (start+1) and (end-1).              */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, char next_char);
/*---------------------------------------------------------------------------*/
/* Find user defined symbol.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Parser.                                                                   */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_FunctDefinition (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Get user defined function.                                                */
/*  FunctDefinition ::= DefFunctDesignator '=' Expression                    */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_DefFunctDesignator (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/* Definition of user defined function.                                      */
/*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'               */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_DefParameterlist (struct parser_data *pdata, int *n_params);
/*---------------------------------------------------------------------------*/
/* Parameter list for user defined function.                                 */
/*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'                */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Expression (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Expression ::= SimpleExpression [ RelationOperator SimpleExpression ]    */ 
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_SimpleExpression (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  SimpleExpression ::= STerm { AddingOperator Term }                       */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_STerm (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  STerm ::= [ '+' | '-' ] Term                                             */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Term (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Term ::= Factor [ MultiplyingOperator Factor ]                           */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Factor (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Factor ::= Base [ '^' Exponent ]                                         */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_Bas_Exp (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  Base ::= Exponent                                                        */
/*  Exponent ::= UnsignedConstant | Identifier | FunctDesignator |           */ 
/*               "not" Base | '(' Expression ')'                             */ 
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_FunctDesignator (struct parser_data *pdata);
/*---------------------------------------------------------------------------*/
/*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'           */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_ActualParameterlist (struct parser_data *pdata, int n_params);
/*---------------------------------------------------------------------------*/
/*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]          */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_next_token (struct parser_data *pdata, int *token, char **symbol);
/*---------------------------------------------------------------------------*/
/* Get next token from list.                                                 */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_fstr_create_node (char *symb, int token,
						 struct ftreenode *left,
						 struct ftreenode *right);
/*---------------------------------------------------------------------------*/
/* Create new node.                                                          */
/*---------------------------------------------------------------------------*/

static struct ftreenode *
_unur_fstr_simplification (char *symb, int token, struct ftreenode *left, struct ftreenode *right);
/*---------------------------------------------------------------------------*/
/* Try to simpify nodes.                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_fstr_reorganize (struct ftreenode *node);
/*---------------------------------------------------------------------------*/
/* Try to reorganize tree at node.                                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Error.                                                                    */                      
/*---------------------------------------------------------------------------*/

static void _unur_fstr_error_scan (const struct parser_data *pdata, const char *symb);
/*---------------------------------------------------------------------------*/
/* Print error message when scanning function string.                        */
/*---------------------------------------------------------------------------*/

static struct ftreenode *_unur_fstr_error_parse ( struct parser_data *pdata, int errno );
/*---------------------------------------------------------------------------*/
/* Print error message when parsing function string.                         */
/*---------------------------------------------------------------------------*/

static const char *_unur_fstr_error_code ( int errno );
/*---------------------------------------------------------------------------*/
/* Print message for error number.                                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the log file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_input ( const char *fstr );
/*---------------------------------------------------------------------------*/
/* Print function string on output stream.                                   */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_token ( const struct parser_data *pdata );
/*---------------------------------------------------------------------------*/
/* Print tokenized string on output stream.                                  */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_tree( const struct parser_data *pdata,
				   const struct ftreenode *root );
/*---------------------------------------------------------------------------*/
/* Print function tree.                                                      */
/*---------------------------------------------------------------------------*/

static void _unur_fstr_debug_show_tree(const struct parser_data *pdata,
				       const struct ftreenode *node,
				       int level, int location);
/*---------------------------------------------------------------------------*/
/* Print function tree by recursion.                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** API                                                                     **/
/*****************************************************************************/

struct ftreenode *
_unur_fstr2tree(char *functstr)
     /*----------------------------------------------------------------------*/
     /* Compute funtion tree from string.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   functstr ... string containing function definition                 */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to root of function tree                                   */
     /*----------------------------------------------------------------------*/
{ 
  struct parser_data *pdata;
  struct ftreenode *root;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_input(functstr);
#endif

  /* initialize parser */
  pdata = _unur_fstr_parser_init(functstr);

  /* tokenize function string */
  _unur_fstr_tokenize(pdata);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_token(pdata);
#endif

  /* exit in case of error */
  if (pdata->errno) {
    _unur_fstr_parser_free(pdata);
    return NULL;
  }

  /* parse list of token */
  root = _unur_FunctDefinition(pdata);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_tree(pdata,root);
#endif

  /* check for possible errors */
  if (pdata->tno < pdata->n_tokens && !pdata->errno)
    _unur_fstr_error_parse(pdata,ERR_UNFINISHED); 
  if (pdata->errno) {
    _unur_fstr_parser_free(pdata);
    _unur_fstr_free(root);
    return NULL;
  }

  /* free working space */
  _unur_fstr_parser_free(pdata);

  /* return pointer to function tree */
  return root; 
} /* end of _unur_fstr2tree() */

/*---------------------------------------------------------------------------*/

double
_unur_fstr_eval_tree(struct ftreenode *root, double x)
     /*----------------------------------------------------------------------*/
     /* Evaluate function tree at x                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*   x    ... argument for which function should be evaluated           */
     /*                                                                      */
     /* return:                                                              */
     /*   result of computation                                              */
     /*----------------------------------------------------------------------*/
{  
  if (root->symbol[0] == '=')    /** TODO: gefaehrlich wegen "==" **/
    return _unur_fstr_eval_node( root->right, x );
  else
    return _unur_fstr_eval_node( root, x );
} /* end of _unur_fstr_eval_tree() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_free (struct ftreenode *root)  
     /*----------------------------------------------------------------------*/
     /* Destroy function tree.                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*----------------------------------------------------------------------*/
{ 
  if( root != NULL ) {
    _unur_fstr_free(root->right); 
    _unur_fstr_free(root->left); 
    free(root); 
  } 
} /* end of _unur_fstr_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Evaluate function                                                       **/
/*****************************************************************************/

double
_unur_fstr_eval_node (struct ftreenode *node, double x)
     /*----------------------------------------------------------------------*/
     /* Evaluate function tree starting from `node' at x                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to node in function tree                          */
     /*   x    ... argument for which function should be evaluated           */
     /*                                                                      */
     /* return:                                                              */
     /*   result of computation                                              */
     /*----------------------------------------------------------------------*/
{
  double val_l, val_r;

  switch (node->type) {
  case S_UCONST:
  case S_SCONST:
    /* node contains constant */
    return node->val;

  case S_UIDENT:
    /* variable */
    return x;

  default:
    /* use evaluation function */
    /* compute values at leaves */
    val_l = (node->left)  ? _unur_fstr_eval_node(node->left, x) : 0. ;
    val_r = (node->right) ? _unur_fstr_eval_node(node->right,x) : 0. ;

    return (*symbol[node->token].vcalc)(val_l,val_r);
  }
} /* end of _unur_fstr_eval_node() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Initialize and destroy parser                                           **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct parser_data *
_unur_fstr_parser_init ( const char *fstr )
     /*----------------------------------------------------------------------*/
     /* Create and initialize parser object for given function string.       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   fstr ... string containing function definition                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to created and initialized parser object                   */
     /*----------------------------------------------------------------------*/
{
  static int symbols_initialized = FALSE;
  struct parser_data *pdata;

  /* first we have to prepare the table of known symbols */
  if (symbols_initialized == FALSE) {
    _unur_fstr_symbols_init();
    symbols_initialized = TRUE;
  }

  /* allocate memory for parser object */
  pdata = _unur_malloc(sizeof(struct parser_data));

  /* make a working copy of the function string,                */
  /* remove all white spaces and convert to lower case letters. */
  pdata->fstr = _unur_parser_prepare_string(fstr);
  pdata->len_fstr = strlen(pdata->fstr);

  /* make arrays to store tokens in string */
  pdata->token = _unur_malloc( (pdata->len_fstr+1) * sizeof(int) );
  pdata->tpos  = _unur_malloc( pdata->len_fstr * sizeof(char *) );
  pdata->tstr  = _unur_malloc( (2*pdata->len_fstr+1) * sizeof(char) );

  /* initialize arrays */
  pdata->n_tokens = 0;
  memset(pdata->token,0,pdata->len_fstr);
  memset(pdata->tpos,0,pdata->len_fstr);
  memset(pdata->tstr,'\0',pdata->len_fstr);

  /* initialize for scanning */
  pdata->scanpos = 0;     /* scan position at beginning */
  pdata->lastpos = -1;
  pdata->errno = 0;

  /* names of user defined symbols */
  pdata->variable_name = NULL;
  pdata->function_name = NULL;

  /* initialize data for parsing */
  pdata->tno = 0;

  /* return pointer to parser object */
  return pdata;

} /* end of _unur_fstr_parser_init() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_symbols_init (void)
     /*----------------------------------------------------------------------*/
     /* Prepare table of known symbols for usage.                            */
     /*                                                                      */
     /* parameters: none                                                     */
     /*----------------------------------------------------------------------*/
{
  int i;
  char *s;

  /* none of these positions found */
  _ros_start = 0;
  _nas_start = 0;
  _ans_start = 0;

  /* find marker in list */
  for (i=0; !_end; i++) {
    s = symbol[i].name;
    if (!_ros_start) {
      if ( strcmp(s,"_ROS") == 0) _ros_start = i;
      continue;
    }
    if (!_nas_start) {
      if ( strcmp(s,"_NAS") == 0) _nas_start = i;
      continue;
    }
    if (!_ans_start) {
      if ( strcmp(s,"_ANS") == 0) _ans_start = i;
      continue;
    }
    if (strcmp(s,"_END") == 0) _end = i;
  }

  /* end of region markers for list */
  _ros_end = _nas_start;
  _nas_end = _ans_start;
  _ans_end = _end;

  /* find location of special symbols */
  s_comma = _unur_fstr_find_symbol(",",_nas_start,_nas_end);
  s_minus = _unur_fstr_find_symbol("-",_nas_start,_nas_end);
  s_plus  = _unur_fstr_find_symbol("+",_nas_start,_nas_end);
  s_mul   = _unur_fstr_find_symbol("*",_nas_start,_nas_end);

} /* end of _unur_fstr_symbols_init() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_parser_free ( struct parser_data *pdata )
     /*----------------------------------------------------------------------*/
     /* Destroy parser object.                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*----------------------------------------------------------------------*/
{
  if (pdata) {
    free(pdata->fstr);
    free(pdata->token);
    free(pdata->tpos);
    free(pdata->tstr);
    free(pdata);
  }

} /* end of _unur_fstr_parser_free() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Scan function string                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_tokenize (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Tokenize function string                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   error code of called subroutines                                   */
     /*----------------------------------------------------------------------*/
{
  int token;

  int n_token = 0;                    /* counter for token      */
  char *symb = pdata->tstr;       /* array to store token   */

  /* locate token in function string and copy into token string */
  while ((token = _unur_fstr_next_symbol(pdata,symb)) != S_NOSYMBOL) {
    pdata->token[n_token] = token;    /* marker for token */
    pdata->tpos[n_token] = symb;  /* name of token    */
    n_token++;                        /* increment counter for tokens */
    symb += pdata->scanpos - pdata->lastpos + 1;  /* skip to unused space in string */ 
  }

  /* store total number of tokens */
  pdata->n_tokens = n_token;

  /* set token pointer to first token */
  pdata->tno = 0;

  /* return error code */
  return pdata->errno;

} /* end of _unur_fstr_tokenize() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_next_symbol (struct parser_data *pdata, char *symb)
     /*----------------------------------------------------------------------*/
     /* Get next symbol in function string.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   symb  ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* return:                                                              */
     /*   error code of called subroutines                                   */
     /*----------------------------------------------------------------------*/
{
  int token;
  int errcode = 0;
  char c;

  /* store position of pointer */
  pdata->lastpos = pdata->scanpos;

  if (pdata->scanpos >= pdata->len_fstr)
    /* end of string */
    return S_NOSYMBOL;
  
  /* next character in string */
  c = pdata->fstr[pdata->scanpos];

  /* new get next symbol in list */
  if ( (c >= '0' && c <= '9') || c == '.') {
    /* Unsigned Constant */
    _unur_fstr_UnsignedConstant(pdata,symb);
    token = s_uconst;
  }

  else if (c >=  'a' && c <= 'z') {
    /* Identifier */
    _unur_fstr_Identifier(pdata,symb);

    if ( ( (token = _unur_fstr_find_symbol(symb,_ans_start,_ans_end)) == 0 ) && 
	 ( (token = _unur_fstr_find_user_defined(pdata,symb,pdata->fstr[pdata->scanpos])) <= 0 ) )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  else if ( c == '<' || c == '>' || c == '=' ) {
    /* Relation Operator */
    _unur_fstr_RelationOperator(pdata,symb);

    if ((token = _unur_fstr_find_symbol(symb,_ros_start,_ros_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  else {
    symb[0] = c; symb[1] = '\0';           /* all other charactors */
    (pdata->scanpos)++;
    if ((token = _unur_fstr_find_symbol(symb,_nas_start,_nas_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }

  /* set errorcode */
  pdata->errno = errcode;

  if (errcode) {
    _unur_fstr_error_scan (pdata,symb);
  }

  return (errcode) ? S_NOSYMBOL : token;

} /* end of _unur_fstr_next_symbol() */

/*---------------------------------------------------------------------------*/

int 
_unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc)
     /*----------------------------------------------------------------------*/
     /* Get Unsigned Constant                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   uc    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   UnsignedConstant ::= UnsignedInteger | UnsignedReal                */
     /*   UnsignedInteger  ::= DigitSequence                                 */
     /*   UnsignedReal ::= UnsignedInteger ['.' DigitSequence] ['e' ScaleFactor] */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  /* store scan position */
  int startpos = pdata->scanpos;

  /* copy digit sequence into uc */
  _unur_fstr_DigitalSequence(pdata,uc);

  if( pdata->fstr[pdata->scanpos] == '.' ) {
    /* decimal point  --> copy point and following digits */
    *(uc + pdata->scanpos - startpos) = '.';
    (pdata->scanpos)++;
    _unur_fstr_DigitalSequence(pdata, uc + pdata->scanpos - startpos);
  }

  if( pdata->fstr[pdata->scanpos] == 'e' ) {
    /* exponent --> copy indicator 'E' and following [sign and] digits */
    *(uc + pdata->scanpos - startpos) = 'e';
    (pdata->scanpos)++;
    _unur_fstr_ScaleFactor(pdata, uc + pdata->scanpos - startpos);
  }

  return 0;
} /* end of _unur_fstr_UnsignedConstant() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds)
     /*----------------------------------------------------------------------*/
     /* Get Digital Sequence                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   ds    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   DigitSequence ::= Digit [ Digit [...] ]                            */
     /*   Digit         ::= '0' | '1' | '2' | ... | '8' | '9'                */
     /*----------------------------------------------------------------------*/
{
  /* copy digit */
  while ( (*ds = pdata->fstr[pdata->scanpos]) >= '0' && *ds <= '9' ) {
     ds++;
     (pdata->scanpos)++;
  }
  /* terminate string */
  *ds = '\0';

  return 0;
} /* end of _unur_fstr_DigitalSequence() */

/*---------------------------------------------------------------------------*/

int 
_unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf)
     /*----------------------------------------------------------------------*/
     /* Get Scale Factor                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   sf    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   ScaleFactor ::= [Sign] DigitSequence                               */
     /*   Sign        ::= '+' | '-'                                          */
     /*----------------------------------------------------------------------*/
{
  /* copy sign */
  if ( (sf[0] = pdata->fstr[pdata->scanpos]) == '+' || sf[0] == '-' ) {
     sf++;
     (pdata->scanpos)++;
  }
  /* copy digital sequence (after sign) */
  _unur_fstr_DigitalSequence(pdata,sf);

  return 0;
} /* _unur_fstr_ScaleFactor() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_Identifier (struct parser_data *pdata, char *id)
     /*----------------------------------------------------------------------*/
     /* Get Identifier                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   id    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   Identifier ::= Letter [ Letter | Digit [...] ]                     */
     /*   Letter     ::= 'a' | 'b' | ... | 'z' | '_'                         */
     /*   Digit      ::= '0' | '1' | ... | '9'                               */
     /*----------------------------------------------------------------------*/
{
  /* copy word */
  while ( ((*id = pdata->fstr[pdata->scanpos]) >= 'a' && *id <= 'z')
	  || *id == '_' 
	  || ( *id >= '0' && *id <= '9')) {
    id++;
    (pdata->scanpos)++;
  }
  /* terminate string */
  *id = '\0';

  return 0;
} /* end of _unur_fstr_Identifier() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_RelationOperator (struct parser_data *pdata, char *ro)
     /*----------------------------------------------------------------------*/
     /* Get Relation Operator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   ro    ... pointer to array for storing symbol                      */
     /*                                                                      */
     /* Syntax:                                                              */
     /*   RelationOperator ::= RelationChar [ RelationChar ]                 */ 
     /*   RelationChar     ::= '<' | '=' | '>'                               */
     /*----------------------------------------------------------------------*/
{
  /* copy relation operator */
  while ((*ro = pdata->fstr[pdata->scanpos]) == '<' || *ro == '=' || *ro == '>') {
    ro++;
    (pdata->scanpos)++;
  }
  /* terminate string */
  *ro = '\0';

  return 0;
} /* _unur_fstr_RelationOperator() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_find_symbol (char *symb, int start, int end)
     /*----------------------------------------------------------------------*/
     /* find symbol in table between position (start+1) and (end-1)          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb  ... symbol to look for                                       */
     /*   start ... starting point for searching                             */
     /*   end   ... endpoint for searching                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   location in table                                                  */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* search for symbol in table */
  for (i = start + 1; i < end; i++)
    if (strcmp(symb,symbol[i].name) == 0) break;

  /* return location if symbol found and 0 otherwise */
  return ((i < end ) ? i : 0);

} /* end of _unur_fstr_find_symbol() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, char next_char)
     /*----------------------------------------------------------------------*/
     /* Find user defined symbol.                                            */
     /* If there are no user defined symbols yet, store it.                  */
     /* If there are already user defined symbols, and this one is new,      */
     /* make return an errorcode.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   symb      ... symbol to look for                                   */
     /*   next_char ... character following given symbol                     */
     /*                                                                      */
     /* return:                                                              */
     /*   location in table                                                  */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* we use next_char to distinguish between variables and functions */

  if (next_char == '(') {
    /* symbol is user defined function */
    if (pdata->function_name == NULL) {
      /* new function --> store name */
      pdata->function_name = symb;
      /* return marker for user defined identifier */
      return s_ufunct;
    }
    else
      /* the identifier name must match with variable name */
      return (strcmp(pdata->function_name,symb) == 0) ? s_ufunct : 0;
  }
 
  else {
    /* symbol is user defined identifier (i.e. a variable) */
    if (pdata->variable_name == NULL) {
      /* new variable --> store name */
      pdata->variable_name = symb;
      /* return marker for user defined identifier */
      return s_uident;
    }
    else
      /* the identifier name must match with variable name */
      return (strcmp(pdata->variable_name,symb) == 0) ? s_uident : 0;
  }

} /* end of _unur_fstr_find_user_defined() */

/*****************************************************************************/
/** Parser                                                                  **/
/*****************************************************************************/

struct ftreenode *
_unur_FunctDefinition (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Get user defined function.                                           */
     /*                                                                      */
     /*  FunctDefinition ::= DefFunctDesignator '=' Expression               */
     /*                                                                      */
     /*                    '='                                               */
     /*                   /   \                                              */
     /*  DefFunctDesignator     Expression                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token; 

  /* left hand side: DefFunctDesignator */
  left = _unur_DefFunctDesignator(pdata);
  if (pdata->errno) return NULL;

  /* next token must be "=" sign */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       strcmp(symb,"=") != 0 )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_EQUAL); 

  /* right hand side: function term */
  right = _unur_Expression(pdata);
  if (pdata->errno) return NULL;

  /* store function in node */
  node = _unur_fstr_create_node(symb,token,left,right); 

  /* return pointer to function */
  return node; 
} /* end of _unur_FunctDefinition() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_DefFunctDesignator (struct parser_data *pdata) 
     /*----------------------------------------------------------------------*/
     /* Definition of user defined function.                                 */
     /*                                                                      */
     /*  DefFunctDesignator ::= Identifier '(' DefParameterlist ')'          */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            DefParameterlist                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              n_params; 
  int              funct, token;

  /* get function identifier */
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) == FALSE ||
       symbol[funct].type != S_UFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT); 

  /* read opening parenthesis '(' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P); 

  /* read the parameter list */
  params = _unur_DefParameterlist(pdata,&n_params);
  if (pdata->errno) return NULL;

  /* read closing parenthesis ')' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       symb[0] != ')' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P); 

  /* store function header in node */
  node = _unur_fstr_create_node(fsymb,funct,NULL,params); 

  /* return pointer to function desigatior */
  return node;
} /* end of _unur_DefFunctDesignator() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_DefParameterlist(struct parser_data *pdata, int *n_params) 
     /*----------------------------------------------------------------------*/
     /* Parameter list for user defined function.                            */
     /*                                                                      */
     /*  DefParameterlist ::= '(' Identifier [ ',' Identifier ] ')'          */
     /*                                                                      */
     /*       Identifier                                 ','                 */
     /*      /          \      or:                      /   \                */
     /*  NULL            NULL      more identifiers tree     Identifier      */
     /*                                                     /          \     */
     /*                                                 NULL            NULL */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   n_params  ... detected number of parameters for defined function   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* read user defined identifier, i.e. a variable */ 
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       symbol[token].type != S_UIDENT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR);

  /* make node for first parameter of function and set */
  /* counter for parameters to 1                       */
  node = _unur_fstr_create_node(symb,token,NULL,NULL); 
  *n_params = 1;

  /* scan token list while we find a list separator `,' */
  while ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
	  symb[0] == ',' ) {

    /* old node becomes left node of `,' node */
    left = node; 

    /* get next variable */
    if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
	 symbol[token].type != S_UIDENT )
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR);

    /* make node for next variable (becomes right node) */
    /* and update counter for parameters                */
    right = _unur_fstr_create_node(symb,token,NULL,NULL); 
    (*n_params)++; 

    /* make node for `,' separator */
    node = _unur_fstr_create_node(",",s_comma,left,right); 
  }

  /* set token pointer to first element after parameter list */
  --(pdata->tno);

  /* return pointer to parameter list */
  return node; 
} /* end of _unur_DefParameterlist() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Expression (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /* Expression ::= SimpleExpression [ RelationOperator SimpleExpression ] */ 
     /*                                                                      */
     /*                                    RelationOperator                  */
     /* SimpleExpression  or:             /                \                 */
     /*                   SimpleExpression                  SimpleExpression */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* read simple expression from function string */
  left = _unur_SimpleExpression(pdata);
  if (pdata->errno) return NULL; 

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) &&
       symbol[token].type == S_REL_OP ) {
    /* relation operator  --> read r.h.s.*/
    right = _unur_SimpleExpression(pdata);
    if (pdata->errno) return NULL; 
    /* create a new node */
    node = _unur_fstr_create_node(symb,token,left,right); 
  }

  else {
    /* we only have a simple expression */
    /* hence we set token pointer to get same token again */
    --(pdata->tno);
    node = left; 
  } 

  /* return pointer to Expression */
  return node; 
} /* end of _unur_Expression() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_SimpleExpression (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  SimpleExpression ::= STerm { AddingOperator Term }                  */
     /*                                                                      */
     /*                                AddingOperator                        */
     /*  STerm  or:                   /              \                       */
     /*                more terms tree                Term                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* get next Term in string */
  node = _unur_STerm(pdata);
  if (pdata->errno) return NULL;

  /* get next token */
  while ( _unur_fstr_next_token(pdata,&token,&symb) &&
	  symbol[token].type == S_ADD_OP) {
    /* get Term after adding operator and        */
    /* use node a right node of new operator node */
    left = node; 

    right = _unur_Term(pdata);
    if (pdata->errno) return NULL; 

    node = _unur_fstr_create_node(symb,token,left,right); 
  }
  
  /* set token pointer to get same token again */
  --(pdata->tno);

  /* return pointer */
  return node; 
} /* end of _unur_SimpleExpression() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_STerm (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  STerm ::= [ '+' | '-' ] Term                                        */
     /*                                                                      */
     /*                        '-'                                           */
     /*  Term  or:            /   \                                          */
     /*                    '0'     Term                                      */
     /*                   /   \                                              */
     /*               NULL     NULL                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb; 
  int              token;

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) &&
       symb[0] == '-' ) {
    /* term with negative sign                 */
    /* "-" is interpreted as binary operator   */
    /* thus "0" is added in front of it         */
    left = _unur_fstr_create_node("0",s_uconst,NULL,NULL); 
    right = _unur_Term(pdata);
    if (pdata->errno) return NULL; 

    node = _unur_fstr_create_node(symb,token,left,right); 
  }

  else {
    /* term with positive sign or without any sign */
    if( symb[0] != '+' ) {
      /* set pointer to previous token again */
      --(pdata->tno);
    }
    node = _unur_Term(pdata);
    if (pdata->errno) return NULL; 
  } 

  /* return pointer to term */
  return node; 
} /* end of _unur_STerm() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Term (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Term ::= Factor [ MultiplyingOperator Factor ]                      */
     /*                                                                      */
     /*                                   MultiplyingOperator                */
     /*  Factor  or:                     /                   \               */
     /*                 more factors tree                     Factor         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* get next factor of multiplication */
  node = _unur_Factor(pdata);
  if (pdata->errno) return NULL;

  /* get next token */
  while ( _unur_fstr_next_token(pdata,&token,&symb) &&
	  symbol[token].type == S_MUL_OP ) {
    /* get Factor after multiplication operator and */
    /* use node a right node of new operator node   */
     left = node; 

     right = _unur_Factor(pdata);
     if (pdata->errno) return NULL;

     node = _unur_fstr_create_node(symb,token,left,right); 
  }

  /* set token pointer to get same token again */
  --(pdata->tno);

  /* return pointer to term */
  return node; 
} /* end of _unur_Term() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Factor (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Factor ::= Base [ '^' Exponent ]                                    */
     /*                                                                      */
     /*                          '^'                                         */
     /*  Bas_Exp  or:           /   \                                        */
     /*                  Bas_Exp     Bas_Exp                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;

  /* get base of factor */
  left = _unur_Bas_Exp(pdata);
  if (pdata->errno) return NULL;

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) &&
       symb[0] == '^' ) {
    /* get exponent of factor */
    right = _unur_Bas_Exp(pdata);
    if (pdata->errno) return NULL;

    /* and create node for '^' operator */
    node = _unur_fstr_create_node(symb,token,left,right); 
  }

  else {
    /* no exponent                         */
    /* set pointer to previous token again */
    --(pdata->tno);
    /* and return base */
    node = left; 
  } 

  /* return pointer to factor */
  return node; 
} /* end of _unur_Factor() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_Bas_Exp (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  Base ::= Exponent                                                   */
     /*  Exponent ::= UnsignedConstant | Identifier | FunctDesignator |      */ 
     /*               "not" Base | '(' Expression ')'                        */ 
     /*                                                                      */
     /*       UnsignedConstant                Identifier                     */
     /*      /                \      or     /          \       or            */
     /*  NULL                  NULL      NULL            NULL                */
     /*                                                                      */
     /*                              "not"                                   */
     /*  FunctDesignator    or      /     \         or   Expression          */
     /*                         NULL       Bas_Exp                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *right; 
  char             *symb;
  int              token;

  /* get next token */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE)
    return _unur_fstr_error_parse(pdata,7); 

  /* constant or and variable */
  if( symbol[token].type==S_UCONST ||
      symbol[token].type==S_UIDENT ||
      symbol[token].type==S_SCONST ) {
    /* make a new end node */
    node = _unur_fstr_create_node(symb,token,NULL,NULL); 
  }

  /* "not" operator */
  else if( strcmp(symb,"not")==0 ) {
    right = _unur_Bas_Exp(pdata);
    if (pdata->errno) return NULL;
    node = _unur_fstr_create_node(symb,token,NULL,right); 
  }

  /* system function */
  else if( symbol[token].type == S_SFUNCT ) {
    /* set pointer to previous token again */
    --(pdata->tno);
    /* and get function */
    node = _unur_FunctDesignator(pdata);
    if (pdata->errno) return NULL;
  }
  
  else if( symb[0] == '(' ) {
    /* opening parenthesis  --> read expression in side parenthesis */
    node = _unur_Expression(pdata); 
    if (pdata->errno) return NULL;

    /* next symbol must be closing parenthesis */
    if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
	 symb[0] != ')' )
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P);
  }
  
  else {
    /* unkown symbol */
    --(pdata->tno);
    return _unur_fstr_error_parse(pdata,ERR_UNKNOWN_SYMBOL);
  } 

  /* return pointer to base or exponent of an expression */
  return node; 
} /* end of _unur_Bas_Exp() */ 

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_FunctDesignator (struct parser_data *pdata)
     /*----------------------------------------------------------------------*/
     /*  FunctDesignator ::= FuncIdentifier '(' ActualParameterlist ')'      */
     /*                                                                      */
     /*       Identifier                                                     */
     /*      /          \                                                    */
     /*  NULL            ActualParameterlist                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              funct, token;
  int              n_params; 

  /* get function identifier for system function */
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) == FALSE ||
       symbol[funct].type != S_SFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT);

  /* get number of parameter for this function */
  n_params = symbol[funct].info;

  /* read opening parenthesis '(' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P);

  /* read the parameter list */
  params = _unur_ActualParameterlist(pdata,n_params);
  if (pdata->errno) return NULL;

  /* read closing parenthesis ')' */
  if ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
       symb[0] != ')' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P);
  
  /* store function in new node */
  node = _unur_fstr_create_node(fsymb,funct,NULL,params); 

  /* return pointer to function */
  return node; 
} /* end of _unur_FunctDesignator() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_ActualParameterlist (struct parser_data *pdata, int n_params) 
     /*----------------------------------------------------------------------*/
     /*  ActualParameterlist ::= ActualParameter [ ',' ActualParameter ]     */
     /*                                                                      */
     /*                                           ','                        */
     /*  Expression  or:                         /   \                       */
     /*                     more expressions tree     Expression             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata     ... pointer to parser object                             */
     /*   n_params  ... number of expected parameters in list                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to function tree                                           */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  int              c_params;   /* counter for parameters */

  /* read first parameter from string ...  */
  node = _unur_Expression(pdata);
  if (pdata->errno) return NULL;

  /* .. and set counter for parameters to 1 */
  c_params = 1; 

  /* scan token list while we find a list separator `,' */
  while ( _unur_fstr_next_token(pdata,&token,&symb) == FALSE ||
	  symb[0] == ',' ) {

    /* update counter for parameters */
    c_params++; 
    if (c_params > n_params)
      return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS);

    /* old node becomes left node of `,' node */
    left = node; 

    /* make node for next variable (becomes right node) */
    right = _unur_Expression(pdata);
    if (pdata->errno) return NULL;
    
    /* make node for `,' separator */
    node = _unur_fstr_create_node(",",s_comma,left,right); 
  }

  /* set token pointer to first element after parameter list */
  --(pdata->tno);

  /* check number of parameters */
  if (c_params < n_params)
    return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS);

  /* return pointer to parameter list */
  return node; 
} /* end of _unur_ActualParameterlist() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_next_token (struct parser_data *pdata, int *token, char **symbol)
     /*----------------------------------------------------------------------*/
     /* Get next token from list.                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata  ... pointer to parser object                                */
     /*   token  ... to store token                                          */
     /*   symbol ... to store symbol for token                               */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... if successful                                                */
     /*   0 ... if there are no tokens any more in the list                  */
     /*----------------------------------------------------------------------*/
{
  if (pdata->tno < pdata->n_tokens) {
    /* return token and increment scan position */
    *token = pdata->token[pdata->tno];
    *symbol = pdata->tpos[pdata->tno];
    ++(pdata->tno);
    return TRUE;
  }
  else {
    /* no more tokens */
    ++(pdata->tno);
    return FALSE;
  }

} /* end of _unur_fstr_next_token() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_create_node (char *symb, int token, struct ftreenode *left, struct ftreenode *right) 
     /*----------------------------------------------------------------------*/
     /* Create new node.                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb   ... symbol string                                           */
     /*   token  ... token for which node should be created                  */
     /*   left   ... pointer to left node                                    */
     /*   right  ... pointer to right node                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new node                                                */
     /*----------------------------------------------------------------------*/
{ 
  struct ftreenode *node; 

  if ( (node = _unur_fstr_simplification(symb,token,left,right)) ) {
    /* node has been simplified  -->  use left node, remove right node */
  }
  else {

    /* make new node */
    node = _unur_malloc(sizeof(struct ftreenode)); 
    node->symbol = symbol[token].name; 
    node->token  = token; 
    node->type   = symbol[token].type; 
    node->left   = left; 
    node->right  = right; 

    /* compute and/or store constants in val field */
    switch (symbol[token].type) {
    case S_UCONST:      /* user defined constant, i.e. a number */
      node->val = atof(symb);  break;
    case S_SCONST:      /* system constant */
      node->val = symbol[token].val;  break;
    default:
      node->val = 0.;
    }
  } 

  /* it might be possible to reorganize the tree */
  _unur_fstr_reorganize(node); 
  
  /* return node */
  return node; 
  
} /* end of _unur_fstr_create_node() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_simplification (char *symb, int token,
			   struct ftreenode *left, struct ftreenode *right) 
     /*----------------------------------------------------------------------*/
     /* Try to simpify nodes                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   symb   ... symbol string                                           */
     /*   token  ... token for which new node should be created              */
     /*   left   ... pointer to left node                                    */
     /*   right  ... pointer to right node                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to simplified node                                         */
     /*   NULL if no simplication is possible                                */
     /*----------------------------------------------------------------------*/
{
  /* some booleans */
  int l_const = left  && (left->type  == S_SCONST || left->type  == S_UCONST); 
  int r_const = right && (right->type == S_SCONST || right->type == S_UCONST); 
  int l_0     = (l_const && left->val == 0.);
  int l_1     = (l_const && left->val == 1.);
  int r_0     = (r_const && right->val == 0.);
  int r_1     = (r_const && right->val == 1.);
  int and, mod;

  char s = symb[0];

  /*             Exp                    Exp
   *            /   \                  /   \ 
   *        NULL    ','       ==>     X     Y
   *               /   \ 
   *              X     Y
   */ 
  if ( left == NULL && right && right->symbol[0] == ',' ) {
    right->token  = token;
    right->symbol = symbol[token].name; 
    right->type   = symbol[token].type;
    return right; 
   }

  /*          Operator            Operator
   *            /   \      or      /   \        ==>     Const (result of computation)
   *       Const     Const     NULL     Const
   */ 
  if ( (l_const || left==NULL) && r_const && s!=',') { 
    /* compute new value */
    right->val   = ( (left) 
		     ? (*symbol[token].vcalc)(left->val,right->val)
		     : (*symbol[token].vcalc)(0.,right->val) );
    right->token = s_uconst;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    if (left) free (left);
    return right; 
  } 

  /*          '+'               '*'
   *         /   \      or     /   \        ==>     X
   *        0     X           1     X
   */ 
  if ( (l_0 && s=='+' ) || (l_1 && s=='*') ) { 
    free (left);
    return right;
  } 

  /*          '+'               '-'  
   *         /   \      or     /   \      or
   *        X     0           X     0
   *
   *          '*'               '/'         
   *         /   \      or     /   \      or
   *        X     1           X     1       
   *
   *          '^'               "mod"
   *         /   \      or     /     \    ==>     X
   *        X     1           X       1
   */ 
  mod = (strcmp(symb,"mod")==0);
  if ( (r_0 && (s=='+' || s=='-')) ||
       (r_1 && (s=='*' || s=='/' || s=='^' || mod)) ) {
    free (right);
    return left;
  }

  /*          '*'               '*'         
   *         /   \      or     /   \      or
   *        0     X           X     0       
   *
   *          '/'               '^'         
   *         /   \      or     /   \      or
   *        0     X           0     X       
   *
   *         "and"             "and"         
   *         /   \      or     /   \      or
   *        0     X           X     0       
   *
   *         "mod"      
   *         /   \                        ==>     0
   *        0     X     
   */
  and = (strcmp(symb,"and")==0);
  if ( l_0 && (s=='*' || s=='/' || s=='^' || and || mod) ) {
    free (right);
    return left;
  }
  if (r_0 && (s=='*' || and ) ) {
    free (left);
    return right;
  }

  /*          '^'               '^'
   *         /   \      or     /   \        ==>     1
   *        X     0           1     X
   */ 
  if (r_0 && s=='^') {
    free (left);
    right->val = 1.;
    return right;
  }
  if (l_1 && s=='^') {
    free (right);
    return left;
  }

  /*              '/'              
   *          ___/   \___
   *         /           \ 
   *        X             X        ==>     1     
   *       / \           / \
   *   NULL   NULL   NULL   NULL
   */ 
  if ( ( left  && left->left==NULL  && left->right==NULL  && 
	 right && right->left==NULL && right->right==NULL &&
	 strcmp(left->symbol,right->symbol)== 0 ) ) {
    free (left);
    right->token = s_uconst;
    right->symbol= symbol[s_uconst].name; 
    right->val   = 1.;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    return right; 
  }

  /* no simplification */
  return NULL; 
} /* _unur_fstr_simplification() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_reorganize (struct ftreenode *node) 
     /*----------------------------------------------------------------------*/
     /* Try to reorganize tree at node                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   node   ... pointer to node in tree                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   1 if successful                                                    */
     /*   0 otherwise                                                        */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *tmp;
  char symb = node->symbol[0];

  /* some booleans */
  int rl_0 = (right && right->left && right->left->type == S_UCONST && right->left->val == 0.);
  int ll_0 = (left  && left->left  && left->left->type  == S_UCONST && left->left->val  == 0.);

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 
  if ( rl_0 && symb=='+' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    node->right  = right->right; 
    free (right->left);
    free (right);
    return 1;
  }

  /*            '-'                  '+'  
   *           /   \                /   \ 
   *          X    '-'      ==>    X     Y
   *              /   \
   *             0     Y
   */ 
  if ( rl_0 && symb=='-' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_plus].name; 
    node->token  = s_plus; 
    node->type   = symbol[s_plus].type; 
    node->right  = right->right; 
    free (right->left);
    free (right);
    return 1;
  }

  /*            '+'                  '-'  
   *           /   \                /   \ 
   *         '-'    Y     ==>      Y     X
   *        /   \
   *       0     X
   */ 
  if ( ll_0 && symb=='+' && left->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    tmp = node->right;
    node->right  = left->right; 
    node->left   = tmp;
    free (left->left);
    free (left);
    _unur_fstr_debug_tree(NULL,node);
    return 1;
  }

  /*            '*'                  '-'  
   *           /   \                /   \ 
   *          X    '-'      ==>    0    '*'
   *              /   \                /   \
   *             0     Y              X     Y
   */ 
  if ( rl_0 && symb=='*' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    right->symbol= symbol[s_mul].name; 
    right->token = s_mul; 
    right->type  = symbol[s_mul].type; 
    tmp = left;
    node->left   = node->right->left;
    node->right  = tmp;
    return 1;
  }

  /* no reorganization */
  return 0;

} /* end of _unur_fstr_reorganize() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for evaluating nodes of the function tree                      **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

double v_dummy  (double l, double r) { return 0.; }
double v_const  (double l, double r) { return 0.; }  /* nothing to do, value in node */

double v_less   (double l, double r) { return (double)(l <  r); }
double v_equal  (double l, double r) { return (double)(l == r); }
double v_greater(double l, double r) { return (double)(l >  r); }
double v_less_or(double l, double r) { return (double)(l <= r); }
double v_unequal(double l, double r) { return (double)(l != r); }
double v_grtr_or(double l, double r) { return (double)(l >= r); }

double v_or     (double l, double r) { return (double)((int)l || (int)r);  }
double v_xor    (double l, double r) { return (double)((l==0 && r!=0) || (l!=0 && r== 0)); }
double v_and    (double l, double r) { return (double)((int)l && (int)r); }
double v_not    (double l, double r) { return (double)(!((int)r))   ; }

double v_plus   (double l, double r) { return (l + r); }
double v_minus  (double l, double r) { return (l - r); }
double v_mul    (double l, double r) { return (l * r); }
double v_div    (double l, double r) { return (l / r); }
double v_power  (double l, double r) { return pow(l,r); }
double v_mod    (double l, double r) { return (double)((int)l % (int)r); }

double v_exp    (double l, double r) { return exp(r); }
double v_ln     (double l, double r) { return log(r); }
double v_log    (double l, double r) { return log(r)/log(l); }
double v_sin    (double l, double r) { return sin(r); }
double v_cos    (double l, double r) { return cos(r); }
double v_tan    (double l, double r) { return tan(r); }
double v_sec    (double l, double r) { return 1./cos(r); }
double v_sqr    (double l, double r) { return sqrt(r); }
double v_abs    (double l, double r) { return abs(r); }
double v_sgn    (double l, double r) { return ((r<0.) ? -1. : ((r>0.) ? 1. : 0.)); }

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Error                                                                   **/
/*****************************************************************************/

void
_unur_fstr_error_scan (const struct parser_data *pdata, const char *symb)
     /*----------------------------------------------------------------------*/
     /* Print error message when scanning function string                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   symb  ... pointer to unknown symbol                                */  
     /*----------------------------------------------------------------------*/
{
  char format[124];
  int wsp;
  
  /* set unuran error code */
  unur_errno = UNUR_ERR_STR_UNKNOWN;

  /* print unknown symbol */
  _unur_stream_printf_simple ( "%s: error: unknown symbol `%s'\n",GENTYPE,symb);

  /* print scanned part of function string to error stream */
  sprintf( format, "%s: %%.%ds\n", GENTYPE,pdata->lastpos);
  _unur_stream_printf_simple ( format,pdata->fstr);

  /* print remaining part of function string including unknown symbol */
  wsp = pdata->lastpos-3;
  wsp = max(wsp,1);
  sprintf( format, "%s: --> %%%d.%ds",GENTYPE,wsp,wsp);
  _unur_stream_printf_simple ( format, "");
  _unur_stream_printf_simple ( "%s\n",pdata->fstr + pdata->lastpos);

  _unur_stream_printf_simple ( "%s:\n",GENTYPE );

} /* end of _unur_fstr_error_scan() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_error_parse ( struct parser_data *pdata, int errno )
     /*----------------------------------------------------------------------*/
     /* Print error message when parsing function string                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   errno ... error number                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   NULL                                                               */
     /*----------------------------------------------------------------------*/
{ 
  int i;

  _unur_stream_printf_simple ( "%s: error: %s\n",GENTYPE,_unur_fstr_error_code(errno));

  _unur_stream_printf_simple ( "%s: ",GENTYPE );
  for (i=0; i<pdata->tno-1; i++)
    _unur_stream_printf_simple ( "%s ",pdata->tpos[i]);

  _unur_stream_printf_simple ( "\n%s:  -->   ",GENTYPE );
  for (i=pdata->tno-1; i<pdata->n_tokens; i++)
    _unur_stream_printf_simple ( "%s ",pdata->tpos[i]);

  _unur_stream_printf_simple ( "\n%s:\n",GENTYPE );

  /* set unuran error code */
  unur_errno = UNUR_ERR_STR_UNKNOWN;
  
  /* set parser error */
  if (!pdata->errno) pdata->errno = errno;

  return NULL; 

} /* _unur_fstr_error_parse() */ 

/*---------------------------------------------------------------------------*/

const char *
_unur_fstr_error_code ( int errno )
     /*----------------------------------------------------------------------*/
     /* Print message for error number                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   errno ... error number                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to message string                                          */
     /*----------------------------------------------------------------------*/
{
  switch (errno) {
  case ERR_UNFINISHED:
    return "incomplete. not all tokens parsed";
  case ERR_UNKNOWN_SYMBOL:
    return "unknown symbol in function string";
  case ERR_EXPECT_EQUAL:
    return "expected symbol: '='";
  case ERR_EXPECT_OPEN_P:
    return "expected symbol: '('";
  case ERR_EXPECT_CLOSE_P:
    return "expected symbol: ')'";
  case ERR_INVALID_N_PARAMS:
    return "invalid number of parameters for function";
  case ERR_EXPECT_FUNCT:
    return "function (name) expected";
  case ERR_EXPECT_VAR:
    return "user identifier (variable name) expected";
  case ERR_MISSING:
    return "more tokens expected";

  default:
    return "";
  }
} /* end of _unur_fstr_error_code() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_fstr_debug_input ( const char *fstr )
     /*----------------------------------------------------------------------*/
     /* Print function string on output stream                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   fstr ... string containing function definition                     */
     /*----------------------------------------------------------------------*/
{
  FILE *log = unur_get_stream();

  fprintf(log,"%s: Input string:\n",GENTYPE);
  fprintf(log,"%s:   %s\n",GENTYPE,fstr);
  fprintf(log,"%s:\n",GENTYPE);

} /* end of _unur_fstr_debug_input() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_debug_token ( const struct parser_data *pdata )
     /*----------------------------------------------------------------------*/
     /* Print tokenized string on output stream                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*----------------------------------------------------------------------*/
{
  FILE *log = unur_get_stream();
  int i;

  fprintf(log,"%s: Tokenized string (token separated by blanks):\n",GENTYPE);
  fprintf(log,"%s:   ",GENTYPE);

  for (i=0; i<pdata->n_tokens; i++)
    fprintf(log,"%s ",pdata->tpos[i]);
  fprintf(log,"\n%s:\n",GENTYPE);

} /* end of _unur_fstr_debug_input() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_debug_tree( const struct parser_data *pdata,
		       const struct ftreenode *root )
     /*----------------------------------------------------------------------*/
     /* Print function tree                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata ... pointer to parser object                                 */
     /*   root  ... pointer to root of tree                                  */
     /*----------------------------------------------------------------------*/

{
  FILE *log = unur_get_stream();

  fprintf(log,"%s: parse tree:  (left nodes above right nodes)\n",GENTYPE); 
  fprintf(log,"%s:\n",GENTYPE);
  _unur_fstr_debug_show_tree(pdata,root,0,0);
  fprintf(log,"%s:\n",GENTYPE);

} /* end of _unur_fstr_debug_tree() */

/*---------------------------------------------------------------------------*/

void
_unur_fstr_debug_show_tree(const struct parser_data *pdata,
			   const struct ftreenode *node,
			   int level, int location)
     /*----------------------------------------------------------------------*/
     /* Print function tree by recursion                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   pdata    ... pointer to parser object                              */
     /*   node     ... pointer to node in tree                               */
     /*   level    ... depth in tree                                         */
     /*   location ... indicates location of node in tree by bit array       */
     /*----------------------------------------------------------------------*/
{ 
  FILE *log = unur_get_stream();
  char *name;
  int i, mask; 

  fprintf(log,"%s: ",GENTYPE); 

  /* draw vertical lines in tree */
  for (i = 0, mask = 1; i < level; i++, mask <<= 1) 
    if (mask & location) 
      fprintf(log,"|   "); 
    else 
      fprintf(log,"    "); 

  /* print node */
  if( node != NULL ) {
    
    /* draw horizontal line in tree */
    (mask & location) ? fprintf(log,"+--") : fprintf(log,"\\__");

    /* print symbol */
    switch (node->type) {
    case S_SCONST:
      fprintf(log,"'%s'\t(const=%g)", node->symbol,node->val);  break;
    case S_UCONST:
      fprintf(log,"'%g'\t(const)", node->val);  break;
    case S_UIDENT:
      name = (pdata) ? pdata->variable_name : "x";
      fprintf(log,"'%s'\t(variable)", name);  break;
    case S_UFUNCT:
      name = (pdata) ? pdata->function_name : "f";
      fprintf(log,"'%s'\t(user function)", name);  break;
    default:
      fprintf(log,"'%s'", node->symbol);
    }
    /* end of line */
    fprintf(log,"\n");

    /* print left and right node */
    if ( node->left || node->right) {
      /* ... unless both leaves are empty */
      _unur_fstr_debug_show_tree(pdata,node->left, level+1,location|(mask<<1)); 
      _unur_fstr_debug_show_tree(pdata,node->right,level+1,location); 
    }
  }

  else {  /* empty leave */
    (mask & location) ? fprintf(log,"+--") : fprintf(log,"\\__");
    fprintf(log,"(void)\n"); 
  } 
} /* end of _unur_fstr_debug_show_tree() */


/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** End                                                                     **/
/*****************************************************************************/






/*****************************************************************************/
/** BAUSTELLE                                                               **/
/*****************************************************************************/


/**************** ****************** ****************** *****************/

char *d_dummy() { return NULL; }
char *d_const() { return NULL; }
char *d_add() { return NULL; }
char *d_mul() { return NULL; }
char *d_div() { return NULL; }
char *d_power() { return NULL; }
char *d_exp() { return NULL; }
char *d_ln() { return NULL; }
char *d_log() { return NULL; }
char *d_sin() { return NULL; }
char *d_cos() { return NULL; }
char *d_tan() { return NULL; }
char *d_sec() { return NULL; }
char *d_sqr() { return NULL; }
char *d_abs() { return NULL; }

/**************** ****************** ****************** *****************/
/**************** ****************** ****************** *****************/
