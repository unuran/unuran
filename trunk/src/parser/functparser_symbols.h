/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_symbols.c                                        *
 *                                                                           *
 *   Table of known symbols for function parser                              *
 *                                                                           *
 *   Parser function string, evaluate function, print programm code.         *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given a string for a function.                                       *
 *         The string is parser.                                             *
 *         A tree representing the function term is generated.               *
 *         A tree for the derivative of the function is generated.           *
 *         The tree is used to evalute the corresponding function for an x.  *
 *         The source code for a program is produced.                        *
 *                                                                           *
 *****************************************************************************
     $Id$
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000 Wolfgang Hoermann and Josef Leydold                  *
 *   Dept. for Statistics, University of Economics, Vienna, Austria          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

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
static double v_plus   (double l, double r); 
static double v_minus  (double l, double r); 
static double v_mul    (double l, double r);
static double v_div    (double l, double r);
static double v_mod    (double l, double r);
static double v_power  (double l, double r);
static double v_const  (double l, double r);
static double v_exp    (double l, double r);
static double v_ln     (double l, double r);
static double v_log    (double l, double r);
static double v_sin    (double l, double r); 
static double v_cos    (double l, double r);
static double v_tan    (double l, double r);
static double v_sec    (double l, double r); 
static double v_sqrt   (double l, double r);
static double v_abs    (double l, double r);
static double v_sgn    (double l, double r);

/*****************************************************************************/
/** Routines for computing derivatives                                      **/
/*****************************************************************************/

static struct ftreenode *d_error ( const struct ftreenode *node, int *error );
static struct ftreenode *d_const ( const struct ftreenode *node, int *error );
static struct ftreenode *d_var   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_add   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_mul   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_div   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_power ( const struct ftreenode *node, int *error );
static struct ftreenode *d_exp   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_ln    ( const struct ftreenode *node, int *error );
static struct ftreenode *d_log   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sin   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_cos   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_tan   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sec   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sqrt  ( const struct ftreenode *node, int *error );
static struct ftreenode *d_abs   ( const struct ftreenode *node, int *error );

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

#define SYMBLENGTH 10           /* Maximum length for symbols                */

/* structure for storing known symbols                                       */
struct symbols { 
  char   name[SYMBLENGTH];       /* name of symbols (e.g. "sin")             */ 
  int    type;                   /* type of symbol */
  int    info;                   /* priority or 
				    number of argument for system function   */
  double val;                    /* value of constant (0. otherwise)         */

  double (*vcalc)(double l, double r);    /* pointer to function
					     for computing value of node     */
  struct ftreenode *(*dcalc)(const struct ftreenode *node, int *error); 
                                          /* pointer to function  
					     for computing derivate          */

  char   name_C      [SYMBLENGTH];  /* name of C function | operator         */ 
  char   name_FORTRAN[SYMBLENGTH];  /* name of FORTRAN function | operator   */ 

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
/*  Constant should have highest priority to avoid enclosing parenthesis     */
/*  in string output.                                                        */
/*                                                                           */
/*  The first entry in the list must be empty!                               */
/*                                                                           */
/*  There must be FOUR markers in the list (in this ordering!):              */
/*     _ROS   ... start of relation operators                                */
/*     _NAS   ... start of non-alphanumeric symbols (single characters!)     */
/*     _ANS   ... start of alphanumeric symbols (constants and functions)    */
/*     _END   ... end of list                                                */
/*                                                                           */
/*  Prefix operators (functions) for C or FORTRAN are indicated by an        */
/*  '@' sign (eg. "@sin"). Operator names without an '@' sign are treated as */
/*  infix operators (eg. "+").                                               */
/*                                                                           */

static struct symbols symbol[] = {   
  /* symbol,	       priority, evaluation routine, C function,             */
  /*       type,          value,           derivative,   FORTRAN function    */

  /* void */
  {""    , S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, "", "" },

  /* user defined symbols: */                        
  /*    constant           */                     
  {"UCONST",S_UCONST , 9, 0.0 , v_const  , d_const, "", "" },
  /*    function           */
  {"UFUNCT",S_UFUNCT , 0, 0.0 , v_dummy  , d_error, "", "" },
  /*    variable           */
  {"VAR" , S_UIDENT  , 9, 0.0 , v_dummy  , d_var  , "", "" },

  /* marker for relation operators */
  {"_ROS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, "", "" },

  /* relation operators  */
  {"<"   , S_REL_OP  , 1, 0.0 , v_less   , d_const, "<" , ".LT." },
  {"="   , S_REL_OP  , 1, 0.0 , v_equal  , d_const, "==", ".EQ." },
  {"=="  , S_REL_OP  , 1, 0.0 , v_equal  , d_const, "==", ".EQ." },
  {">"   , S_REL_OP  , 1, 0.0 , v_greater, d_const, ">" , ".GT." },
  {"<="  , S_REL_OP  , 1, 0.0 , v_less_or, d_const, "<=", ".LE." },
  {"<>"  , S_REL_OP  , 1, 0.0 , v_unequal, d_const, "!=", ".NE." },
  {"!="  , S_REL_OP  , 1, 0.0 , v_unequal, d_const, "!=", ".NE." },
  {">="  , S_REL_OP  , 1, 0.0 , v_grtr_or, d_const, ">=", ".GE." },

  /* marker for non-alphanumeric symbols */
  {"_NAS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, "", "" },

  /* special symbols */
  {"("   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, "", "" },
  {")"   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, "", "" },
  {","   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, ",", "," },

  /* arithmetic operators */
  {"+"   , S_ADD_OP  , 2, 0.0 , v_plus   , d_add  , "+"   , "+"  },
  {"-"   , S_ADD_OP  , 2, 0.0 , v_minus  , d_add  , "-"   , "-"  },
  {"*"   , S_MUL_OP  , 4, 0.0 , v_mul    , d_mul  , "*"   , "*"  },
  {"/"   , S_MUL_OP  , 4, 0.0 , v_div    , d_div  , "/"   , "/"  },
  {"^"   , S_HPR_OP  , 5, 0.0 , v_power  , d_power, "@pow", "**" },

  /* marker for alphanumeric symbols */
  {"_ANS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, "", "" },

  /* logical operators: removed */

  /* system constants */
  {"pi"  , S_SCONST  , 9, M_PI, v_const  , d_const, "", "" },
  {"e"   , S_SCONST  , 9, M_E , v_const  , d_const, "", "" },

  /* system functions */
  {"mod" , S_SFUNCT  , 2, 0.0 , v_mod    , d_const, "%", "xxxx" },
  {"exp" , S_SFUNCT  , 1, 0.0 , v_exp    , d_exp  , "@exp", "xxxx" },
  {"ln"  , S_SFUNCT  , 1, 0.0 , v_ln     , d_ln   , "@log", "xxxx" },
  {"log" , S_SFUNCT  , 2, 0.0 , v_log    , d_log  , "xxxx", "xxxx" },
  {"sin" , S_SFUNCT  , 1, 0.0 , v_sin    , d_sin  , "@sin", "xxxx" },
  {"cos" , S_SFUNCT  , 1, 0.0 , v_cos    , d_cos  , "@cos", "xxxx" },
  {"tan" , S_SFUNCT  , 1, 0.0 , v_tan    , d_tan  , "xxxx", "xxxx" },
  {"sec" , S_SFUNCT  , 1, 0.0 , v_sec    , d_sec  , "xxxx", "xxxx" },
  {"sqrt", S_SFUNCT  , 1, 0.0 , v_sqrt   , d_sqrt , "@sqrt", "xxxx" },
  {"abs" , S_SFUNCT  , 1, 0.0 , v_abs    , d_abs  , "@fabs", "xxxx" },
  {"sgn" , S_SFUNCT  , 1, 0.0 , v_sgn    , d_const, "xxxx", "xxxx" },

  /* marker for end-of-table */
  {"_END", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, "", "" },
};

/*  location of special symbols in table                                     */
/** DO NOT CHANGE THEIR POSITION IN TABLE **/
static int s_uconst = 1;      /* user defined constant                       */
static int s_ufunct = 2;      /* user defined function                       */
static int s_uident = 3;      /* user defined variable                       */

/* location of special symbols */
static int s_comma, s_minus, s_plus, s_mul, s_div, s_power;

/* Marker for different regions in symbol table                              */
static int _ros_start, _ros_end;    /* relation symbols                      */
static int _nas_start, _nas_end;    /* non-alphanumeric symbols              */
static int _ans_start, _ans_end;    /* alphanumeric symbols                  */
static int _end;                    /* end of list                           */

