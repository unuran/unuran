/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_codegen_FORTRAN.c                                *
 *                                                                           *
 *   Make FORTRAN code for function given by its tree.                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
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
/* Codes for special FORTRAN functions (must be added to code)               */
/* (use bits in integer)                                                     */

enum {
  F_FUNCT_ERROR = 0x10000000u,      /* error                                 */

  F_FUNCT_SGN   = 0x00000001u,      /* sign function                         */
  F_FUNCT_SEC   = 0x00000002u,      /* secant function                       */

  F_FUNCT_LT    = 0x00000010u,      /* '<'  function                         */
  F_FUNCT_LE    = 0x00000020u,      /* '<=' function                         */
  F_FUNCT_GT    = 0x00000040u,      /* '>'  function                         */
  F_FUNCT_GE    = 0x00000080u,      /* '>=' function                         */
  F_FUNCT_EQ    = 0x00000100u,      /* '==' function                         */
  F_FUNCT_NE    = 0x00000200u       /* '!=' function                         */
};

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** API                                                                     **/
/*****************************************************************************/

/*--------------------------------------------------------------------------*/

int 
_unur_fstr_tree2FORTRAN ( FILE *out, const struct ftreenode *root,
			  const char *variable, const char *funct_name )
     /*----------------------------------------------------------------------*/
     /* Produce string from function tree.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out        ... output stream                                       */
     /*   root       ... pointer to root of function tree                    */
     /*   variable   ... pointer to name of variable                         */
     /*   funct_name ... pointer to name of FORTRAN function                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_string output = {NULL, 0, 0};
  unsigned rcode;
  int line;

  /* check arguments */
  _unur_check_NULL( GENTYPE, root, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, symbol[root->token].node2F, UNUR_ERR_NULL );

  /* make body of FORTRAN routine */
  rcode = symbol[root->token].node2F (&output,root,variable);
  if (rcode & F_FUNCT_ERROR) { 
    if (output.text) free(output.text);
    return UNUR_ERR_GEN_DATA;
  }

  /* print FORTRAN routine */
  fprintf (out,"      DOUBLE PRECISION FUNCTION %.6s(x)\n\n", funct_name);
  fprintf (out,"      IMPLICIT DOUBLE PRECISION (A-Z)\n");

  /* print code for special statement functions (if necessary) */ 
  _unur_fstr_F_specfunct (out,rcode);

  fprintf (out,"C\n");
  fprintf (out,"C     compute PDF\n");
  fprintf (out,"C\n");

  /* print body */
  fprintf (out,"      %.6s = \n", funct_name);
  for (line = 0; line < (output.length-1)/60 + 1; line++) {
    fprintf (out,"     $   %.60s\n", output.text+60*line);
  }

  fprintf (out,"      RETURN\n");
  fprintf (out,"\n");
  fprintf (out,"      END\n");
  fprintf (out,"\n");

  /* free memory */
  free(output.text);

  return UNUR_SUCCESS;

} /* end of _unur_fstr_tree2FORTRAN() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Create function string and source code.                                 **/
/*****************************************************************************/

/* F_xxxx ( struct unur_string *output, const struct ftreenode *node, const char *variable ) */
/*---------------------------------------------------------------------------*/
/* Produce string from function subtree rooted at node.                      */
/*                                                                           */
/* parameters:                                                               */
/*   output   ... pointer to string for output                               */
/*   node     ... pointer to root of function tree                           */
/*   variable ... pointer to name of variable                                */
/*                                                                           */
/* return:                                                                   */
/*   Flags for special functions (see above)                                 */
/*---------------------------------------------------------------------------*/

unsigned
F_error ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Error (This should not happen).                                      */
     /*----------------------------------------------------------------------*/
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return F_FUNCT_ERROR;
} /* end of F_error() */

/*---------------------------------------------------------------------------*/

unsigned
F_const ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for constant                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_F( output, NULL, node->val );
  return 0u;
} /* end of F_const() */

/*---------------------------------------------------------------------------*/

unsigned
F_var ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for variable                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_F( output, variable, 0. );
  return 0u;
} /* end of F_var() */

/*---------------------------------------------------------------------------*/

unsigned
F_prefix_generic ( struct unur_string *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
     /*----------------------------------------------------------------------*/
     /* print prefix operator (function). generic version                    */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  /* node '(' left ',' right ')' */
  _unur_fstr_print_F( output, symb, 0. );
  _unur_fstr_print_F( output, "(", 0. );

  if (left) {
    rcode |= symbol[left->token].node2F (output,left,variable);
    _unur_fstr_print_F( output, ",", 0. );
  }
  if (right) {
    rcode |= symbol[right->token].node2F (output,right,variable);
  }

  _unur_fstr_print_F( output, ")", 0. );

  return rcode;
} /* end of F_prefix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
F_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for prefix operator (function)                                */
     /*----------------------------------------------------------------------*/
{
  return F_prefix_generic(output,node,variable,symbol[node->token].name);
} /* end of F_prefix() */

/*---------------------------------------------------------------------------*/

unsigned
F_lt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '<' function                                              */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_LT | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelLT"));
} /* end of F_lt() */

/*---------------------------------------------------------------------------*/

unsigned
F_le ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '<=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_LE | F_prefix_generic(output,node,variable,"RelLE"));
} /* end of F_le() */

/*---------------------------------------------------------------------------*/

unsigned
F_gt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '>' function                                              */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_GT | F_FUNCT_LE | F_prefix_generic(output,node,variable,"RelGT"));
} /* end of F_gt() */

/*---------------------------------------------------------------------------*/

unsigned
F_ge ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '>=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelGE"));
} /* end of F_ge() */

/*---------------------------------------------------------------------------*/

unsigned
F_eq ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '==' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_EQ | F_FUNCT_LE | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelEQ"));
} /* end of F_eq() */

/*---------------------------------------------------------------------------*/

unsigned
F_ne ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '!=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_NE | F_FUNCT_LE | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelNE"));
} /* end of F_ne() */

/*---------------------------------------------------------------------------*/

unsigned
F_sec ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for secant function                                           */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_SEC | F_prefix_generic(output,node,variable,"sec"));
} /* end of F_sec() */

/*---------------------------------------------------------------------------*/

unsigned
F_sgn ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for sign function                                            */
     /*----------------------------------------------------------------------*/
{
  return (F_FUNCT_SGN | F_prefix_generic(output,node,variable,"sgn"));
} /* end of F_sgn() */

/*---------------------------------------------------------------------------*/

unsigned
F_infix_generic ( struct unur_string *output, const struct ftreenode *node,
		  const char *variable, const char *symb )
     /*----------------------------------------------------------------------*/
     /* print infix operator (binary operator). generic version              */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return F_FUNCT_ERROR;

  /* '(' left node right ')' */
  _unur_fstr_print_F( output, "(", 0. );
  rcode |= symbol[left->token].node2F (output,left,variable);
  _unur_fstr_print_F( output, symb, 0. );
  rcode |= symbol[right->token].node2F (output,right,variable);
  _unur_fstr_print_F( output, ")", 0. );

  return rcode;

} /* end of F_infix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
F_infix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for infix operator (binary operator).                         */
     /*----------------------------------------------------------------------*/
{
  return F_infix_generic(output,node,variable,symbol[node->token].name);
} /* end of F_infix() */

/*---------------------------------------------------------------------------*/

/* unsigned */
/* F_equal ( struct unur_string *output, const struct ftreenode *node, const char *variable ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* string for equality operator                                         *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   return F_infix_generic(output,node,variable,"=="); */
/* } /\* end of F_equal() *\/ */

/*---------------------------------------------------------------------------*/

/* unsigned */
/* F_unequal ( struct unur_string *output, const struct ftreenode *node, const char *variable ) */
/*      /\*----------------------------------------------------------------------*\/ */
/*      /\* string for inequality operator                                       *\/ */
/*      /\*----------------------------------------------------------------------*\/ */
/* { */
/*   return F_infix_generic(output,node,variable,"!="); */
/* } /\* end of F_unequal() *\/ */

/*---------------------------------------------------------------------------*/

unsigned
F_minus ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for minus operator                                            */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return F_FUNCT_ERROR;

  /* '(' [left] '-' right ')' */
  _unur_fstr_print_F( output, "(", 0. );
  if (!(left->type == S_UCONST && left->val == 0.))
    /* there is no need to print "0 - ..." */
    rcode |= symbol[left->token].node2F (output,left,variable);
  _unur_fstr_print_F( output, "-", 0. );
  rcode |= symbol[right->token].node2F (output,right,variable);
  _unur_fstr_print_F( output, ")", 0. );

  return rcode;
} /* end of F_minus() */

/*---------------------------------------------------------------------------*/

unsigned
F_power ( struct unur_string *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for power function                                            */
     /*----------------------------------------------------------------------*/
{
  return F_infix_generic(output,node,variable,"**");
} /* end of F_power() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Special functions                                                       **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_F_specfunct ( FILE *out, unsigned flags )
     /*----------------------------------------------------------------------*/
     /* Print FORTRAN code for special functions                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out   ... output stream                                            */
     /*   flags ... bit array with flags for special functions               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (flags & F_FUNCT_SEC) {
    fprintf (out,"      sec(a)=1.d0/cos(a)\n");
  }

  if (flags & F_FUNCT_SGN) {
    fprintf (out,"      sgn(a)=sign(1.d0,a)\n");
  }

  if (flags & F_FUNCT_LE) {
    fprintf (out,"      RelLE(a,b)=sign(0.5d0,b-a)+0.5d0\n");
  }
  if (flags & F_FUNCT_GE) {
    fprintf (out,"      RelGE(a,b)=sign(0.5d0,a-b)+0.5d0\n");
  }
  if (flags & F_FUNCT_LT) {
    fprintf (out,"      RelLT(a,b)=1.d0-RelGE(a,b)\n");
  }
  if (flags & F_FUNCT_GT) {
    fprintf (out,"      RelGT(a,b)=1.d0-RelLE(a,b)\n");
  }
  if (flags & F_FUNCT_EQ) {
    fprintf (out,"      RelEQ(a,b)=RelGE(a,b)*RelLE(a,b)\n");
  }
  if (flags & F_FUNCT_NE) {
    fprintf (out,"      RelNE(a,b)=1.d0-RelGE(a,b)*RelLE(a,b)\n");
  }

  return UNUR_SUCCESS;
} /* end of _unur_fstr_F_specfunct() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_print_F ( struct unur_string *output, const char *symb, const double number )
     /*----------------------------------------------------------------------*/
     /* Print string or number as FORTRAN code into output string.           */
     /* The number is only printed if symb is the NULL pointer.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   output ... pointer to string for output                            */
     /*   symb   ... string to be printed                                    */
     /*   number ... constant to be printed                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  char buf[128];
  char *here_is_e;

  if (symb) {
    /* copy symbol into string */
    _unur_string_appendtext( output, symb );
  }
  else {
    /* make string */
    sprintf(buf,"%.20e",number);
    /* replace `e' by `D' */
    here_is_e = strchr(buf, 'e');
    if (here_is_e) *here_is_e = 'D';
    /* copy number symbol into output */
    _unur_string_appendtext( output, buf );
  }

  return UNUR_SUCCESS;
} /* end of _unur_fstr_print_F() */

/*---------------------------------------------------------------------------*/
