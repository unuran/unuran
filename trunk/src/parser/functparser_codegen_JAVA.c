/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_codegen_JAVA.c                                   *
 *                                                                           *
 *   Make JAVA code for function given by its tree.                          *
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
/* Codes for special JAVA functions (must be added to code)                  */
/* (use bits in integer)                                                     */

enum {
  J_FUNCT_ERROR = 0x80000000u,      /* error                                 */

  J_FUNCT_SGN   = 0x00000001u,      /* sign function                         */
  J_FUNCT_SEC   = 0x00000002u,      /* secant function                       */

  J_FUNCT_LT    = 0x00000010u,      /* '<'  function                         */
  J_FUNCT_LE    = 0x00000020u,      /* '<=' function                         */
  J_FUNCT_GT    = 0x00000040u,      /* '>'  function                         */
  J_FUNCT_GE    = 0x00000080u,      /* '>=' function                         */
  J_FUNCT_EQ    = 0x00000100u,      /* '==' function                         */
  J_FUNCT_NE    = 0x00000200u,      /* '!=' function                         */
};

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** API                                                                     **/
/*****************************************************************************/

/*--------------------------------------------------------------------------*/
/* we can use the print function for C symbols                              */

#define _unur_fstr_print_J(output,symb,number)   _unur_fstr_print_C((output),(symb),(number)) 

/*--------------------------------------------------------------------------*/

int 
_unur_fstr_tree2JAVA ( FILE *out, const struct ftreenode *root,
		       const char *variable, const char *funct_name )
     /*----------------------------------------------------------------------*/
     /* Produce string from function tree.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out        ... output stream                                       */
     /*   root       ... pointer to root of function tree                    */
     /*   variable   ... pointer to name of variable                         */
     /*   funct_name ... pointer to name of JAVA function                    */
     /*                                                                      */
     /* return:                                                              */
     /*   1 ... success                                                      */
     /*   0 ... failure                                                      */
     /*----------------------------------------------------------------------*/
{
  struct concat output = {NULL, 0, 0};
  unsigned rcode;

  /* check arguments */
  _unur_check_NULL( GENTYPE,root,0 );
  _unur_check_NULL( GENTYPE,symbol[root->token].node2J,0 );

  /* make body of JAVA routine */
  rcode = symbol[root->token].node2J (&output,root,variable);
  if (rcode & J_FUNCT_ERROR) { 
    if (output.string) free(output.string);
    return 0;
  }

  /* print code for special functions (if necessary) */ 
  _unur_fstr_J_specfunct (out,rcode);

  /* print JAVA routine */
  *(output.string + output.length) = '\0';
  fprintf (out,"\tstatic double %s (double %s)\n",funct_name,variable );
  fprintf (out,"\t{\n\t\treturn (%s);\n\t}\n",output.string);

  /* free memory */
  free(output.string);

  return 1;

} /* end of _unur_fstr_tree2JAVA() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Create function string and source code.                                 **/
/*****************************************************************************/



/* J_xxxx ( struct concat *output, const struct ftreenode *node, const char *variable ) */
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
J_error ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Error (This should not happen).                                      */
     /*----------------------------------------------------------------------*/
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return J_FUNCT_ERROR;
} /* end of J_error() */

/*---------------------------------------------------------------------------*/

unsigned
J_const ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for constant                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_J( output, NULL, node->val );
  return 0u;
} /* end of J_const() */

/*---------------------------------------------------------------------------*/

unsigned
J_var ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for variable                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_J( output, variable, 0 );
  return 0u;
} /* end of J_var() */

/*---------------------------------------------------------------------------*/

unsigned
J_prefix_generic ( struct concat *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
     /*----------------------------------------------------------------------*/
     /* print prefix operator (function). generic version                    */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  /* node '(' left ',' right ')' */
  _unur_fstr_print_J( output, symb, 0 );
  _unur_fstr_print_J( output, "(", 0 );

  if (left) {
    rcode |= symbol[left->token].node2J (output,left,variable);
    _unur_fstr_print_J( output, ",", 0 );
  }
  if (right) {
    rcode |= symbol[right->token].node2J (output,right,variable);
  }

  _unur_fstr_print_J( output, ")", 0 );

  return rcode;
} /* end of J_prefix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
J_prefix ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for prefix operator (function)                                */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_J( output, "Math.", 0 );
  return J_prefix_generic(output,node,variable,symbol[node->token].name);
} /* end of J_prefix() */

/*---------------------------------------------------------------------------*/

unsigned
J_lt ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '<' function                                              */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_LT | J_prefix_generic(output,node,variable,"RelLT"));
} /* end of J_lt() */

/*---------------------------------------------------------------------------*/

unsigned
J_le ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '<=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_LE | J_prefix_generic(output,node,variable,"RelLE"));
} /* end of J_le() */

/*---------------------------------------------------------------------------*/

unsigned
J_gt ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '>' function                                              */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_GT | J_prefix_generic(output,node,variable,"RelGT"));
} /* end of J_gt() */

/*---------------------------------------------------------------------------*/

unsigned
J_ge ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '>=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_GE | J_prefix_generic(output,node,variable,"RelGE"));
} /* end of J_ge() */

/*---------------------------------------------------------------------------*/

unsigned
J_eq ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '==' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_EQ | J_prefix_generic(output,node,variable,"RelEQ"));
} /* end of J_eq() */

/*---------------------------------------------------------------------------*/

unsigned
J_ne ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for '!=' function                                             */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_NE | J_prefix_generic(output,node,variable,"RelNE"));
} /* end of J_ne() */

/*---------------------------------------------------------------------------*/

unsigned
J_power ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for power function                                            */
     /*----------------------------------------------------------------------*/
{
  return J_prefix_generic(output,node,variable,"Math.pow");
} /* end of J_power() */

/*---------------------------------------------------------------------------*/

unsigned
J_sec ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for secant function                                           */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_SEC | J_prefix_generic(output,node,variable,"sec"));
} /* end of J_sec() */

/*---------------------------------------------------------------------------*/

unsigned
J_sgn ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for sign function                                            */
     /*----------------------------------------------------------------------*/
{
  return (J_FUNCT_SGN | J_prefix_generic(output,node,variable,"sgn"));
} /* end of J_sgn() */

/*---------------------------------------------------------------------------*/

unsigned
J_infix_generic ( struct concat *output, const struct ftreenode *node,
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
    return J_FUNCT_ERROR;

  /* '(' left node right ')' */
  _unur_fstr_print_J( output, "(", 0 );
  rcode |= symbol[left->token].node2J (output,left,variable);
  _unur_fstr_print_J( output, symb, 0 );
  rcode |= symbol[right->token].node2J (output,right,variable);
  _unur_fstr_print_J( output, ")", 0 );

  return rcode;

} /* end of J_infix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
J_infix ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for infix operator (binary operator).                         */
     /*----------------------------------------------------------------------*/
{
  return J_infix_generic(output,node,variable,symbol[node->token].name);
} /* end of J_infix() */

/*---------------------------------------------------------------------------*/

unsigned
J_minus ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for minus operator                                            */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return J_FUNCT_ERROR;

  /* '(' [left] '-' right ')' */
  _unur_fstr_print_J( output, "(", 0 );
  if (!(left->type == S_UCONST && left->val == 0.))
    /* there is no need to print "0 - ..." */
    rcode |= symbol[left->token].node2J (output,left,variable);
  _unur_fstr_print_J( output, "-", 0 );
  rcode |= symbol[right->token].node2J (output,right,variable);
  _unur_fstr_print_J( output, ")", 0 );

  return rcode;
} /* end of J_minus() */

/*---------------------------------------------------------------------------*/

unsigned
J_mod ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for mudolus function                                          */
     /*----------------------------------------------------------------------*/
{
  return J_infix_generic(output,node,variable,"%");
} /* end of J_mod() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Special functions                                                       **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_J_specfunct ( FILE *out, unsigned flags )
     /*----------------------------------------------------------------------*/
     /* Print FORTRAN code for special functions                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out   ... output stream                                            */
     /*   flags ... bit array with flags for special functions               */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  if (flags & J_FUNCT_SGN) {
    _unur_fstr_J_sgn(out);
  }

  if (flags & J_FUNCT_SEC) {
    fprintf(out,"\tstatic double sec (double x) { ");
    fprintf(out,"return (1./Math.cos(x)); }\n\n");
  }

  if (flags & J_FUNCT_LE) {
    fprintf(out,"\tstatic double RelLE (double x, double y) { ");
    fprintf(out,"return ((x<=y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_GE) {
    fprintf(out,"\tstatic double RelGE (double x, double y) { ");
    fprintf(out,"return ((x>=y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_LT) {
    fprintf(out,"\tstatic double RelLT (double x, double y) { ");
    fprintf(out,"return ((x<y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_GT) {
    fprintf(out,"\tstatic double RelGT (double x, double y) { ");
    fprintf(out,"return ((x>y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_EQ) {
    fprintf(out,"\tstatic double RelEQ (double x, double y) { ");
    fprintf(out,"return ((x==y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_NE) {
    fprintf(out,"\tstatic double RelNE (double x, double y) { ");
    fprintf(out,"return ((x!=y) ? 1. : 0.); }\n\n");
  }

  return 1;
} /* end of _unur_fstr_J_specfunct() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_J_sgn ( FILE *out )
     /*----------------------------------------------------------------------*/
     /* Print FORTRAN code for special functions                                   */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"\tstatic double sgn (double x)\n\t{\n");
  fprintf(out,"\t\tif (x<0.)  return -1.;\n");
  fprintf(out,"\t\tif (x>0.)  return  1.;\n");
  fprintf(out,"\t\t/* else */ return  0.;\n");
  fprintf(out,"\t}\n\n");

  return 1;
} /* end of _unur_fstr_J_sgn() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
