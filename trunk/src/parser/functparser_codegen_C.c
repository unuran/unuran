/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_codegen_C.c                                      *
 *                                                                           *
 *   Make C code for function given by its tree.                             *
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
/** API                                                                     **/
/*****************************************************************************/

/*--------------------------------------------------------------------------*/

int 
_unur_fstr_tree2C ( FILE *out, const struct ftreenode *root,
		    const char *variable, const char *function )
     /*----------------------------------------------------------------------*/
     /* Produce string from function tree.                                   */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*   function ... pointer to name of function                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to output string (should be freed when not used any more)  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE,root,0 );
  _unur_check_NULL( GENTYPE,symbol[root->token].node2C,0 );

  /* make body of C routine */
  symbol[root->token].node2C (out,root,variable);

  return 1;

} /* end of _unur_fstr_tree2C() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Create function string and source code.                                 **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
C_error ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  return 0;
} /* end of C_error() */

/*---------------------------------------------------------------------------*/

int
C_const ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"%.16g",node->val );
  return 1;
} /* end of C_const() */

/*---------------------------------------------------------------------------*/

int
C_var ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"%s",variable);
  return 1;
} /* end of C_var() */

/*---------------------------------------------------------------------------*/

int
C_prefix ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
  /* prefix operators (functions) */
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  /* node '(' left ',' right ')' */
  fprintf(out,"%s(",symbol[node->token].name);

  if (left) {
    symbol[left->token].node2C (out,left,variable);
    fprintf(out,",");
  }
  if (right) {
    symbol[right->token].node2C (out,right,variable);
  }

  fprintf(out,")");

} /* end of C_prefix() */

/*---------------------------------------------------------------------------*/

int
C_power ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
  /* prefix operators (functions) */
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* 'pow(' left ',' right ')' */
  fprintf(out,"pow(");
  symbol[left->token].node2C (out,left,variable);
  fprintf(out,",");
  symbol[right->token].node2C (out,right,variable);
  fprintf(out,")");

} /* end of C_power() */

/*---------------------------------------------------------------------------*/

int
C_infix ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
  /* infix operators (functions) */
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* '(' left node right ')' */

  /* enclosing parenthesis */
  fprintf(out,"(");

  /* left branch */
  symbol[left->token].node2C (out,left,variable);

  /* symbol for node */
  fprintf(out,"%s",symbol[node->token].name);
  
  /* right branch */
  symbol[right->token].node2C (out,right,variable);
  
  /* enclosing parenthesis */
  fprintf(out,")");

} /* end of C_infix() */

/*---------------------------------------------------------------------------*/

int
C_equal ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* '(' left '==' right ')' */
  fprintf(out,"(");
  symbol[left->token].node2C (out,left,variable);
  fprintf(out,"==");
  symbol[right->token].node2C (out,right,variable);
  fprintf(out,")");

} /* end of C_equal() */

/*---------------------------------------------------------------------------*/

int
C_unequal ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* '(' left '==' right ')' */
  fprintf(out,"(");
  symbol[left->token].node2C (out,left,variable);
  fprintf(out,"!=");
  symbol[right->token].node2C (out,right,variable);
  fprintf(out,")");

} /* end of C_unequal() */

/*---------------------------------------------------------------------------*/

int
C_minus ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
  /* infix operators (functions) */
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* '(' [left] '-' right ')' */

  fprintf(out,"(");
  if (!(left->type == S_UCONST && left->val == 0.))
    /* there is no need to print "0 - ..." */
    symbol[left->token].node2C (out,left,variable);
  fprintf(out,"-");
  symbol[right->token].node2C (out,right,variable);
  fprintf(out,")");

} /* end of C_minus() */

/*---------------------------------------------------------------------------*/

int
C_mod ( FILE *out, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
  /* infix operators (functions) */
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return 0;

  /* '(' [left] '%' right ')' */

  fprintf(out,"(");
  symbol[left->token].node2C (out,left,variable);
  fprintf(out,"%%");
  symbol[right->token].node2C (out,right,variable);
  fprintf(out,")");

} /* end of C_mod() */

/*---------------------------------------------------------------------------*/

