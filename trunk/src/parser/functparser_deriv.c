/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_deriv.c                                          *
 *                                                                           *
 *   Compute function tree for derivative.                                   *
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

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_make_derivative ( const struct ftreenode *root )
     /*----------------------------------------------------------------------*/
     /* Make function tree for derivate of given function (tree).            */ 
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to tree of derivative                                      */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *deriv = NULL;

  /* check arguments */
  _unur_check_NULL( GENTYPE,root,NULL );

  deriv = (root) ? (*symbol[root->token].dcalc)(root) : NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into log file */
  if (_unur_default_debugflag)
    _unur_fstr_debug_deriv(root,deriv);
#endif

  return deriv;
} /* end of _unur_fstr_make_derivative() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for computing derivatives                                      **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_dummy (const struct ftreenode *node)
{
  return NULL;
} /* end of d_error() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_const (const struct ftreenode *node)
     /* (const)' = 0                                                         */
     /*                                                                      */
     /*       Const                   0.                                     */
     /*      /     \       ==>       /  \                                    */
     /*  NULL       NULL         NULL    NULL                                */
{
  return _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL);
} /* end of d_const() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_var (const struct ftreenode *node)
     /* x' = 1                                                               */
     /*                                                                      */
     /*       Var                   1.                                       */
     /*      /   \       ==>       /  \                                      */
     /*  NULL     NULL         NULL    NULL                                  */
{
  return _unur_fstr_create_node(NULL,1.,s_uconst,NULL,NULL);
} /* end of d_var() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_add (const struct ftreenode *node)
     /* summation rule:  (l+r)' = l' + r'                                    */
     /*                                                                      */
     /*    Op             Op                                                 */
     /*   /  \    ==>    /  \         (Op ::= '+' | '-')                     */
     /*  X    Y         X'   Y'                                              */
{
  struct ftreenode *left = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *d_left, *d_right;

  /* token for operator */
  int op = node->token;

  /* derivative of both branches */
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make subtree */
  return _unur_fstr_create_node(node->symbol,0.,op,d_left,d_right);
} /* end of d_add() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_mul (const struct ftreenode *node)
     /* product rule:  (l*r)' = l'*r + l*r'                                  */
     /*                                                                      */
     /*    '*'              __'+'__                                          */
     /*   /   \            /       \                                         */
     /*  X     Y   ==>   '*'       '*'                                       */
     /*                 /   \     /   \                                      */
     /*                X'    Y   X     Y'                                    */
{
  struct ftreenode *left = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_left, *br_right;

  /* derivative of both branches */
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of the branches of node */
  left  = _unur_fstr_dup_tree(node->left);
  right = _unur_fstr_dup_tree(node->right);

  /* make subtree */
  br_left  = _unur_fstr_create_node("*",0.,s_mul,d_left,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,left,d_right);

  /* return subtree */
  return _unur_fstr_create_node("+",0.,s_plus,br_left,br_right);
} /* end of d_mul() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_div (const struct ftreenode *node)
     /* Quotient rule:  (l/r)' = (l'*r-l*r')/r^2                             */
     /*                                                                      */
     /*    '/'              ___'/'___                                        */
     /*   /   \            /         \                                       */
     /*  X     Y   ==>   '-'         '^'                                     */
     /*                 /   \       /   \                                    */
     /*               '*'   '*'    Y     2                                   */
     /*               / \   / \                                              */
     /*              X'  Y X   Y'                                            */
     /*                                                                      */
{
  struct ftreenode *left = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_left, *br_right, *two;
  struct ftreenode *numerator, *denominator; 

  /* derivative of both branches */
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of the branches of node */
  left  = _unur_fstr_dup_tree(node->left);
  right = _unur_fstr_dup_tree(node->right);

  /* make nominator */
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   /* const 2 */
  denominator = _unur_fstr_create_node("^",0.,s_power,right,two);
  
  /* make numerator */
  right = _unur_fstr_dup_tree(node->right);    /* we need another copy */
  br_left  = _unur_fstr_create_node("*",0.,s_mul,d_left,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,left,d_right);
  numerator= _unur_fstr_create_node("-",0.,s_minus,br_left,br_right);

  /* subtree */
  return _unur_fstr_create_node("/",0.,s_div,numerator,denominator);
} /* end of d_div() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_power (const struct ftreenode *node)
     /* case: r constant                                                     */
     /* (l^r)' = r * l^(r-1) * l'                                            */
     /*                                                                      */
     /*    '^'                  '*'                                          */
     /*   /   \       ==>      /   \                                         */
     /*  X     Y              X'   '*'                                       */
     /*       / \                 /   \                                      */
     /*   NULL   NULL            Y    '^'                                    */
     /*                              /   \                                   */
     /*                             X    (Y-1)                               */
     /*                                                                      */
     /* case: l constant                                                     */
     /* (l^r)' = r' * l^r * ln(l)                                            */
     /*                                                                      */
     /*        '^'                '*'                                        */
     /*       /   \     ==>      /   \                                       */
     /*      X     Y            Y'   _'*'_                                   */
     /*     / \                     /     \                                  */
     /*  NULL  NULL              "ln"     '^'                                */
     /*                          /  \    /   \                               */
     /*                       NULL   X  X     Y                              */
     /*                                                                      */
     /* otherwise:                                                           */
     /* (l^r)' = l^(r-1) * ( r * l' + l * ln(l) * r' )                       */
     /*                                                                      */
     /*    '^'                  ______'*'_____                               */
     /*   /   \       ==>      /              \                              */
     /*  X     Y             '^'            __'+'__                          */
     /*                     /   \          /       \                         */
     /*                    X    '-'      '*'       '*'                       */
     /*                        /   \    /   \     /   \                      */
     /*                       Y     1  Y     X'  X    '*'                    */
     /*                                              /   \                   */
     /*                                             Y'   "ln"                */
     /*                                                 /    \               */
     /*                                             NULL      X              */
     /*                                                                      */
{
  struct ftreenode *left = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_right;
  struct ftreenode *dup_node;

  if (right->type == S_UCONST || right->type == S_SCONST ) {
     /*    '^'                  '*'                                          */
     /*   /   \       ==>      /   \                                         */
     /*  X     Y              X'   '*'                                       */
     /*       / \                 /   \                                      */
     /*   NULL   NULL            Y    '^'                                    */
     /*                              /   \                                   */
     /*                             X    (Y-1)                               */
    /* derivative of left branch */
    d_left  = (left)  ? (*symbol[left->token].dcalc) (left)  : NULL;
    /* make a copy of the branches of node */
    left  = _unur_fstr_dup_tree(node->left);
    right = _unur_fstr_dup_tree(node->right);
    /* make right branch */
    br_right = _unur_fstr_create_node(NULL,right->val-1,s_uconst,NULL,NULL);
    br_right = _unur_fstr_create_node("^",0.,s_power,left,br_right);
    br_right = _unur_fstr_create_node("*",0.,s_mul,right,br_right);
    /* subtree */
    return _unur_fstr_create_node("*",0.,s_mul,d_left,br_right);
  }

  else if (left->type == S_UCONST || left->type == S_SCONST ) {
     /*        '^'                '*'                                        */
     /*       /   \     ==>      /   \                                       */
     /*      X     Y            Y'   _'*'_                                   */
     /*     / \                     /     \                                  */
     /*  NULL  NULL              "ln"     '^'                                */
     /*                          /  \    /   \                               */
     /*                       NULL   X  X     Y                              */
     /*                                                                      */
    /* find symbol "ln" */
    int s_ln = _unur_fstr_find_symbol("ln",_ans_start,_ans_end);

    /* derivative of right branch */
    d_right = (right) ? (*symbol[right->token].dcalc) (right)  : NULL;
    /* make copies of branches */
    left = _unur_fstr_dup_tree(node->left);
    dup_node = _unur_fstr_dup_tree(node);
    /* make right branch */
    br_right = _unur_fstr_create_node("ln",0.,s_ln,NULL,left);
    br_right = _unur_fstr_create_node("*",0.,s_mul,br_right,dup_node);
    /* subtree */
    return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
  }

  else {
     /*    '^'                  ______'*'_____                               */
     /*   /   \       ==>      /              \                              */
     /*  X     Y             '^'            __'+'__                          */
     /*                     /   \          /       \                         */
     /*                    X    '-'      '*'       '*'                       */
     /*                        /   \    /   \     /   \                      */
     /*                       Y     1  Y     X'  X    '*'                    */
     /*                                              /   \                   */
     /*                                             Y'   "ln"                */
     /*                                                 /    \               */
     /*                                             NULL      X              */
     /*                                                                      */
    /** TODO **/
    return NULL;
  }
} /* end of d_power() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_exp (const struct ftreenode *node)
     /* (exp(r))' = r' * exp(r)                                              */
     /*                                                                      */
     /*     "exp"            '*'                                             */
     /*     /   \    ==>    /   \                                            */
     /*  NULL    X         X'  "exp"                                         */
     /*                        /   \                                         */
     /*                     NULL    X                                        */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of whole subrtree at node */
  right = _unur_fstr_dup_tree(node);

  /* subtree */
  return _unur_fstr_create_node("*",0.,s_mul,d_right,right);
} /* ebd of d_exp() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_ln (const struct ftreenode *node)
     /* (ln(r))' = r'/r                                                      */
     /*                                                                      */
     /*     "ln"            '/'                                              */
     /*     /  \    ==>    /   \                                             */
     /*  NULL   X         X'    X                                            */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* subtree */
  return _unur_fstr_create_node("/",0.,s_div,d_right,right);
} /* end of d_ln() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_log (const struct ftreenode *node)
     /* case: l = const                                                      */
     /* (log(l,r))' = r' / (r * ln(l))                                       */
     /*                                                                      */
     /*     "log"            '/'                                             */
     /*     /   \    ==>    /   \                                            */
     /*    X     Y         Y'   '*'                                          */
     /*                        /   \                                         */
     /*                       Y    "ln"                                      */
     /*                            /  \                                      */
     /*                         NULL   X                                     */
     /*                                                                      */
     /* otherwise:                                                           */
     /*  (TODO)                                                              */
{
  struct ftreenode *left = node->left;
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right, *sub_right;

  if (left->type == S_UCONST || left->type == S_SCONST ) {
    /* find symbol "ln" */
    int s_ln = _unur_fstr_find_symbol("ln",_ans_start,_ans_end);
    /* derivative of right branch */
    d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;
    /* make a copies of both branches of node */
    left  = _unur_fstr_dup_tree(node->left);
    right = _unur_fstr_dup_tree(node->right);
    /* right branch of new tree */
    sub_right = _unur_fstr_create_node("ln",0.,s_ln,NULL,left);
    br_right = _unur_fstr_create_node("*",0.,s_mul,right,sub_right);
    /* subtree */
    return _unur_fstr_create_node("/",0.,s_div,d_right,br_right);
  }

  else {
    /** TODO **/
    return NULL;
  }
} /* end of d_log() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_sin (const struct ftreenode *node)
     /* (sin(r))' = r' * cos(r)                                              */
     /*                                                                      */
     /*     "sin"            '*'                                             */
     /*     /   \    ==>    /   \                                            */
     /*  NULL    X         X'  "cos"                                         */
     /*                        /   \                                         */
     /*                     NULL    X                                        */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;

  /* find symbol "cos" */
  int s_cos = _unur_fstr_find_symbol("cos",_ans_start,_ans_end);

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* right branch of new tree */
  br_right = _unur_fstr_create_node("cos",0.,s_cos,NULL,right);

  /* subtree */
  return _unur_fstr_create_node(NULL,0.,s_mul,d_right,br_right);
} /* end of d_sin() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_cos (const struct ftreenode *node)
     /* (cos(r))' = -r' * sin(r)                                             */
     /*                                                                      */
     /*     "cos"            ____'*'____                                     */
     /*     /   \    ==>    /           \                                    */
     /*  NULL    X        '-'          "sin"                                 */
     /*                  /   \         /   \                                 */
     /*                 0     X'    NULL    X                                */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_left, *br_right;
  struct ftreenode *zero;

  /* find symbol "cos" */
  int s_sin = _unur_fstr_find_symbol("sin",_ans_start,_ans_end);

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* branches of new tree */
  br_right = _unur_fstr_create_node("sin",0.,s_sin,NULL,right);

  zero = _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL);
  br_left  = _unur_fstr_create_node("-",0.,s_minus,zero,d_right);

  /* subtree */
  return _unur_fstr_create_node("*",0.,s_mul,br_left,br_right);
} /* end of d_cos() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_tan (const struct ftreenode *node)
     /* (tan(r))' = r' * (sec(r))^2                                          */
     /*                                                                      */
     /*     "tan"            '*'                                             */
     /*     /   \    ==>    /   \                                            */
     /*  NULL    X         X'   '^'                                          */
     /*                        /   \                                         */
     /*                    "sec"    2.                                       */
     /*                    /   \                                             */
     /*                 NULL    X                                            */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right, *sub_right;
  struct ftreenode *two;

  /* find symbol "sec" */
  int s_sec = _unur_fstr_find_symbol("sec",_ans_start,_ans_end);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* right branch of new tree */
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   /* const 2 */
  sub_right = _unur_fstr_create_node("sec",0.,s_sec,NULL,right);
  br_right = _unur_fstr_create_node("^",0.,s_power,sub_right,two);

  /* subtree */
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} /* end of d_tan() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_sec (const struct ftreenode *node)
     /* (sec(r))' = r' * tan(r) * sec(r) */
     /*                                                                      */
     /*     "sec"            '*'                                             */
     /*     /   \    ==>    /   \                                            */
     /*  NULL    X         X'   '*'                                          */
     /*                        /   \                                         */
     /*                   "tan"     "sec"                                    */
     /*                   /   \     /   \                                    */
     /*                NULL    X  NULL   X                                   */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right, *sub_right, *dup_node;

  /* find symbols "sec" and "tan" */
  int s_tan = _unur_fstr_find_symbol("tan",_ans_start,_ans_end);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* make a copy of the whole subtree at node */
  dup_node = _unur_fstr_dup_tree(node);

  /* right branch of new tree */
  sub_right = _unur_fstr_create_node("tan",0.,s_tan,NULL,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,sub_right,dup_node);

  /* subtree */
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} /* end of d_sec() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_sqrt (const struct ftreenode *node)
     /* (sqrt(r))' = r'/(2*sqrt(r))                                          */
     /*                                                                      */
     /*     "sqrt"            '/'                                            */
     /*     /    \    ==>    /   \                                           */
     /*  NULL     X         X'   '*'                                         */
     /*                         /   \                                        */
     /*                        2.   "sqrt"                                   */
     /*                             /    \                                   */
     /*                          NULL     X                                  */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;
  struct ftreenode *two;

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* make a copy of the whole subtree at node */
  right = _unur_fstr_dup_tree(node);

  /* right branch of new tree */
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   /* const 2 */
  br_right = _unur_fstr_create_node("*",0.,s_mul,two,right);

  /* subtree */
  return _unur_fstr_create_node("/",0.,s_div,d_right,br_right);
} /* end of d_sqrt() */

/*---------------------------------------------------------------------------*/

struct ftreenode *
d_abs (const struct ftreenode *node)
     /* (abs(r))' = r' * sgn(x)                                              */
     /*                                                                      */
     /*     "abs"            '*'                                             */
     /*     /   \    ==>    /   \                                            */
     /*  NULL    X         X'  "sgn"                                         */
     /*                        /   \                                         */
     /*                     NULL    X                                        */
{
  struct ftreenode *right = node->right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;

  /* find symbol "cos" */
  int s_sgn = _unur_fstr_find_symbol("sgn",_ans_start,_ans_end);

  /* make a copy of the right branch of node */
  right = _unur_fstr_dup_tree(node->right);

  /* derivative of right branch */
  d_right = (right) ? (*symbol[right->token].dcalc)(right) : NULL;

  /* right branch of new tree */
  br_right = _unur_fstr_create_node("sgn",0.,s_sgn,NULL,right);

  /* subtree */
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} /* end of d_abs() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary rRoutines                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct ftreenode *
_unur_fstr_dup_tree (const struct ftreenode *root)
     /*----------------------------------------------------------------------*/
     /* Duplicate function tree rooted at root                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   root ... pointer to root of function tree                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to duplicated tree                                         */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *dup;

  if (root==NULL) return NULL;

  dup = _unur_malloc(sizeof(struct ftreenode));
  memcpy(dup,root,sizeof(struct ftreenode));
  if (root->left)  dup->left  = _unur_fstr_dup_tree(root->left);
  if (root->right) dup->right = _unur_fstr_dup_tree(root->right);

  return dup;

} /* end of _unur_fstr_dup_tree() */

/*---------------------------------------------------------------------------*/

