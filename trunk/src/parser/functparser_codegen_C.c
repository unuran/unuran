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

/*---------------------------------------------------------------------------*/

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

  /* make body of C routine */
  _unur_fstr_node2C(out,root,variable,function);

  return 1;

} /* end of _unur_fstr_tree2C() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Create function string and source code.                                 **/
/*****************************************************************************/

int
_unur_fstr_node2C ( FILE *out, const struct ftreenode *node,
		    const char *variable, const char *function )
     /*----------------------------------------------------------------------*/
     /* Produce string from function subtree rooted at node.                 */
     /* As a side effect a string is allocated.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out      ... output stream                                         */
     /*   root     ... pointer to root of function tree                      */
     /*   variable ... pointer to name of variable                           */
     /*   function ... pointer to name of function                           */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */
  const char *symb;                /* symbol for node (or NULL for constant) */           
  int type = node->type;                   /* type of symbol                 */
  int priority = symbol[node->token].info; /* priority of symbol             */
  int operator, parenthesis;               /* booleans                       */

  /* get symbol for node */
  switch (type) {
  case S_UIDENT:    /* variable */
    symb = variable;  break;
  case S_UFUNCT:    /* function */
    symb = function;  break;
  case S_UCONST:    /* node contains constant */
    /* use value in node instead of symbol */
    symb = NULL;      break;  
  case S_SCONST:
  default:
    symb = node->symbol;
  }

  if (type == S_SFUNCT || type == S_UFUNCT) {
    /* node '(' left ',' right ')' */
    _unur_fstr_print_C( out, symb, 0 );
    _unur_fstr_print_C( out, "(", 0. );
    if (left) {
      _unur_fstr_node2C(out,left,variable,function);
      _unur_fstr_print_C( out, ",", 0. );
    }
    if (right) {
      _unur_fstr_node2C(out,right,variable,function);
    }
    _unur_fstr_print_C( out, ")", 0. );
  }

  else if (symb && symb[0] == ',') {
    /* left ',' right */
    _unur_fstr_print_C( out, ",", 0. );
    if (left) {
      _unur_fstr_node2C(out,left,variable,function);
      _unur_fstr_print_C( out, ",", 0. );
    }
    if (right) {
      _unur_fstr_node2C(out,right,variable,function);
    }
  }    

  else {
    /* check whether enclosing blanks are required for typography  */
    operator = (type==S_REL_OP || type==S_ADD_OP || type==S_MUL_OP);

    /* left branch */
    if (left) {
      /* always use parenthesis ... */
      parenthesis = 1;
      /* ... except ... */
      if (left->type == S_SCONST || left->type == S_UCONST || 
	  left->type == S_SFUNCT || left->type == S_UFUNCT || 
	  ( left->type == S_UIDENT && left->val >= 0. ) ||
	  ( priority < symbol[left->token].info && !isalpha(node->symbol[0]) ) ||
	  ( priority == symbol[left->token].info && (type == S_ADD_OP ) ) )
	parenthesis = 0;
      if (parenthesis) _unur_fstr_print_C( out, "(", 0. );

      if (left->type == S_UCONST && left->val == 0. && node->symbol[0] == '-')
	/* there is no need to print "0 - ..." */ ;
      else
	_unur_fstr_node2C(out,left,variable,function);

      if (parenthesis) _unur_fstr_print_C( out, ")", 0. );
    }

    /* symbol for node */
    if (operator) _unur_fstr_print_C( out, " ", 0. );
    _unur_fstr_print_C( out, symb, node->val );
    if (operator) _unur_fstr_print_C( out, " ", 0. );

    /* right branch */
    if (right) {
      /* always use parenthesis ... */
      parenthesis = 1;
      /* ... except ... */
      if (right->type == S_SCONST || right->type == S_UCONST ||
	  right->type == S_SFUNCT || right->type == S_UFUNCT || 
	  ( right->type == S_UIDENT && right->val >= 0. ) ||
	  ( priority < symbol[right->token].info && !isalpha(node->symbol[0]) ) )
	parenthesis = 0;
      if (parenthesis) _unur_fstr_print_C( out, "(", 0. );
      _unur_fstr_node2C(out,right,variable,function);
      if (parenthesis) _unur_fstr_print_C( out, ")", 0. );
    }
  }

  return 1;
} /* end of _unur_fstr_node2C() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_print_C ( FILE *out, const char *symb, const double number )
     /*----------------------------------------------------------------------*/
     /* Print string or number into output stream.                           */
     /* The number is only printed if symb is the NULL pointer.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out    ... output stream                                           */
     /*   symb   ... string to be printed                                    */
     /*   number ... constant to be printed                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  if (symb)
    /* print symbol into output stream */
    fprintf(out,"%s",symb);
  else
    /* print number into output stream */
    fprintf(out,"%.16g",number);

  return 1;
} /* end of _unur_fstr_print_C() */

/*---------------------------------------------------------------------------*/
