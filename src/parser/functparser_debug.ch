/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_debug.c                                          *
 *                                                                           *
 *   Debugging tools for function parser.                                    *
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

  /* check arguments */
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);

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

  /* check arguments */
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
  CHECK_NULL(root,RETURN_VOID);   COOKIE_CHECK(root,CK_FSTR_TNODE,RETURN_VOID);

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

  /* check arguments */
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);

  fprintf(log,"%s: ",GENTYPE); 

  /* draw vertical lines in tree */
  for (i = 0, mask = 1; i < level; i++, mask <<= 1) 
    if (mask & location) 
      fprintf(log,"|   "); 
    else 
      fprintf(log,"    "); 
  
  /* print node */
  if( node != NULL ) {
    COOKIE_CHECK(node,CK_FSTR_TNODE,RETURN_VOID);
    
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

void
_unur_fstr_debug_deriv (const struct ftreenode *funct, const struct ftreenode *deriv)
     /*----------------------------------------------------------------------*/
     /* Print function and its derivative                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   funct ... pointer to function tree                                 */
     /*   deriv ... pointer to derivative                                    */
     /*----------------------------------------------------------------------*/
{
  FILE *log = unur_get_stream();
  char *str;

  /* check arguments */
  CHECK_NULL(funct,RETURN_VOID);  COOKIE_CHECK(funct,CK_FSTR_TNODE,RETURN_VOID);
  CHECK_NULL(deriv,RETURN_VOID);  COOKIE_CHECK(deriv,CK_FSTR_TNODE,RETURN_VOID);

  fprintf(log,"%s: Derivative df/dx of \n",GENTYPE);
  str = _unur_fstr_tree2string(funct,"x","f",TRUE);
  fprintf(log,"%s:  f(x) = %s\n",GENTYPE,str);
  free (str);

  if (deriv) {
    str = _unur_fstr_tree2string(deriv,"x","df",TRUE);
    fprintf(log,"%s:  f'(x) = %s\n",GENTYPE,str);
    free (str);
  }
  else {
    fprintf(log,"%s:  f'(x) = (unknown)\n",GENTYPE);
  }

  /*    _unur_fstr_debug_tree(NULL,deriv); */

  fprintf(log,"%s:\n",GENTYPE);

} /* end of _unur_fstr_debug_deriv() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
