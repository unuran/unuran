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
/* Codes for special C functions (must be added to code)                     */
/* (use bits in integer)                                                     */

enum {
  C_FUNCT_ERROR = 0x10000000u,      /* error                                 */

  C_FUNCT_SGN   = 0x00000001u,      /* sign function                         */
  C_FUNCT_SEC   = 0x00000002u       /* secant function                       */
};

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** API                                                                     **/
/*****************************************************************************/

/*--------------------------------------------------------------------------*/

int 
_unur_fstr_tree2C ( FILE *out, const struct ftreenode *root,
		    const char *variable, const char *funct_name )
     /*----------------------------------------------------------------------*/
     /* Produce string from function tree.                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out        ... output stream                                       */
     /*   root       ... pointer to root of function tree                    */
     /*   variable   ... pointer to name of variable                         */
     /*   funct_name ... pointer to name of C function                       */
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
  _unur_check_NULL( GENTYPE,symbol[root->token].node2C,0 );

  /* make body of C routine */
  rcode = symbol[root->token].node2C (&output,root,variable);
  if (rcode & C_FUNCT_ERROR) { 
    if (output.string) free(output.string);
    return 0;
  }

  /* print code for special functions (if necessary) */ 
  _unur_fstr_C_specfunct (out,rcode);

  /* print C routine */
  *(output.string + output.length) = '\0';
  fprintf (out,"static double %s (double %s)\n",funct_name,variable );
  fprintf (out,"{\n\treturn (%s);\n}\n",output.string);

  /* free memory */
  free(output.string);

  return 1;

} /* end of _unur_fstr_tree2C() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Create function string and source code.                                 **/
/*****************************************************************************/

/* C_xxxx ( struct concat *output, const struct ftreenode *node, const char *variable ) */
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
C_error ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* Error (This should not happen).                                      */
     /*----------------------------------------------------------------------*/
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return C_FUNCT_ERROR;
} /* end of C_error() */

/*---------------------------------------------------------------------------*/

unsigned
C_const ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for constant                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_C( output, NULL, node->val );
  return 0u;
} /* end of C_const() */

/*---------------------------------------------------------------------------*/

unsigned
C_var ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for variable                                                  */
     /*----------------------------------------------------------------------*/
{
  _unur_fstr_print_C( output, variable, 0 );
  return 0u;
} /* end of C_var() */

/*---------------------------------------------------------------------------*/

unsigned
C_prefix_generic ( struct concat *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
     /*----------------------------------------------------------------------*/
     /* print prefix operator (function). generic version                    */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  /* node '(' left ',' right ')' */
  _unur_fstr_print_C( output, symb, 0 );
  _unur_fstr_print_C( output, "(", 0 );

  if (left) {
    rcode |= symbol[left->token].node2C (output,left,variable);
    _unur_fstr_print_C( output, ",", 0 );
  }
  if (right) {
    rcode |= symbol[right->token].node2C (output,right,variable);
  }

  _unur_fstr_print_C( output, ")", 0 );

  return rcode;
} /* end of C_prefix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
C_prefix ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for prefix operator (function)                                */
     /*----------------------------------------------------------------------*/
{
  return C_prefix_generic(output,node,variable,symbol[node->token].name);
} /* end of C_prefix() */

/*---------------------------------------------------------------------------*/

unsigned
C_power ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for power function                                            */
     /*----------------------------------------------------------------------*/
{
  return C_prefix_generic(output,node,variable,"pow");
} /* end of C_power() */

/*---------------------------------------------------------------------------*/

unsigned
C_sec ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for secant function                                           */
     /*----------------------------------------------------------------------*/
{
  return (C_FUNCT_SEC | C_prefix_generic(output,node,variable,"_acg_sec"));
} /* end of C_sec() */

/*---------------------------------------------------------------------------*/

unsigned
C_abs ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for absolute value function                                   */
     /*----------------------------------------------------------------------*/
{
  return C_prefix_generic(output,node,variable,"fabs");
} /* end of C_power() */

/*---------------------------------------------------------------------------*/

unsigned
C_sgn ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for sign function                                            */
     /*----------------------------------------------------------------------*/
{
  return (C_FUNCT_SGN | C_prefix_generic(output,node,variable,"_acg_sgn"));
} /* end of C_sgn() */

/*---------------------------------------------------------------------------*/

unsigned
C_infix_generic ( struct concat *output, const struct ftreenode *node,
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
    return C_FUNCT_ERROR;

  /* '(' left node right ')' */
  _unur_fstr_print_C( output, "(", 0 );
  rcode |= symbol[left->token].node2C (output,left,variable);
  _unur_fstr_print_C( output, symb, 0 );
  rcode |= symbol[right->token].node2C (output,right,variable);
  _unur_fstr_print_C( output, ")", 0 );

  return rcode;

} /* end of C_infix_generic() */

/*---------------------------------------------------------------------------*/

unsigned
C_infix ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for infix operator (binary operator).                         */
     /*----------------------------------------------------------------------*/
{
  return C_infix_generic(output,node,variable,symbol[node->token].name);
} /* end of C_infix() */

/*---------------------------------------------------------------------------*/

unsigned
C_equal ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for equality operator                                         */
     /*----------------------------------------------------------------------*/
{
  return C_infix_generic(output,node,variable,"==");
} /* end of C_equal() */

/*---------------------------------------------------------------------------*/

unsigned
C_unequal ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for inequality operator                                       */
     /*----------------------------------------------------------------------*/
{
  return C_infix_generic(output,node,variable,"!=");
} /* end of C_unequal() */

/*---------------------------------------------------------------------------*/

unsigned
C_minus ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for minus operator                                            */
     /*----------------------------------------------------------------------*/
{
  unsigned rcode = 0u;

  struct ftreenode *left  = node->left;    /* left branch of node            */
  struct ftreenode *right = node->right;   /* right branch of node           */

  if (left==NULL || right==NULL)
    /* error */
    return C_FUNCT_ERROR;

  /* '(' [left] '-' right ')' */
  _unur_fstr_print_C( output, "(", 0 );
  if (!(left->type == S_UCONST && left->val == 0.))
    /* there is no need to print "0 - ..." */
    rcode |= symbol[left->token].node2C (output,left,variable);
  _unur_fstr_print_C( output, "-", 0 );
  rcode |= symbol[right->token].node2C (output,right,variable);
  _unur_fstr_print_C( output, ")", 0 );

  return rcode;
} /* end of C_minus() */

/*---------------------------------------------------------------------------*/

unsigned
C_mod ( struct concat *output, const struct ftreenode *node, const char *variable )
     /*----------------------------------------------------------------------*/
     /* string for mudolus function                                          */
     /*----------------------------------------------------------------------*/
{
  return C_infix_generic(output,node,variable,"%");
} /* end of C_mod() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Special functions                                                       **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_C_specfunct ( FILE *out, unsigned flags )
     /*----------------------------------------------------------------------*/
     /* Print C code for special functions                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out   ... output stream                                            */
     /*   flags ... bit array with flags for special functions               */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  if (flags & C_FUNCT_SGN) {
    _unur_fstr_C_sgn(out);
  }
  if (flags & C_FUNCT_SEC) {
    _unur_fstr_C_sec(out);
  }

  return 1;
} /* end of _unur_fstr_C_specfunct() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_C_sgn ( FILE *out )
     /*----------------------------------------------------------------------*/
     /* Print C code for sign function                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out   ... output stream                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"#ifndef _ACG_FUNCT_SGN\n");
  fprintf(out,"#define _ACG_FUNCT_SGN\n");

  fprintf(out,"static double _acg_sgn(double x)\n{\n");
  fprintf(out,"\treturn ((x<0.) ? -1. : ((x>0.) ? 1. : 0.));\n");
  fprintf(out,"}\n");

  fprintf(out,"#endif /* _ACG_FUNCT_SGN */\n\n");

  return 1;
} /* end of _unur_fstr_C_sgn() */

/*---------------------------------------------------------------------------*/

int
_unur_fstr_C_sec ( FILE *out )
     /*----------------------------------------------------------------------*/
     /* Print C code for secant function                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   out   ... output stream                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  fprintf(out,"#ifndef _ACG_FUNCT_SEC\n");
  fprintf(out,"#define _ACG_FUNCT_SEC\n");

  fprintf(out,"static double _acg_sec(double x)\n{\n");
  fprintf(out,"\tdouble cosx = cos(x);\n");
  fprintf(out,"\treturn ((cosx == 0.) ? HUGE_VAL : 1./cosx) ;\n");
  fprintf(out,"}\n");

  fprintf(out,"#endif /* _ACG_FUNCT_SEC */\n\n");

  return 1;
} /* end of _unur_fstr_C_sec() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Auxilliary routines                                                     **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_fstr_print_C ( struct concat *output, const char *symb, const double number )
     /*----------------------------------------------------------------------*/
     /* Print string or number as C code into output string.                 */
     /* The number is only printed if symb is the NULL pointer.              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   output ... pointer to string for output                            */
     /*   symb   ... string to be printed                                    */
     /*   number ... constant to be printed                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   1 on success                                                       */
     /*----------------------------------------------------------------------*/
{
  size_t len;

  /* (possible) length of output string */
  len = (symb) ? strlen(symb) : 64;
  
  /* Resize the allocated memory if necessary */
  if (output->length + len + 1 > output->allocated) {
    if (output->string == NULL) {
      output->allocated = 100;
      output->string = _unur_malloc( 100*sizeof(char) );
    }
    else {
      output->allocated = (output->allocated + len) * 2;
      output->string = _unur_realloc( output->string, output->allocated );
    }
  }

  if (symb)
    /* copy symbol into output */
    memcpy( output->string+output->length, symb, len );
  else
    /* copy number symbol into output */
    len = sprintf(output->string+output->length,"%.20e",number);

  /* update length of output string */
  output->length += len;

  return 1;
} /* end of _unur_fstr_print_C() */

/*---------------------------------------------------------------------------*/
