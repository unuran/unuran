/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      functparser_eval.c                                           *
 *                                                                           *
 *   Evaluate function tree for given argument x.                            *
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

double
_unur_fstr_eval_tree (const struct ftreenode *root, const double x)
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
  /* check arguments */
  CHECK_NULL(root,0.);   COOKIE_CHECK(root,CK_FSTR_TNODE,0.);

  return _unur_fstr_eval_node( root, x );
} /* end of _unur_fstr_eval_tree() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Routines for evaluating nodes of the function tree                      **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#define CHECK_INF(x)    if(_unur_FP_is_infinity((x))) return INFINITY;
#define CHECK_INFS(l,r) do { CHECK_INF((l)); CHECK_INF((r)); } while(0)
/*---------------------------------------------------------------------------*/

double v_dummy  (double l, double r) { return 0.; }
double v_const  (double l, double r) { return 0.; }  /* nothing to do, value in node */

double v_less   (double l, double r) { CHECK_INFS(l,r); return (double)(l <  r); }
double v_equal  (double l, double r) { CHECK_INFS(l,r); return (double)(l == r); }
double v_greater(double l, double r) { CHECK_INFS(l,r); return (double)(l >  r); }
double v_less_or(double l, double r) { CHECK_INFS(l,r); return (double)(l <= r); }
double v_unequal(double l, double r) { CHECK_INFS(l,r); return (double)(l != r); }
double v_grtr_or(double l, double r) { CHECK_INFS(l,r); return (double)(l >= r); }

double v_plus   (double l, double r) { CHECK_INFS(l,r); return (l + r); }
double v_minus  (double l, double r) { CHECK_INFS(l,r); return (l - r); }
double v_mul    (double l, double r) { CHECK_INFS(l,r); return (l * r); }
double v_div    (double l, double r) { CHECK_INFS(l,r); if (r==0.) return INFINITY;
                                                        return (l / r); }
double v_power  (double l, double r) { CHECK_INFS(l,r); return pow(l,r); }

double v_mod    (double l, double r) { CHECK_INFS(l,r); return (double)((int)l % (int)r); }
double v_exp    (double l, double r) { CHECK_INF(r);    return exp(r); }
double v_log    (double l, double r) { CHECK_INF(r);    if(r<=0.) return INFINITY; return log(r); }
double v_sin    (double l, double r) { CHECK_INF(r);    return sin(r); }
double v_cos    (double l, double r) { CHECK_INF(r);    return cos(r); }
double v_tan    (double l, double r) { CHECK_INF(r);    return tan(r); }
double v_sec    (double l, double r) { double cosr; CHECK_INF(r); if((cosr=cos(r))==0.) return INFINITY;
                                                        return 1./cosr; }
double v_sqrt   (double l, double r) { CHECK_INF(r);    if(r<0.) return INFINITY;
                                                        return sqrt(r); }
double v_abs    (double l, double r) { CHECK_INF(r);    return fabs(r); }
double v_sgn    (double l, double r) { CHECK_INF(r);    return ((r<0.) ? -1. : ((r>0.) ? 1. : 0.)); }

/*---------------------------------------------------------------------------*/
#undef CHECK_INF
#undef CHECK_INFS
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/** Evaluate function                                                       **/
/*****************************************************************************/

double
_unur_fstr_eval_node (const struct ftreenode *node, const double x)
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

  /* check arguments */
  CHECK_NULL(node,0.);   COOKIE_CHECK(node,CK_FSTR_TNODE,0.);

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
