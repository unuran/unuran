/************************************************************************
 * MODUL pars.h                                                        *
 *                                                                      *
 *                                                                      * 
 ************************************************************************/

#include "scanpars.h"

/************************************************************************/
extern  void             _unur_fstr_init(void);                                          /* initialisiert hash table */
extern  struct treenode *_unur_fstr2tree(char *function, int *errcodep, int *errposp);  /* berrechnet aus string einen tree */
extern  double            _unur_fstr_eval_tree(struct treenode *E_root,double argument);     /* berechnet Funktionswert einer Funktion gegeben als tree */
extern  void             _unur_fstr_free(struct treenode *root);                          /* gibt Speicher fuer tree wieder frei */
extern  char            *Ntree2string(struct treenode *tree_root, char *ret_str);    /* berechnet aus tree string */
extern  struct treenode *_unur_fstr_make_derivative(struct treenode *root);                          /* berechnet 'Ableitungsbaum'*/

extern  double _unur_fstr_dev_eval_tree(struct treenode *E_root, double argument);
