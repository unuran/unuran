/************************************************************************
 * MODUL pars.h                                                        *
 *                                                                      *
 *                                                                      * 
 ************************************************************************/

#include "scanpars.h"



/************************************************************************/
extern  void             _unur_fstr_init(void);                                          /* initialisiert hash table */
extern  struct treenode *_unur_fstr2tree(char *function, int *errcodep, int *errposp);  /* berrechnet aus string einen tree */
extern  float            _unur_fstr_eval_tree(struct treenode *E_root,float argument);     /* berechnet Funktionswert einer Funktion gegeben als tree */
extern  void             _unur_fstr_free(struct treenode *root);                          /* gibt Speicher fuer tree wieder frei */
extern  char            *Ntree2string(struct treenode *tree_root, char *ret_str);    /* berechnet aus tree string *)
extern  struct treenode *_unur_fstr_make_derivative(struct treenode *root);                          /* berechnet 'Ableitungsbaum'*/

#if 0
/************************************************************************/
extern  void           void);                                          /* initialisiert hash table */
extern  struct treenode (char *function, int *errcodep, int *errposp);  /* berrechnet aus string einen tree */
extern  float          (struct treenode *E_root,float argument);     /* berechnet Funktionswert einer Funktion gegeben als tree */
extern  void            (struct treenode *root);                          /* gibt Speicher fuer tree wieder frei */
extern  char            *_unur_fstr_get_C_code(struct treenode *tree_root, char *ret_str);    /* berechnet aus tree string *)
extern  struct treenode *_unur_fstr_make_derivative(struct treenode *root);                          /* berechnet 'Ableitungsbaum'*/
#endif
