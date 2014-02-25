#include "common.h"
#ifndef DIST_VARS_H
#include "dist_vars.h"  /* defines the type Tree */
#endif
CPPEXTERN  Tree   *make_tree(char *struc); /* make input for tree_edit_distance */
CPPEXTERN  float   tree_edit_distance(Tree *T1, Tree *T2);
/* compare to structures using tree editing */
CPPEXTERN  void    print_tree(Tree *t);    /* mainly for debugging */
CPPEXTERN  void    free_tree(Tree *t);     /* free space allocated by make_tree */
