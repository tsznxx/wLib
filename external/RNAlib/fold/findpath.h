#include "common.h"
#ifndef FIND_PATH_H
#define FIND_PATH_H
typedef struct path {
  double en;
  char *s;
} path_t;

CPPEXTERN int find_saddle (char *seq, char *struc1, char *struc2, int max);
CPPEXTERN path_t* get_path(char *seq, char *s1, char* s2, int maxkeep);

#endif
