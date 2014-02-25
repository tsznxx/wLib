#include "common.h"
#ifndef DUPLEX_H
#define DUPLEX_H

typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;

CPPEXTERN duplexT duplexfold(const char *s1, const char *s2);
CPPEXTERN duplexT *duplex_subopt(const char *s1, const char *s2, int delta, int w);
CPPEXTERN duplexT aliduplexfold(const char *s1[], const char *s2[]);
CPPEXTERN duplexT *aliduplex_subopt(const char *s1[], const char *s2[], int delta, int w);

#endif
