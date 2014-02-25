#include "common.h"
#ifndef ALN_UTIL_H
#define ALN_UTIL_H

CPPEXTERN int read_clustal(FILE *clust, char *AlignedSeqs[], char *names[]);
CPPEXTERN /*@only@*/ /*@notnull@*/ char *consensus(const char *AS[]);
CPPEXTERN /*@only@*/ /*@notnull@*/ char *consens_mis(const char *AS[]);

#endif
