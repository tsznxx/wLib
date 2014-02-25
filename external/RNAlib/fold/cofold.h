#include "common.h"
/* function from fold.c */
#include "subopt.h"
CPPEXTERN float  cofold(char *sequence, char *structure); 
/* calculate energy of string on structure */
CPPEXTERN void   free_co_arrays(void);          /* free arrays for mfe folding */
CPPEXTERN void   initialize_cofold(int length); /* allocate arrays for folding */
CPPEXTERN void   update_cofold_params(void);    /* recalculate parameters */
CPPEXTERN SOLUTION *zukersubopt(const char *string);
