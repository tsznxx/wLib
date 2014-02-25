#include "common.h"
/* prototypes from inverse.c */
CPPEXTERN char *symbolset;    /* alphabet default is "AUGC" */
CPPEXTERN float inverse_fold(char *start, const char *target);  
/* find sequences with predefined structure.
   the found sequence is written to start,
   return value is
      energy_of_struct(start, target) - fold(start, structure),
   i.e. 0. if search was successful; */
   
CPPEXTERN float inverse_pf_fold(char *start, const char *target);
/*  inverse folding maximising the frequency of target in the
    ensemble of structures, final sequence is written to start, returns 
       energy_of_struct(start, target) - part_func(start, structure)
*/
CPPEXTERN float final_cost; /* when to stop inverse_pf_fold() */
CPPEXTERN int give_up;   /* default 0: try to minimize structure distance even if 
			 no exact solution can be found */

